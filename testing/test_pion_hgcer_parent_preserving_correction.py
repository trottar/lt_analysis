"""Focused contracts for the detached E.7.1 parent-preserving map."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))

import pion_hgcer_parent_preserving_correction as correction


def _fingerprint(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")).hexdigest()


def _e7_cell(t_edges, delta_edges, t_index, delta_index, scale=None):
    available = scale is not None
    return {
        "t_index": t_index,
        "t_low": t_edges[t_index],
        "t_high": t_edges[t_index + 1],
        "delta_index": delta_index,
        "delta_low": delta_edges[delta_index],
        "delta_high": delta_edges[delta_index + 1],
        "method_a": {
            "present": True,
            "candidate": 4.0 if available else 1.0,
            "low": 3.0 if available else 0.8,
            "high": 5.0 if available else 1.2,
            "status": "marginal" if available else "available",
        },
        "method_b": {
            "present": available,
            "candidate": 9.0 if available else None,
            "uncertainty": 1.0 if available else None,
            "status": "available_multi_region" if available else "unavailable",
        },
        "comparison": {
            "availability": "both_comparable" if available else "a_only",
            "ratio_B_over_A": 2.25 if available else None,
            "log_ratio_B_over_A": math.log(2.25) if available else None,
            "diagnostic_interval_relation": "disjoint" if available else "not_evaluable",
        },
        "prototype_status": "available" if available else "unavailable",
        "prototype_reason": None if available else "method_b_unavailable",
        "prototype_log_scale": math.log(scale) if available else None,
        "prototype_relative_scale": scale,
    }


def _finalize_prototype(payload):
    payload["fingerprint_inputs"] = {
        "schema_version": "pion_hgcer_ab_combination_prototype/v1",
        "method": "equal_weight_log_geometric_mean",
        "log_weight_method_a": 0.5,
        "log_weight_method_b": 0.5,
        "source_phase_d_checkpoint_schema": payload["source_phase_d_checkpoint_schema"],
        "source_checkpoint_payload_fingerprint": payload[
            "source_checkpoint_payload_fingerprint"
        ],
        "phase_a_contract_fingerprint": payload["phase_a_contract_fingerprint"],
        "coordinate_fingerprint": payload["coordinate_fingerprint"],
        "method_a_comparison_fingerprint": payload[
            "method_a_comparison_fingerprint"
        ],
        "method_b_comparison_fingerprint": payload[
            "method_b_comparison_fingerprint"
        ],
        "source_method_a_comparison_payload_fingerprint": payload[
            "source_method_a_comparison_payload_fingerprint"
        ],
        "source_method_b_comparison_payload_fingerprint": payload[
            "source_method_b_comparison_payload_fingerprint"
        ],
        "ab_comparison_fingerprint": payload["ab_comparison_fingerprint"],
        "host_state": payload["host_state"],
        "source_target_state": payload["source_target_state"],
        "t_edges": payload["t_edges"],
        "delta_edges": payload["delta_edges"],
        "cells": payload["cells"],
    }
    payload["fingerprint"] = _fingerprint(payload["fingerprint_inputs"])
    return payload


def _prototype(t_edges=None, delta_edges=None, scales=None):
    t_edges = list(t_edges or (0.0, 1.0, 2.0))
    delta_edges = list(delta_edges or (-2.0, 0.0, 2.0))
    scales = list(scales or (2.0, 0.5, 3.0, None))
    cells = []
    index = 0
    for t_index in range(len(t_edges) - 1):
        for delta_index in range(len(delta_edges) - 1):
            cells.append(
                _e7_cell(t_edges, delta_edges, t_index, delta_index, scales[index])
            )
            index += 1
    return _finalize_prototype({
        "schema_version": "pion_hgcer_ab_combination_prototype/v1",
        "status": "available",
        "available": True,
        "reason": None,
        "method": "equal_weight_log_geometric_mean",
        "log_weight_method_a": 0.5,
        "log_weight_method_b": 0.5,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "production_application_performed": False,
        "prototype_statistical_optimality_claimed": False,
        "prototype_uncertainty_model": "not_defined",
        "source_phase_d_checkpoint_schema": "pion_hgcer_phase_d_checkpoint/v1",
        "source_checkpoint_payload_fingerprint": "phase-c-source",
        "phase_a_contract_fingerprint": "phase-a-contract",
        "coordinate_fingerprint": "coordinate-contract",
        "method_a_comparison_fingerprint": "method-a-representation",
        "method_b_comparison_fingerprint": "method-b-representation",
        "source_method_a_comparison_payload_fingerprint": "method-a-source",
        "source_method_b_comparison_payload_fingerprint": "method-b-source",
        "ab_comparison_fingerprint": "ab-representation",
        "host_state": "proton_cleaned",
        "source_target_state": "post_proton_noRF",
        "t_edges": t_edges,
        "delta_edges": delta_edges,
        "cells": cells,
    })


def _record(t_edges, delta_edges, t_index, delta_index, contribution, *, nommcuts=True):
    record = {
        "nommcuts": nommcuts,
        "canonical_t_index": t_index,
        "canonical_t_lower_edge": t_edges[t_index],
        "canonical_t_upper_edge": t_edges[t_index + 1],
        "delta_index": delta_index,
        "delta_lower_edge": None if delta_index is None else delta_edges[delta_index],
        "delta_upper_edge": None if delta_index is None else delta_edges[delta_index + 1],
        "signed_source_coefficient": contribution,
        "baseline_pion_weight_w0": 1.0,
        "signed_baseline_event_contribution": contribution,
        "noRF_provenance": "noRF",
    }
    return record


def _phase_a(prototype, records=None):
    t_edges = prototype["t_edges"]
    delta_edges = prototype["delta_edges"]
    if records is None:
        records = [
            _record(t_edges, delta_edges, 0, 0, 2.0),
            _record(t_edges, delta_edges, 0, 1, 6.0),
            _record(t_edges, delta_edges, 0, None, 3.0),
            _record(t_edges, delta_edges, 0, 0, 1000.0, nommcuts=False),
            _record(t_edges, delta_edges, 1, 0, 4.0),
            _record(t_edges, delta_edges, 1, 1, 5.0),
        ]
    payload = {
        "schema_version": "pion_hgcer_event_contract/v1",
        "fingerprint_schema_version": "pion_hgcer_event_contract_fingerprint/v2",
        "status": "available",
        "available": True,
        "immutable_record_contract": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "rf_restoration_applied": False,
        "pion_closure": {"passed": True},
        "host_closure": {"passed": True},
        "coordinate_fingerprint": prototype["coordinate_fingerprint"],
        "pion_event_population_fingerprint": "phase-a-pion-population",
        "canonical_t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "pion_records": records,
        "host_state": prototype["host_state"],
        "source_target_state": prototype["source_target_state"],
        "fingerprint_inputs": {"contract": "phase-a-fixture"},
    }
    payload["contract_fingerprint"] = _fingerprint(payload["fingerprint_inputs"])
    prototype["phase_a_contract_fingerprint"] = payload["contract_fingerprint"]
    _finalize_prototype(prototype)
    return payload


def _valid_inputs():
    prototype = _prototype()
    phase_a = _phase_a(prototype)
    return prototype, phase_a


def _cell(payload, t_index, delta_index):
    return next(
        row for row in payload["cells"]
        if row["t_index"] == t_index and row["delta_index"] == delta_index
    )


class ParentPreservingCorrectionTests(unittest.TestCase):
    def test_multi_cell_signed_normalization_closure_and_outside_delta_audit(self):
        prototype, phase_a = _valid_inputs()
        before_prototype = json.dumps(prototype, sort_keys=True, allow_nan=False)
        before_phase_a = json.dumps(phase_a, sort_keys=True, allow_nan=False)

        result = correction.build_pion_hgcer_parent_preserving_correction(
            prototype, phase_a
        )

        self.assertTrue(result["available"])
        self.assertEqual(
            result["schema_version"],
            correction.PION_HGCER_PARENT_PRESERVING_CORRECTION_SCHEMA_VERSION,
        )
        self.assertTrue(result["non_authoritative"])
        self.assertFalse(result["production_objects_mutated"])
        self.assertFalse(result["event_application_performed"])
        self.assertEqual(len(result["cells"]), 4)
        self.assertEqual(len(result["parents"]), 2)
        first = result["parents"][0]
        self.assertEqual(first["parent_status"], "available_parent_preserved")
        self.assertAlmostEqual(first["normalization_divisor"], 7.0 / 8.0)
        self.assertAlmostEqual(_cell(result, 0, 0)["C_final"], 16.0 / 7.0)
        self.assertAlmostEqual(_cell(result, 0, 1)["C_final"], 4.0 / 7.0)
        self.assertAlmostEqual(first["canonical_baseline_before"], 8.0)
        self.assertAlmostEqual(first["canonical_baseline_after"], 8.0)
        self.assertAlmostEqual(first["outside_delta_baseline"], 3.0)
        self.assertAlmostEqual(first["full_parent_baseline_before"], 11.0)
        self.assertAlmostEqual(first["full_parent_baseline_after"], 11.0)
        self.assertTrue(first["closure_passed"])
        self.assertEqual(_cell(result, 0, 0)["baseline_record_count"], 1)
        self.assertEqual(_cell(result, 0, 0)["baseline_signed_sum"], 2.0)
        self.assertEqual(_cell(result, 0, 0)["baseline_absolute_support"], 2.0)
        self.assertEqual(_cell(result, 0, 0)["baseline_sumw2"], 4.0)
        self.assertEqual(_cell(result, 0, 0)["baseline_neff"], 1.0)
        self.assertEqual(result["parents"][1]["parent_status"], "identity_single_refinable_cell")
        self.assertEqual(_cell(result, 1, 0)["e7_raw_relative_scale"], 3.0)
        self.assertEqual(_cell(result, 1, 0)["C_final"], 1.0)
        self.assertEqual(
            _cell(result, 1, 0)["C_final_status"], "identity_single_refinable_cell"
        )
        self.assertEqual(_cell(result, 1, 1)["C_final"], 1.0)
        self.assertEqual(_cell(result, 1, 1)["C_final_status"], "identity_unmodified")
        self.assertEqual(json.dumps(prototype, sort_keys=True, allow_nan=False), before_prototype)
        self.assertEqual(json.dumps(phase_a, sort_keys=True, allow_nan=False), before_phase_a)
        json.dumps(result, sort_keys=True, allow_nan=False)

    def test_nommcuts_false_does_not_change_signed_parent_normalization(self):
        prototype, phase_a = _valid_inputs()
        baseline = correction.build_pion_hgcer_parent_preserving_correction(
            prototype, phase_a
        )
        phase_a["pion_records"].append(
            _record(prototype["t_edges"], prototype["delta_edges"], 0, 1, -99999.0, nommcuts=False)
        )
        modified = correction.build_pion_hgcer_parent_preserving_correction(
            prototype, phase_a
        )
        self.assertTrue(modified["available"])
        self.assertEqual(
            baseline["parents"][0]["normalization_divisor"],
            modified["parents"][0]["normalization_divisor"],
        )
        self.assertEqual(_cell(baseline, 0, 1)["C_final"], _cell(modified, 0, 1)["C_final"])
        self.assertEqual(_cell(modified, 0, 1)["baseline_signed_sum"], 6.0)

    def test_negative_baseline_cell_uses_signed_algebra_not_absolute_support(self):
        prototype, phase_a = _valid_inputs()
        prototype["cells"][0]["prototype_relative_scale"] = 0.5
        prototype["cells"][0]["prototype_log_scale"] = math.log(0.5)
        prototype["cells"][1]["prototype_relative_scale"] = 2.0
        prototype["cells"][1]["prototype_log_scale"] = math.log(2.0)
        _finalize_prototype(prototype)
        phase_a["pion_records"] = [
            _record(prototype["t_edges"], prototype["delta_edges"], 0, 0, -2.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 0, 1, 3.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 1, 0, 4.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 1, 1, 5.0),
        ]
        result = correction.build_pion_hgcer_parent_preserving_correction(prototype, phase_a)
        parent = result["parents"][0]
        self.assertTrue(result["available"])
        self.assertEqual(parent["negative_baseline_cell_count"], 1)
        self.assertAlmostEqual(parent["normalization_divisor"], 5.0)
        self.assertAlmostEqual(_cell(result, 0, 0)["C_final"], 0.5 / 5.0)
        self.assertAlmostEqual(_cell(result, 0, 1)["C_final"], 2.0 / 5.0)
        self.assertAlmostEqual(parent["canonical_baseline_after"], 1.0)

    def test_no_refinable_and_invalid_parent_are_local_identity_cases(self):
        prototype, phase_a = _valid_inputs()
        prototype["cells"][0]["prototype_status"] = "unavailable"
        prototype["cells"][0]["prototype_reason"] = "method_b_unavailable"
        prototype["cells"][0]["prototype_relative_scale"] = None
        prototype["cells"][0]["prototype_log_scale"] = None
        prototype["cells"][0]["comparison"]["availability"] = "a_only"
        prototype["cells"][0]["comparison"]["ratio_B_over_A"] = None
        prototype["cells"][0]["comparison"]["log_ratio_B_over_A"] = None
        prototype["cells"][0]["comparison"]["diagnostic_interval_relation"] = "not_evaluable"
        prototype["cells"][0]["method_b"] = {
            "present": False,
            "candidate": None,
            "uncertainty": None,
            "status": "unavailable",
        }
        _finalize_prototype(prototype)
        phase_a["pion_records"] = [
            _record(prototype["t_edges"], prototype["delta_edges"], 0, 0, 2.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 0, 1, 3.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 1, 0, -2.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 1, 1, 3.0),
        ]
        result = correction.build_pion_hgcer_parent_preserving_correction(prototype, phase_a)
        self.assertTrue(result["available"])
        self.assertEqual(result["parents"][0]["parent_status"], "identity_single_refinable_cell")
        self.assertEqual(
            result["parents"][1]["parent_status"], "identity_single_refinable_cell"
        )
        self.assertTrue(all(row["C_final"] == 1.0 for row in result["cells"]))

    def test_no_refinable_parent_preserves_every_cell_at_identity(self):
        prototype = _prototype(scales=(None, None, None, None))
        phase_a = _phase_a(prototype)
        result = correction.build_pion_hgcer_parent_preserving_correction(prototype, phase_a)
        self.assertTrue(result["available"])
        self.assertEqual(
            [parent["parent_status"] for parent in result["parents"]],
            ["identity_no_refinable_cells", "identity_no_refinable_cells"],
        )
        for cell in result["cells"]:
            self.assertEqual(cell["C_final"], 1.0)
            self.assertEqual(cell["C_final_status"], "identity_unmodified")
            self.assertEqual(cell["e7_raw_relative_scale"], None)

    def test_nonpositive_refinable_signed_sum_and_raw_weighted_sum_fail_closed(self):
        prototype, phase_a = _valid_inputs()
        phase_a["pion_records"] = [
            _record(prototype["t_edges"], prototype["delta_edges"], 0, 0, 2.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 0, 1, -2.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 1, 0, 4.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 1, 1, 5.0),
        ]
        result = correction.build_pion_hgcer_parent_preserving_correction(prototype, phase_a)
        self.assertEqual(
            result["parents"][0]["parent_reason"],
            "refinable_baseline_signed_sum_nonpositive",
        )
        self.assertTrue(all(_cell(result, 0, index)["C_final"] == 1.0 for index in (0, 1)))

        prototype, phase_a = _valid_inputs()
        prototype["cells"][0]["prototype_relative_scale"] = 10.0
        prototype["cells"][0]["prototype_log_scale"] = math.log(10.0)
        prototype["cells"][1]["prototype_relative_scale"] = 0.1
        prototype["cells"][1]["prototype_log_scale"] = math.log(0.1)
        _finalize_prototype(prototype)
        phase_a["pion_records"] = [
            _record(prototype["t_edges"], prototype["delta_edges"], 0, 0, -2.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 0, 1, 3.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 1, 0, 4.0),
            _record(prototype["t_edges"], prototype["delta_edges"], 1, 1, 5.0),
        ]
        result = correction.build_pion_hgcer_parent_preserving_correction(prototype, phase_a)
        self.assertEqual(
            result["parents"][0]["parent_reason"],
            "refinable_raw_scaled_signed_sum_nonpositive",
        )
        self.assertTrue(all(_cell(result, 0, index)["C_final"] == 1.0 for index in (0, 1)))

    def test_forced_closure_failure_retains_attempt_and_resets_parent(self):
        prototype, phase_a = _valid_inputs()
        with mock.patch.object(correction, "_closure_tolerance", return_value=-1.0):
            result = correction.build_pion_hgcer_parent_preserving_correction(
                prototype, phase_a
            )
        parent = result["parents"][0]
        self.assertEqual(parent["parent_reason"], "parent_preservation_closure_failed")
        self.assertIsNone(parent["normalization_divisor"])
        self.assertIsNotNone(parent["attempted_normalization"]["normalization_divisor"])
        self.assertFalse(parent["attempted_normalization"]["closure_passed"])
        self.assertTrue(all(_cell(result, 0, index)["C_final"] == 1.0 for index in (0, 1)))

    def test_contract_provenance_grid_and_record_failures_are_unavailable(self):
        prototype, phase_a = _valid_inputs()
        cases = []
        broken = deepcopy(prototype)
        broken["schema_version"] = "wrong"
        cases.append((broken, phase_a, "e7_prototype_contract_invalid"))
        broken = deepcopy(prototype)
        broken["cells"] = broken["cells"][:-1]
        _finalize_prototype(broken)
        cases.append((broken, phase_a, "e7_prototype_cell_grid_invalid"))
        broken = deepcopy(prototype)
        broken["cells"][0]["prototype_relative_scale"] = 0.0
        _finalize_prototype(broken)
        cases.append((broken, phase_a, "e7_prototype_cell_contract_invalid"))
        broken_phase_a = deepcopy(phase_a)
        broken_phase_a["coordinate_fingerprint"] = "wrong"
        cases.append((prototype, broken_phase_a, "e7_phase_a_provenance_mismatch"))
        broken_phase_a = deepcopy(phase_a)
        broken_phase_a["pion_records"][0]["canonical_t_upper_edge"] = 123.0
        cases.append((prototype, broken_phase_a, "phase_a_pion_record_contract_invalid:0"))
        broken_phase_a = deepcopy(phase_a)
        broken_phase_a["pion_records"][0]["delta_upper_edge"] = 123.0
        cases.append((prototype, broken_phase_a, "phase_a_pion_record_contract_invalid:0"))
        broken_phase_a = deepcopy(phase_a)
        broken_phase_a["pion_records"][2]["delta_lower_edge"] = 0.0
        cases.append((prototype, broken_phase_a, "phase_a_pion_record_contract_invalid:2"))
        for source_prototype, source_phase_a, reason in cases:
            with self.subTest(reason=reason):
                result = correction.build_pion_hgcer_parent_preserving_correction(
                    source_prototype, source_phase_a
                )
                self.assertFalse(result["available"])
                self.assertEqual(result["reason"], reason)

    def test_authority_geometry_and_contributing_record_contracts_fail_globally(self):
        prototype, phase_a = _valid_inputs()
        cases = []
        broken = deepcopy(prototype)
        broken["non_authoritative"] = False
        cases.append((broken, phase_a, "e7_prototype_contract_invalid"))
        broken_phase_a = deepcopy(phase_a)
        broken_phase_a["pion_closure"]["passed"] = False
        cases.append((prototype, broken_phase_a, "phase_a_contract_invalid"))
        broken_phase_a = deepcopy(phase_a)
        broken_phase_a["source_target_state"] = "wrong"
        cases.append((prototype, broken_phase_a, "phase_a_contract_invalid"))
        broken_phase_a = deepcopy(phase_a)
        broken_phase_a["canonical_t_edges"] = [0.0, 1.5, 2.0]
        cases.append((prototype, broken_phase_a, "e7_phase_a_provenance_mismatch"))
        for field, value in (
            ("canonical_t_index", -1),
            ("delta_index", -1),
            ("signed_baseline_event_contribution", float("inf")),
            ("signed_source_coefficient", float("inf")),
            ("baseline_pion_weight_w0", float("inf")),
            ("signed_source_coefficient", 3.0),
        ):
            broken_phase_a = deepcopy(phase_a)
            broken_phase_a["pion_records"][0][field] = value
            cases.append((prototype, broken_phase_a, "phase_a_pion_record_contract_invalid:0"))
        for source_prototype, source_phase_a, reason in cases:
            with self.subTest(reason=reason):
                result = correction.build_pion_hgcer_parent_preserving_correction(
                    source_prototype, source_phase_a
                )
                self.assertFalse(result["available"])
                self.assertEqual(result["reason"], reason)

    def test_fingerprint_determinism_and_source_change_sensitivity(self):
        prototype, phase_a = _valid_inputs()
        first = correction.build_pion_hgcer_parent_preserving_correction(prototype, phase_a)
        second = correction.build_pion_hgcer_parent_preserving_correction(
            deepcopy(prototype), deepcopy(phase_a)
        )
        self.assertEqual(first["fingerprint"], second["fingerprint"])
        changed = deepcopy(phase_a)
        changed["pion_records"][0]["signed_baseline_event_contribution"] = 3.0
        changed["pion_records"][0]["signed_source_coefficient"] = 3.0
        updated = correction.build_pion_hgcer_parent_preserving_correction(prototype, changed)
        self.assertTrue(updated["available"])
        self.assertNotEqual(first["fingerprint"], updated["fingerprint"])

    def test_artifact_filename_and_deterministic_round_trip(self):
        prototype, phase_a = _valid_inputs()
        result = correction.build_pion_hgcer_parent_preserving_correction(prototype, phase_a)
        setting = {
            "kinematic_token": "Q4p4W2p74",
            "Q2": 4.4,
            "W": 2.74,
            "epsilon_setting": "low",
            "epsilon_filename_token": "lowe",
            "phi_setting": "Left",
            "particle_type": "kaon",
        }
        artifact = correction.build_pion_hgcer_parent_preserving_correction_artifact(
            setting=setting, correction=result
        )
        basename = correction.pion_hgcer_parent_preserving_correction_filename(
            "Left", "kaon", "Q4p4W2p74", "lowe"
        )
        self.assertEqual(
            basename,
            "Left_kaon_pion-background_hgcer_parent-preserving-correction_Q4p4W2p74_lowe.json",
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / basename
            correction.write_pion_hgcer_parent_preserving_correction_json(path, artifact)
            raw = path.read_text(encoding="utf-8")
            self.assertTrue(raw.endswith("\n"))
            self.assertEqual(json.loads(raw), artifact)

    def test_module_remains_detached_from_production_and_reconstruction_helpers(self):
        source = (REPO_ROOT / "src" / "cuts" / "pion_hgcer_parent_preserving_correction.py").read_text(
            encoding="utf-8"
        )
        for forbidden in (
            "import ROOT", "find_canonical_bin", "build_pion_hgcer_method_a",
            "build_pion_hgcer_method_b", "build_pion_hgcer_ab_comparison",
            "simc_shape_pion_weight_from_value", "calculate_yield", "interp1d",
            "gaussian_filter", "numpy.clip", "np.clip", "math.sqrt",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
