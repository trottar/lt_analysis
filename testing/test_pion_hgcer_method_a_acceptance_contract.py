"""Focused F.1 detached Method-A acceptance-contract tests."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))

import pion_hgcer_method_a_acceptance_contract as acceptance


T_EDGES = [0.0, 1.0, 2.0]
DELTA_EDGES = [-10.0, 0.0, 10.0]
PHI_EDGES = [-20.0, 0.0, 20.0]


class _AmbiguousTruthArray:
    """Indexable test stand-in for NumPy's ambiguous array truth value."""

    def __init__(self, values):
        self._values = list(values)

    def __len__(self):
        return len(self._values)

    def __getitem__(self, index):
        return self._values[index]

    def __bool__(self):
        raise ValueError(
            "The truth value of an array with more than one element is ambiguous."
        )


class _ArrayScalar:
    """Test-only NumPy-like scalar that requires explicit .item() detachment."""

    def __init__(self, value):
        self._value = value

    def item(self):
        return self._value


class _BrokenArrayScalar:
    def item(self):
        raise RuntimeError("synthetic child scalar failure")


def _phase_record(source, entry, t_index, delta_index, delta, npe, *, nommcuts=True):
    return {
        "source_label": source,
        "entry_index": entry,
        "coordinate_fingerprint": "coordinates",
        "canonical_t_index": t_index,
        "canonical_t_lower_edge": T_EDGES[t_index],
        "canonical_t_upper_edge": T_EDGES[t_index + 1],
        "analysis_abs_t": T_EDGES[t_index] + 0.25,
        "analysis_MM": 1.1 + 0.01 * entry,
        "signed_source_coefficient": -0.25 if entry == 1 else 1.0,
        "allcuts": True,
        "nommcuts": nommcuts,
        "SHMS_delta": delta,
        "delta_index": delta_index,
        "delta_lower_edge": None if delta_index is None else DELTA_EDGES[delta_index],
        "delta_upper_edge": None if delta_index is None else DELTA_EDGES[delta_index + 1],
        "P_hgcer_npeSum": npe,
        "P_hgcer_xAtCer": 1.0 + entry,
        "P_hgcer_yAtCer": -1.0 - entry,
        "baseline_pion_weight_w0": 0.5 + 0.1 * entry,
        "signed_baseline_event_contribution": (-0.25 if entry == 1 else 1.0) * (0.5 + 0.1 * entry),
        "noRF_provenance": "noRF",
    }


def _parent(record, phi_index):
    entry = record["entry_index"]
    return {
        "source_label": record["source_label"],
        "entry_index": entry,
        "coordinate_fingerprint": record["coordinate_fingerprint"],
        "rf_state": "noRF",
        "source_tree_name": "PionEvents_noRF",
        "t_index": record["canonical_t_index"],
        "adj_t": record["analysis_abs_t"],
        "adj_MM": record["analysis_MM"],
        "coefficient": record["signed_source_coefficient"],
        "allcuts": record["allcuts"],
        "nommcuts": record["nommcuts"],
        "ssdelta": record["SHMS_delta"],
        "delta_index": record["delta_index"],
        "P_hgcer_npeSum": record["P_hgcer_npeSum"],
        "P_hgcer_xAtCer": record["P_hgcer_xAtCer"],
        "P_hgcer_yAtCer": record["P_hgcer_yAtCer"],
        "ssxptar": -0.04 + 0.01 * entry,
        "ssyptar": 0.02 - 0.005 * entry,
        "hsxptar": 0.001 * entry,
        "hsyptar": -0.002 * entry,
        "phi_degrees": -10.0 if phi_index == 0 else 10.0 if phi_index == 1 else 30.0,
        "phi_index": phi_index,
    }


def _child(parent):
    return {
        "source_label": parent["source_label"],
        "entry_index": parent["entry_index"],
        "coefficient": parent["coefficient"],
        "coordinate_fingerprint": parent["coordinate_fingerprint"],
        "adj_t": parent["adj_t"],
        "adj_MM": parent["adj_MM"],
        "ssdelta": parent["ssdelta"],
        "P_hgcer_npeSum": parent["P_hgcer_npeSum"],
        "P_hgcer_xAtCer": parent["P_hgcer_xAtCer"],
        "P_hgcer_yAtCer": parent["P_hgcer_yAtCer"],
        "ssxptar": parent["ssxptar"],
        "ssyptar": parent["ssyptar"],
        "hsxptar": parent["hsxptar"],
        "hsyptar": parent["hsyptar"],
        "phi_degrees": parent["phi_degrees"],
        "allcuts": parent["allcuts"],
        "nommcuts": parent["nommcuts"],
        "t_index": parent["t_index"],
        "phi_index": parent["phi_index"],
        "Q2": 4.4,
        "W": 2.74,
        "epsilon": 0.3,
        "theta_cm_deg": 12.0,
    }


def _fixture():
    records = [
        _phase_record("prompt", 0, 0, 0, -5.0, 1.5),
        _phase_record("prompt", 1, 0, 1, 5.0, 2.0),
        _phase_record("prompt", 4, 0, 0, -5.0, 1.2),
        _phase_record("prompt", 5, 0, 1, 5.0, 2.5),
        _phase_record("prompt", 2, 1, 0, -5.0, 2.1),
        _phase_record("rand", 3, 1, 1, 5.0, 3.0),
    ]
    parents = [
        _parent(record, phi_index)
        for record, phi_index in zip(records, (0, 1, 1, 1, 1, None))
    ]
    children = [_child(parent) for parent in parents if parent["phi_index"] is not None]
    fields = tuple(children[0])
    child_columns = {field: [row[field] for row in children] for field in fields}
    phase = {
        "schema_version": "pion_hgcer_event_contract/v1",
        "fingerprint_schema_version": "pion_hgcer_event_contract_fingerprint/v2",
        "status": "available", "available": True,
        "contract_fingerprint": "phase-contract",
        "pion_event_population_fingerprint": "phase-pion-events",
        "coordinate_fingerprint": "coordinates",
        "canonical_t_edges": list(T_EDGES), "delta_edges": list(DELTA_EDGES),
        "pion_records": records, "host_state": "proton_cleaned",
        "source_target_state": "post_proton_noRF", "rf_restoration_applied": False,
        "immutable_record_contract": True, "production_objects_mutated": False,
        "refinement_applied": False,
    }
    cache = {
        "coordinate_fingerprint": "coordinates", "delta_edges": list(DELTA_EDGES),
        "records": parents, "child_event_cache": {"prompt": child_columns},
    }
    return phase, cache


def _runtime_like_child_scalars(cache):
    for section in cache["child_event_cache"].values():
        for field, values in section.items():
            section[field] = [_ArrayScalar(value) for value in values]


class MethodAAcceptanceContractTests(unittest.TestCase):
    def test_runtime_like_child_scalars_preserve_parity_and_fingerprints(self):
        baseline_phase, baseline_cache = _fixture()
        baseline = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            baseline_phase, baseline_cache, phi_edges=PHI_EDGES
        )
        phase, cache = _fixture()
        _runtime_like_child_scalars(cache)
        runtime_like = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=PHI_EDGES
        )
        self.assertTrue(runtime_like["available"], runtime_like.get("reason"))
        self.assertEqual(runtime_like["status"], "available")
        self.assertEqual(runtime_like["diagnostic_stage"], "complete")
        for field in (
            "records", "summary", "feature_metadata",
            "event_population_fingerprint",
            "acceptance_feature_metadata_fingerprint",
            "child_assignment_projection_fingerprint", "fingerprint_inputs", "fingerprint",
        ):
            self.assertEqual(runtime_like[field], baseline[field])

    def test_runtime_like_boolean_parity_remains_exact(self):
        phase, cache = _fixture()
        _runtime_like_child_scalars(cache)
        self.assertTrue(acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=PHI_EDGES
        )["available"])
        for field, reason in (
            ("allcuts", "parent_child_parity_mismatch:allcuts"),
            ("nommcuts", "parent_child_parity_mismatch:nommcuts"),
        ):
            with self.subTest(field=field):
                candidate_phase, candidate_cache = _fixture()
                _runtime_like_child_scalars(candidate_cache)
                candidate_cache["child_event_cache"]["prompt"][field][0] = _ArrayScalar(False)
                result = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
                    candidate_phase, candidate_cache, phi_edges=PHI_EDGES
                )
                self.assertFalse(result["available"])
                self.assertEqual(result["reason"], reason)

    def test_malformed_runtime_like_child_scalar_fails_locally(self):
        phase, cache = _fixture()
        _runtime_like_child_scalars(cache)
        cache["child_event_cache"]["prompt"]["entry_index"][0] = _BrokenArrayScalar()
        result = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=PHI_EDGES
        )
        self.assertFalse(result["available"])
        self.assertEqual(result["status"], "unavailable")
        self.assertEqual(result["reason"], "pion_cache_child_scalar_invalid")

    def test_ambiguous_truth_edge_arrays_preserve_geometry_and_fingerprints(self):
        baseline_phase, baseline_cache = _fixture()
        baseline = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            baseline_phase, baseline_cache, phi_edges=PHI_EDGES
        )
        phase, cache = _fixture()
        phase["canonical_t_edges"] = _AmbiguousTruthArray(T_EDGES)
        phase["delta_edges"] = _AmbiguousTruthArray(DELTA_EDGES)
        cache["delta_edges"] = _AmbiguousTruthArray(DELTA_EDGES)
        runtime_like = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=_AmbiguousTruthArray(PHI_EDGES)
        )
        self.assertTrue(runtime_like["available"], runtime_like.get("reason"))
        self.assertEqual(runtime_like["status"], "available")
        self.assertEqual(runtime_like["diagnostic_stage"], "complete")
        for field in (
            "records", "summary", "feature_metadata", "t_edges", "delta_edges",
            "phi_edges", "event_population_fingerprint",
            "acceptance_feature_metadata_fingerprint",
            "child_assignment_projection_fingerprint", "fingerprint_inputs", "fingerprint",
        ):
            self.assertEqual(runtime_like[field], baseline[field])

    def test_array_like_geometry_validation_and_unavailable_path_are_safe(self):
        cases = (
            ("canonical_t_edges", [0.0, 1.0, 1.0], "canonical_t_edges_not_strictly_increasing"),
            ("delta_edges", [-10.0, float("nan"), 10.0], "delta_edges_nonfinite"),
            ("canonical_t_edges", [2.0, 1.0, 0.0], "canonical_t_edges_not_strictly_increasing"),
        )
        for field, edges, reason in cases:
            with self.subTest(field=field, reason=reason):
                phase, cache = _fixture()
                phase[field] = _AmbiguousTruthArray(edges)
                result = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
                    phase, cache, phi_edges=_AmbiguousTruthArray(PHI_EDGES)
                )
                self.assertFalse(result["available"])
                self.assertEqual(result["status"], "unavailable")
                self.assertEqual(result["reason"], reason)
        phase, cache = _fixture()
        cache["coordinate_fingerprint"] = "bad-coordinates"
        result = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=_AmbiguousTruthArray(PHI_EDGES)
        )
        self.assertFalse(result["available"])
        self.assertEqual(result["reason"], "pion_cache_coordinate_fingerprint_mismatch")
        self.assertEqual(result["phi_edges"], PHI_EDGES)

    def test_early_unavailable_preserves_array_like_phi_geometry(self):
        phase, cache = _fixture()
        phase["schema_version"] = "wrong"
        result = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=_AmbiguousTruthArray(PHI_EDGES)
        )
        self.assertFalse(result["available"])
        self.assertEqual(result["status"], "unavailable")
        self.assertEqual(result["reason"], "phase_a_contract_unavailable")
        self.assertEqual(result["diagnostic_stage"], "validation")
        self.assertEqual(result["phi_edges"], PHI_EDGES)

    def test_complete_join_outside_phi_and_deterministic_fingerprints(self):
        phase, cache = _fixture()
        before = deepcopy((phase, cache))
        first = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=PHI_EDGES
        )
        second = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=PHI_EDGES
        )
        self.assertTrue(first["available"], first.get("reason"))
        self.assertEqual(first["fingerprint"], second["fingerprint"])
        self.assertEqual(first["event_population_fingerprint"], second["event_population_fingerprint"])
        self.assertFalse(first["method_b_numerical_dependency"])
        self.assertEqual(len(first["records"]), 6)
        self.assertEqual(first["records"][-1]["phi_status"], "outside_phi")
        self.assertIsNone(first["records"][-1]["phi_index"])
        self.assertEqual(first["summary"]["prompt_low_count"], 3)
        self.assertEqual(first["summary"]["prompt_control_count"], 2)
        first["records"][0]["SHMS_xptar"] = 99.0
        self.assertEqual((phase, cache), before)

    def test_strict_boundaries_and_identity_parity_fail_locally(self):
        phase, cache = _fixture()
        baseline = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=PHI_EDGES
        )
        self.assertTrue(baseline["available"])
        cases = (
            (lambda p, _c: p["pion_records"].append(deepcopy(p["pion_records"][0])), "phase_a_pion_identity_duplicate"),
            (lambda _p, c: c["records"].pop(), "phase_a_parent_match_missing"),
            (lambda _p, c: (
                c["records"][0].update(ssxptar=float("nan")),
                c["child_event_cache"]["prompt"]["ssxptar"].__setitem__(0, float("nan")),
            ), "acceptance_feature_nonfinite:ssxptar"),
            (lambda _p, c: c["child_event_cache"]["prompt"]["ssdelta"].__setitem__(0, 99.0), "parent_child_parity_mismatch:ssdelta"),
            (lambda _p, c: c["records"][3].update(phi_index=0), "inside_phi_child_missing"),
        )
        for mutate, reason in cases:
            with self.subTest(reason=reason):
                candidate_phase, candidate_cache = _fixture()
                mutate(candidate_phase, candidate_cache)
                result = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
                    candidate_phase, candidate_cache, phi_edges=PHI_EDGES
                )
                self.assertFalse(result["available"])
                self.assertEqual(result["reason"], reason)

    def test_phase_a_and_parent_cache_no_rf_provenance_are_required(self):
        cases = (
            (
                lambda phase, _cache: phase["pion_records"][0].pop(
                    "noRF_provenance"
                ),
                "phase_a_noRF_provenance_invalid",
            ),
            (
                lambda phase, _cache: phase["pion_records"][0].update(
                    noRF_provenance="withRF"
                ),
                "phase_a_noRF_provenance_invalid",
            ),
            (
                lambda _phase, cache: cache["records"][0].pop("rf_state"),
                "pion_cache_parent_fields_missing",
            ),
            (
                lambda _phase, cache: cache["records"][0].update(
                    rf_state="withRF"
                ),
                "pion_cache_parent_noRF_provenance_invalid",
            ),
        )
        for mutate, reason in cases:
            with self.subTest(reason=reason):
                phase, cache = _fixture()
                mutate(phase, cache)
                result = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
                    phase, cache, phi_edges=PHI_EDGES
                )
                self.assertFalse(result["available"])
                self.assertEqual(result["reason"], reason)

    def test_fingerprint_covers_primary_fields_and_artifact_is_json_safe(self):
        phase, cache = _fixture()
        baseline = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            phase, cache, phi_edges=PHI_EDGES
        )
        changed_phase, changed_cache = _fixture()
        changed_cache["records"][0]["ssxptar"] += 0.01
        changed_cache["child_event_cache"]["prompt"]["ssxptar"][0] += 0.01
        changed = acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
            changed_phase, changed_cache, phi_edges=PHI_EDGES
        )
        self.assertTrue(changed["available"])
        self.assertNotEqual(baseline["event_population_fingerprint"], changed["event_population_fingerprint"])
        setting = {
            "kinematic_token": "Q4p4W2p74", "Q2": 4.4, "W": 2.74,
            "epsilon_setting": "low", "epsilon_filename_token": "lowe",
            "phi_setting": "Left", "particle_type": "kaon",
        }
        artifact = acceptance.build_pion_hgcer_method_a_acceptance_contract_artifact(
            setting=setting, contract=baseline
        )
        self.assertEqual(
            acceptance.pion_hgcer_method_a_acceptance_contract_filename(
                "Left", "kaon", "Q4p4W2p74", "lowe"
            ),
            "Left_kaon_pion-background_hgcer_method-a-acceptance-contract_Q4p4W2p74_lowe.json",
        )
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "acceptance.json"
            acceptance.write_pion_hgcer_method_a_acceptance_contract_json(target, artifact)
            self.assertTrue(target.read_text(encoding="utf-8").endswith("\n"))
            self.assertEqual(json.loads(target.read_text(encoding="utf-8"))["contract"]["fingerprint"], baseline["fingerprint"])


if __name__ == "__main__":
    unittest.main()
