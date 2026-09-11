"""Focused F.1 v2 detached dual-population acceptance-contract tests."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))
sys.path.insert(0, str(REPO_ROOT / "src" / "utility"))

import pion_hgcer_method_a_acceptance_contract as acceptance
import pion_hgcer_refinement_method_a as method_a


T_EDGES = [0.0, 1.0, 2.0]
DELTA_EDGES = [-10.0, 0.0, 10.0]
PHI_EDGES = [-20.0, 0.0, 20.0]
COORDINATE_FINGERPRINT = "coordinates"


class _AmbiguousTruthArray:
    """Indexable test stand-in for NumPy's ambiguous array truth value."""

    def __init__(self, values):
        self._values = list(values)

    def __len__(self):
        return len(self._values)

    def __getitem__(self, index):
        return self._values[index]

    def __bool__(self):
        raise ValueError("The truth value of an array is ambiguous.")


class _ArrayScalar:
    """Test-only NumPy-like scalar requiring explicit .item() detachment."""

    def __init__(self, value):
        self._value = value

    def item(self):
        return self._value


class _BrokenArrayScalar:
    def item(self):
        raise RuntimeError("synthetic child scalar failure")


def _coefficient(source):
    return 1.0 if source == "prompt" else -1.0


def _response_record(source, entry, t_index, delta_index, delta, npe, *, nommcuts=True):
    coefficient = _coefficient(source)
    return {
        "side": "pion", "source_label": source, "entry_index": entry,
        "coefficient": coefficient, "diagnostic_weight": coefficient,
        "raw_MM": 1.0 + 0.01 * entry, "raw_t": T_EDGES[t_index] + 0.25,
        "analysis_MM": 1.1 + 0.01 * entry,
        "analysis_t": T_EDGES[t_index] + 0.25,
        "coordinate_fingerprint": COORDINATE_FINGERPRINT,
        "canonical_t_index": t_index, "ssdelta": delta,
        "delta_index": delta_index, "P_hgcer_npeSum": npe,
        "P_hgcer_xAtCer": 1.0 + entry, "P_hgcer_yAtCer": -1.0 - entry,
        "allcuts": True, "nommcuts": nommcuts, "Q2": 4.4, "W": 2.74,
        "epsilon": 0.3, "phi": -10.0 if entry % 2 == 0 else 10.0,
        "proton_cleaning_factor": None, "rf_applied_to_diagnostic": False,
    }


def _acceptance_row(response):
    entry = response["entry_index"]
    return {
        "side": "pion", "source_label": response["source_label"],
        "entry_index": entry,
        "coordinate_fingerprint": response["coordinate_fingerprint"],
        "analysis_t": response["analysis_t"],
        "canonical_t_index": response["canonical_t_index"],
        "ssdelta": response["ssdelta"], "delta_index": response["delta_index"],
        "P_hgcer_npeSum": response["P_hgcer_npeSum"],
        "ssxptar": -0.04 + 0.01 * entry,
        "ssyptar": 0.02 - 0.005 * entry,
        "diagnostic_weight": response["diagnostic_weight"],
        "allcuts": response["allcuts"], "nommcuts": response["nommcuts"],
        "rf_applied_to_diagnostic": False,
    }


def _phase_record(response):
    entry = response["entry_index"]
    t_index, delta_index = response["canonical_t_index"], response["delta_index"]
    baseline_weight = 0.5 + 0.1 * entry
    return {
        "source_label": response["source_label"], "entry_index": entry,
        "coordinate_fingerprint": response["coordinate_fingerprint"],
        "canonical_t_index": t_index,
        "canonical_t_lower_edge": T_EDGES[t_index],
        "canonical_t_upper_edge": T_EDGES[t_index + 1],
        "analysis_abs_t": response["analysis_t"], "analysis_MM": response["analysis_MM"],
        "signed_source_coefficient": response["coefficient"],
        "allcuts": response["allcuts"], "nommcuts": response["nommcuts"],
        "SHMS_delta": response["ssdelta"], "delta_index": delta_index,
        "delta_lower_edge": DELTA_EDGES[delta_index],
        "delta_upper_edge": DELTA_EDGES[delta_index + 1],
        "P_hgcer_npeSum": response["P_hgcer_npeSum"],
        "P_hgcer_xAtCer": response["P_hgcer_xAtCer"],
        "P_hgcer_yAtCer": response["P_hgcer_yAtCer"],
        "baseline_pion_weight_w0": baseline_weight,
        "signed_baseline_event_contribution": response["coefficient"] * baseline_weight,
        "noRF_provenance": "noRF",
    }


def _parent(phase_record, acceptance_row, phi_index):
    return {
        "source_label": phase_record["source_label"],
        "entry_index": phase_record["entry_index"],
        "coordinate_fingerprint": phase_record["coordinate_fingerprint"],
        "rf_state": "noRF", "source_tree_name": "PionEvents_noRF",
        "t_index": phase_record["canonical_t_index"],
        "adj_t": phase_record["analysis_abs_t"], "adj_MM": phase_record["analysis_MM"],
        "coefficient": phase_record["signed_source_coefficient"],
        "allcuts": phase_record["allcuts"], "nommcuts": phase_record["nommcuts"],
        "ssdelta": phase_record["SHMS_delta"], "delta_index": phase_record["delta_index"],
        "P_hgcer_npeSum": phase_record["P_hgcer_npeSum"],
        "P_hgcer_xAtCer": phase_record["P_hgcer_xAtCer"],
        "P_hgcer_yAtCer": phase_record["P_hgcer_yAtCer"],
        "ssxptar": acceptance_row["ssxptar"], "ssyptar": acceptance_row["ssyptar"],
        "hsxptar": 0.001 * phase_record["entry_index"],
        "hsyptar": -0.002 * phase_record["entry_index"],
        "phi_degrees": -10.0 if phi_index == 0 else 10.0 if phi_index == 1 else 30.0,
        "phi_index": phi_index,
    }


def _child(parent):
    return {
        "source_label": parent["source_label"], "entry_index": parent["entry_index"],
        "coefficient": parent["coefficient"],
        "coordinate_fingerprint": parent["coordinate_fingerprint"],
        "adj_t": parent["adj_t"], "adj_MM": parent["adj_MM"],
        "ssdelta": parent["ssdelta"], "P_hgcer_npeSum": parent["P_hgcer_npeSum"],
        "P_hgcer_xAtCer": parent["P_hgcer_xAtCer"],
        "P_hgcer_yAtCer": parent["P_hgcer_yAtCer"],
        "ssxptar": parent["ssxptar"], "ssyptar": parent["ssyptar"],
        "hsxptar": parent["hsxptar"], "hsyptar": parent["hsyptar"],
        "phi_degrees": parent["phi_degrees"], "allcuts": parent["allcuts"],
        "nommcuts": parent["nommcuts"], "t_index": parent["t_index"],
        "phi_index": parent["phi_index"], "Q2": 4.4, "W": 2.74,
        "epsilon": 0.3, "theta_cm_deg": 12.0,
    }


def _diagnostic(records):
    labels = {row["source_label"] for row in records}
    provenance = {
        label: {
            "source_role": label, "coefficient": _coefficient(label),
            "tree_name": "PionEvents_{}_noRF".format(label), "rf_state": "noRF",
            "proton_factor_scope": "none",
        }
        for label in labels
    }
    return {
        "status": "available", "non_authoritative": True,
        "production_side_effect_free": True, "rf_restoration_applied": False,
        "coordinate_fingerprint": COORDINATE_FINGERPRINT,
        "config_fingerprint": "part1-config", "t_edges": list(T_EDGES),
        "delta_edges": list(DELTA_EDGES), "source_provenance": {"pion": provenance},
        "records": {"pion": tuple(records), "kaon": ()},
        "phase_e_acceptance_records": {
            "pion": tuple(_acceptance_row(row) for row in records), "kaon": (),
        },
    }


def _fixture():
    response_rows = [
        _response_record("prompt", 0, 0, 0, -5.0, 1.5),
        _response_record("prompt", 1, 0, 1, 5.0, 3.0),
        _response_record("prompt", 2, 1, 0, -5.0, 0.0),
        _response_record("prompt", 3, 1, 1, 5.0, 1.2, nommcuts=False),
        _response_record("rand", 4, 1, 1, 5.0, 3.5),
    ]
    diagnostic = _diagnostic(response_rows)
    application_response_rows = [response_rows[1], response_rows[4]]
    phase_rows = [_phase_record(row) for row in application_response_rows]
    acceptance_by_identity = {
        (row["source_label"], row["entry_index"]): row
        for row in diagnostic["phase_e_acceptance_records"]["pion"]
    }
    parents = [
        _parent(row, acceptance_by_identity[(row["source_label"], row["entry_index"])], phi)
        for row, phi in zip(phase_rows, (0, 1))
    ]
    children = [_child(parent) for parent in parents]
    child_cache = {}
    for child in children:
        section = child_cache.setdefault(child["source_label"], {name: [] for name in child})
        for name, value in child.items():
            section[name].append(value)
    phase = {
        "schema_version": "pion_hgcer_event_contract/v1",
        "fingerprint_schema_version": "pion_hgcer_event_contract_fingerprint/v2",
        "status": "available", "available": True,
        "contract_fingerprint": "phase-contract",
        "pion_event_population_fingerprint": "phase-pion-events",
        "physical_pion_control_mask_fingerprint": "physical-npe-gt-2",
        "coordinate_fingerprint": COORDINATE_FINGERPRINT,
        "canonical_t_edges": list(T_EDGES), "delta_edges": list(DELTA_EDGES),
        "pion_records": phase_rows, "host_state": "proton_cleaned",
        "source_target_state": "post_proton_noRF", "rf_restoration_applied": False,
        "immutable_record_contract": True, "production_objects_mutated": False,
        "refinement_applied": False,
    }
    cache = {
        "coordinate_fingerprint": COORDINATE_FINGERPRINT,
        "delta_edges": list(DELTA_EDGES), "records": parents,
        "child_event_cache": child_cache,
    }
    method = method_a.build_pion_hgcer_method_a(diagnostic, phase)
    assert method["available"], method
    return diagnostic, method, phase, cache


def _build(diagnostic, method, phase, cache, *, phi_edges=PHI_EDGES):
    return acceptance.build_pion_hgcer_method_a_acceptance_event_contract(
        diagnostic, method, phase, cache, phi_edges=phi_edges
    )


def _runtime_like_child_scalars(cache):
    for section in cache["child_event_cache"].values():
        for field, values in section.items():
            section[field] = [_ArrayScalar(value) for value in values]


class MethodAAcceptanceContractTests(unittest.TestCase):
    def test_dual_population_method_a_closure_and_actual_method_a_integration(self):
        diagnostic, method, phase, cache = _fixture()
        before = deepcopy((diagnostic, method, phase, cache))
        result = _build(diagnostic, method, phase, cache)
        self.assertTrue(result["available"], result.get("reason"))
        self.assertEqual(result["schema_version"], "pion_hgcer_method_a_acceptance_event_contract/v2")
        self.assertEqual(result["fingerprint_schema_version"], "pion_hgcer_method_a_acceptance_event_contract_fingerprint/v2")
        self.assertEqual(len(result["method_a_training_records"]), 2)
        self.assertEqual(len(result["application_records"]), 2)
        self.assertEqual([row["response_class"] for row in result["method_a_training_records"]], ["low", "control"])
        self.assertEqual(
            {(row["source_label"], row["entry_index"]) for row in result["method_a_training_records"]},
            {("prompt", 0), ("prompt", 1)},
        )
        self.assertEqual(
            {(row["source_label"], row["entry_index"]) for row in result["application_records"]},
            {("prompt", 1), ("rand", 4)},
        )
        self.assertTrue(all(row["P_hgcer_npeSum"] > 2.0 for row in result["application_records"]))
        training = result["method_a_training_summary"]
        self.assertEqual(
            (training["prompt_positive_count"], training["prompt_low_count"], training["prompt_control_count"]),
            (2, 1, 1),
        )
        self.assertTrue(training["method_a_closure_passed"])
        self.assertEqual(training["response_threshold_audit"]["observed_nonpositive_response_record_count"], 1)
        self.assertFalse(training["response_threshold_audit"]["absolute_leakage_probability_claimed"])
        self.assertEqual((diagnostic, method, phase, cache), before)

    def test_training_and_application_fingerprints_are_independent(self):
        diagnostic, method, phase, cache = _fixture()
        baseline = _build(diagnostic, method, phase, cache)
        changed_diagnostic = deepcopy(diagnostic)
        changed_diagnostic["records"]["pion"] = tuple(
            dict(row, P_hgcer_xAtCer=row["P_hgcer_xAtCer"] + 0.5)
            if row["entry_index"] == 0 else row
            for row in changed_diagnostic["records"]["pion"]
        )
        training_changed = _build(changed_diagnostic, method, phase, cache)
        self.assertTrue(training_changed["available"])
        self.assertNotEqual(baseline["method_a_training_population_fingerprint"], training_changed["method_a_training_population_fingerprint"])
        self.assertEqual(baseline["application_population_fingerprint"], training_changed["application_population_fingerprint"])
        changed_phase, changed_cache = deepcopy(phase), deepcopy(cache)
        changed_phase["pion_records"][0]["P_hgcer_xAtCer"] += 0.25
        changed_cache["records"][0]["P_hgcer_xAtCer"] += 0.25
        changed_cache["child_event_cache"]["prompt"]["P_hgcer_xAtCer"][0] += 0.25
        application_changed = _build(diagnostic, method, changed_phase, changed_cache)
        self.assertTrue(application_changed["available"])
        self.assertEqual(baseline["method_a_training_population_fingerprint"], application_changed["method_a_training_population_fingerprint"])
        self.assertNotEqual(baseline["application_population_fingerprint"], application_changed["application_population_fingerprint"])
        reassigned_phase, reassigned_cache = deepcopy(phase), deepcopy(cache)
        reassigned_cache["records"][0].update(phi_degrees=10.0, phi_index=1)
        reassigned_cache["child_event_cache"]["prompt"]["phi_degrees"][0] = 10.0
        reassigned_cache["child_event_cache"]["prompt"]["phi_index"][0] = 1
        reassigned = _build(diagnostic, method, reassigned_phase, reassigned_cache)
        self.assertTrue(reassigned["available"])
        self.assertNotEqual(
            baseline["application_child_assignment_projection_fingerprint"],
            reassigned["application_child_assignment_projection_fingerprint"],
        )

    def test_response_acceptance_join_and_parity_fail_closed(self):
        cases = (
            (lambda d, _m, _p, _c: d["records"].update(pion=d["records"]["pion"] + (d["records"]["pion"][0],)), "part1_response_records_identity_duplicate"),
            (lambda d, _m, _p, _c: d["phase_e_acceptance_records"].update(pion=d["phase_e_acceptance_records"]["pion"][:-1]), "part1_acceptance_match_missing"),
            (lambda d, _m, _p, _c: d["phase_e_acceptance_records"].update(pion=d["phase_e_acceptance_records"]["pion"] + (dict(d["phase_e_acceptance_records"]["pion"][0], entry_index=99),)), "part1_acceptance_match_extra"),
            (lambda d, _m, _p, _c: d["phase_e_acceptance_records"]["pion"][0].update(allcuts=False), "part1_response_acceptance_parity_mismatch:allcuts"),
            (lambda d, _m, _p, _c: d["phase_e_acceptance_records"]["pion"][0].update(nommcuts=False), "part1_response_acceptance_parity_mismatch:nommcuts"),
            (lambda d, _m, _p, _c: d["phase_e_acceptance_records"]["pion"][0].update(diagnostic_weight=9.0), "part1_response_acceptance_parity_mismatch:diagnostic_weight"),
            (lambda d, _m, _p, _c: d["records"]["pion"][0].update(rf_applied_to_diagnostic=True), "part1_response_acceptance_parity_mismatch:rf_applied_to_diagnostic"),
            (lambda d, _m, _p, _c: d["phase_e_acceptance_records"]["pion"][0].update(ssxptar=float("nan")), "acceptance_feature_nonfinite:ssxptar"),
            (lambda d, _m, _p, _c: d["records"]["pion"][0].update(P_hgcer_yAtCer=float("nan")), "acceptance_feature_nonfinite:P_hgcer_yAtCer"),
        )
        for mutate, reason in cases:
            with self.subTest(reason=reason):
                diagnostic, method, phase, cache = _fixture()
                mutate(diagnostic, method, phase, cache)
                result = _build(diagnostic, method, phase, cache)
                self.assertFalse(result["available"])
                self.assertEqual(result["reason"], reason)

    def test_method_a_provenance_and_closure_fail_closed(self):
        mutations = (
            lambda m: m.update(schema_version="wrong"),
            lambda m: m.update(available=False),
            lambda m: m.update(coordinate_fingerprint="wrong"),
            lambda m: m.update(phase_a_contract_fingerprint="wrong"),
            lambda m: m.update(response_population_definition="wrong"),
            lambda m: m["cells"][0].update(prompt_low_count=99),
            lambda m: m["cells"][0].update(partition_closure_passed=False),
            lambda m: m["summary"].update(prompt_positive_nommcuts_records=99),
        )
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                diagnostic, method, phase, cache = _fixture()
                changed = deepcopy(method)
                mutate(changed)
                self.assertFalse(_build(diagnostic, changed, phase, cache)["available"])

    def test_application_population_rejects_low_npe_and_preserves_baseline_fields(self):
        diagnostic, method, phase, cache = _fixture()
        baseline = _build(diagnostic, method, phase, cache)
        self.assertTrue(baseline["available"])
        self.assertEqual(baseline["application_records"][0]["baseline_pion_weight_w0"], phase["pion_records"][0]["baseline_pion_weight_w0"])
        phase, cache = deepcopy(phase), deepcopy(cache)
        phase["pion_records"][0]["P_hgcer_npeSum"] = 2.0
        cache["records"][0]["P_hgcer_npeSum"] = 2.0
        cache["child_event_cache"]["prompt"]["P_hgcer_npeSum"][0] = 2.0
        result = _build(diagnostic, method, phase, cache)
        self.assertFalse(result["available"])
        self.assertEqual(result["reason"], "application_population_npe_not_physical_control")

    def test_array_and_scalar_container_invariance_and_scalar_failure(self):
        diagnostic, method, phase, cache = _fixture()
        baseline = _build(diagnostic, method, phase, cache)
        changed_diagnostic, changed_phase, changed_cache = deepcopy(diagnostic), deepcopy(phase), deepcopy(cache)
        changed_diagnostic["t_edges"] = _AmbiguousTruthArray(T_EDGES)
        changed_diagnostic["delta_edges"] = _AmbiguousTruthArray(DELTA_EDGES)
        changed_phase["canonical_t_edges"] = _AmbiguousTruthArray(T_EDGES)
        changed_phase["delta_edges"] = _AmbiguousTruthArray(DELTA_EDGES)
        changed_cache["delta_edges"] = _AmbiguousTruthArray(DELTA_EDGES)
        _runtime_like_child_scalars(changed_cache)
        changed = _build(changed_diagnostic, method, changed_phase, changed_cache, phi_edges=_AmbiguousTruthArray(PHI_EDGES))
        self.assertTrue(changed["available"], changed.get("reason"))
        for name in (
            "method_a_training_records", "application_records",
            "method_a_training_population_fingerprint", "application_population_fingerprint",
            "application_child_assignment_projection_fingerprint", "fingerprint",
        ):
            self.assertEqual(changed[name], baseline[name])
        bad_diagnostic, bad_method, bad_phase, bad_cache = _fixture()
        _runtime_like_child_scalars(bad_cache)
        bad_cache["child_event_cache"]["prompt"]["entry_index"][0] = _BrokenArrayScalar()
        bad = _build(bad_diagnostic, bad_method, bad_phase, bad_cache)
        self.assertFalse(bad["available"])
        self.assertEqual(bad["reason"], "pion_cache_child_scalar_invalid")

    def test_artifact_json_safety_and_filename_are_v2(self):
        diagnostic, method, phase, cache = _fixture()
        contract = _build(diagnostic, method, phase, cache)
        setting = {
            "kinematic_token": "Q4p4W2p74", "Q2": 4.4, "W": 2.74,
            "epsilon_setting": "low", "epsilon_filename_token": "lowe",
            "phi_setting": "Left", "particle_type": "kaon",
        }
        artifact = acceptance.build_pion_hgcer_method_a_acceptance_contract_artifact(setting=setting, contract=contract)
        self.assertEqual(artifact["schema_version"], "pion_hgcer_method_a_acceptance_event_contract_artifact/v2")
        self.assertEqual(
            acceptance.pion_hgcer_method_a_acceptance_contract_filename("Left", "kaon", "Q4p4W2p74", "lowe"),
            "Left_kaon_pion-background_hgcer_method-a-acceptance-contract_Q4p4W2p74_lowe.json",
        )
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "acceptance.json"
            acceptance.write_pion_hgcer_method_a_acceptance_contract_json(target, artifact)
            self.assertTrue(target.read_text(encoding="utf-8").endswith("\n"))
            self.assertEqual(json.loads(target.read_text(encoding="utf-8"))["contract"]["fingerprint"], contract["fingerprint"])


if __name__ == "__main__":
    unittest.main()
