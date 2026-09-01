"""Focused Phase-A event-contract tests with optional PyROOT closure."""

from __future__ import annotations

import copy
import json
import math
import subprocess
import sys
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

import pion_hgcer_event_contract as event_contract
import pion_hgcer_refinement_plots as refinement_plots

try:
    import ROOT
except ImportError:  # pragma: no cover - reported as a gated skip.
    ROOT = None


class _Axis:
    def __init__(self, edges):
        self.edges = [float(value) for value in edges]

    def GetXmin(self):
        return self.edges[0]

    def GetXmax(self):
        return self.edges[-1]

    def GetBinLowEdge(self, index):
        index = int(index)
        if index <= 0:
            return self.edges[0] - (self.edges[1] - self.edges[0])
        if index <= len(self.edges):
            return self.edges[index - 1]
        return self.edges[-1]

    def GetBinUpEdge(self, index):
        index = int(index)
        if index <= 0:
            return self.edges[0]
        if index < len(self.edges):
            return self.edges[index]
        return self.edges[-1] + (self.edges[-1] - self.edges[-2])

    def FindBin(self, value):
        value = float(value)
        if value < self.edges[0]:
            return 0
        if value >= self.edges[-1]:
            return len(self.edges)
        for index, (low, high) in enumerate(
            zip(self.edges[:-1], self.edges[1:]), start=1
        ):
            if low <= value < high:
                return index
        return len(self.edges)


class _Histogram:
    counter = 0

    def __init__(self, name, edges=(0.0, 1.0, 2.0)):
        self.name = str(name)
        self.axis = _Axis(edges)
        self.contents = [0.0] * (len(edges) + 1)
        self.sumw2 = [0.0] * (len(edges) + 1)

    def GetName(self):
        return self.name

    def GetNbinsX(self):
        return len(self.contents) - 2

    def GetXaxis(self):
        return self.axis

    def GetBinContent(self, index):
        return self.contents[int(index)]

    def GetBinError(self, index):
        return math.sqrt(max(0.0, self.sumw2[int(index)]))

    def Fill(self, value, weight=1.0):
        index = self.axis.FindBin(value)
        weight = float(weight)
        self.contents[index] += weight
        self.sumw2[index] += weight * weight

    def Clone(self, name):
        cloned = _Histogram(name, self.axis.edges)
        cloned.contents = list(self.contents)
        cloned.sumw2 = list(self.sumw2)
        return cloned

    def Reset(self):
        self.contents = [0.0] * len(self.contents)
        self.sumw2 = [0.0] * len(self.sumw2)

    def Sumw2(self):
        return None

    def SetDirectory(self, _directory):
        return None

    def Add(self, other, scale=1.0):
        scale = float(scale)
        for index in range(len(self.contents)):
            self.contents[index] += scale * other.contents[index]
            self.sumw2[index] += scale * scale * other.sumw2[index]

    def Integral(self, first=None, last=None):
        if first is None:
            first, last = 1, self.GetNbinsX()
        return sum(self.contents[int(first):int(last) + 1])

    def InheritsFrom(self, name):
        return name in ("TH1", "TObject")


def _sum_histograms(name, histograms):
    result = histograms[0].Clone(name)
    result.Reset()
    for histogram in histograms:
        result.Add(histogram)
    return result


def _fill_template(name, records, reference, weights, selection):
    histogram = _Histogram(name, reference.axis.edges)
    for record in records:
        if not record[selection]:
            continue
        w0 = event_contract.simc_shape_pion_weight_from_value(
            record["adj_MM"], reference, weights
        )
        contribution = float(record["coefficient"]) * float(w0)
        if contribution:
            histogram.Fill(record["adj_MM"], contribution)
    return histogram


def _fill_host(name, records, selection):
    histogram = _Histogram(name)
    for record in records:
        if record[selection]:
            contribution = record["coefficient"] * record["final_cleaned_factor"]
            if contribution:
                histogram.Fill(record["adj_mm"], contribution)
    return histogram


class PhaseAEventContractTests(unittest.TestCase):
    def _fixture(self, *, identity_host=False):
        coordinate = "coordinate-phase-a"
        pair_id = "pair-phase-a"
        pair_hash = "pair-hash-phase-a"
        t_edges = [0.0, 1.0, 2.0]
        delta_edges = [-10.0, 0.0, 10.0]
        pion_by_t = [
            [
                {
                    "source_label": "prompt", "entry_index": 0,
                    "coefficient": 2.0, "source_tree_name": "Cut_Pion_prompt_noRF",
                    "rf_state": "noRF", "raw_MM": 0.25, "adj_MM": 0.25,
                    "raw_t": 1.4, "adj_t": 0.4, "ssdelta": -5.0,
                    "t_index": 0, "delta_index": 0, "P_hgcer_npeSum": 2.5,
                    "P_hgcer_xAtCer": 1.25, "P_hgcer_yAtCer": -2.5,
                    "allcuts": True, "nommcuts": True,
                },
                {
                    "source_label": "rand", "entry_index": 1,
                    "coefficient": -0.5, "source_tree_name": "Cut_Pion_rand_noRF",
                    "rf_state": "noRF", "raw_MM": 1.25, "adj_MM": 1.25,
                    "raw_t": 0.35, "adj_t": 0.35, "ssdelta": 12.0,
                    "t_index": 0, "delta_index": None, "P_hgcer_npeSum": 3.0,
                    "P_hgcer_xAtCer": 2.0, "P_hgcer_yAtCer": -1.0,
                    "allcuts": False, "nommcuts": True,
                },
            ],
            [
                {
                    "source_label": "dummy", "entry_index": 2,
                    "coefficient": -1.0, "source_tree_name": "Cut_Pion_prompt_noRF",
                    "rf_state": "noRF", "raw_MM": 0.75, "adj_MM": 0.75,
                    "raw_t": 1.2, "adj_t": 1.2, "ssdelta": 10.0,
                    "t_index": 1, "delta_index": 1, "P_hgcer_npeSum": 5.0,
                    "P_hgcer_xAtCer": -3.0, "P_hgcer_yAtCer": 4.0,
                    "allcuts": True, "nommcuts": True,
                },
                {
                    "source_label": "dummy_rand", "entry_index": 3,
                    "coefficient": 0.25, "source_tree_name": "Cut_Pion_rand_noRF",
                    "rf_state": "noRF", "raw_MM": 1.75, "adj_MM": 1.75,
                    "raw_t": 1.7, "adj_t": 1.7, "ssdelta": 0.0,
                    "t_index": 1, "delta_index": 1, "P_hgcer_npeSum": 8.0,
                    "P_hgcer_xAtCer": 0.0, "P_hgcer_yAtCer": 0.5,
                    "allcuts": True, "nommcuts": True,
                },
            ],
        ]
        references = [_Histogram("weight_reference_t1"), _Histogram("weight_reference_t2")]
        weights = [[0.0, 2.0, 3.0, 0.0], [0.0, 4.0, 5.0, 0.0]]
        parents = []
        full_templates = []
        cut_templates = []
        for t_index in range(2):
            full = _fill_template(
                "pion_full_t{}".format(t_index + 1), pion_by_t[t_index],
                references[t_index], weights[t_index], "nommcuts",
            )
            cut = _fill_template(
                "pion_cut_t{}".format(t_index + 1), pion_by_t[t_index],
                references[t_index], weights[t_index], "allcuts",
            )
            weight_payload = {
                "H_pion_control_model": references[t_index],
                "weights": list(weights[t_index]),
                "diagnostics": {"model_variant": "final", "clip_min": 0.0},
            }
            parent = {
                "t_bin_index": t_index,
                "t_edges": t_edges[t_index:t_index + 2],
                "pion_parent_id": "pion-parent-t{}".format(t_index + 1),
                "parent_physics_identity": {
                    "Q2": "4p4", "W": "2p74", "phi_setting": "Left",
                    "t_bin_index": t_index,
                },
                "parent_fit_configuration_hash": "fit-hash-t{}".format(t_index + 1),
                "coordinate_fingerprint": coordinate,
                "canonical_interval_pair_id": pair_id,
                "canonical_interval_pair_hash": pair_hash,
                "fit_result": {
                    "analysis_scope": "t_bin{}".format(t_index + 1),
                    "_phase_a_weight_payload": weight_payload,
                },
                "final_diagnostic_application_result": {
                    **weight_payload,
                    "H_pion_subtraction_template_MM": cut,
                    "H_pion_subtraction_template_MM_nosub": full,
                },
                "final_diagnostic_application_status": {
                    "final_status": "applied_component"
                },
            }
            parents.append(parent)
            full_templates.append(full)
            cut_templates.append(cut)

        source_accounting = {}
        for records in pion_by_t:
            for record in records:
                source_accounting.setdefault(record["source_label"], {
                    "coefficient": record["coefficient"],
                    "tree_name": record["source_tree_name"],
                    "rf_state": "noRF",
                })
        pion_cache = {
            "coordinate_fingerprint": coordinate,
            "physical_pion_control_mask_fingerprint": "pion-mask-phase-a",
            "delta_edges": delta_edges,
            "source_accounting": source_accounting,
            "by_t": [
                {
                    "t_index": index,
                    "t_edges": t_edges[index:index + 2],
                    "coordinate_fingerprint": coordinate,
                    "records": tuple(records),
                }
                for index, records in enumerate(pion_by_t)
            ],
        }

        host_records = {
            "prompt": [
                {
                    "entry_index": 0, "adj_mm": 0.5, "adj_t": 0.3,
                    "delta_value": -5.0, "allcuts": True, "nommcuts": True,
                    "cleaned_factor": 0.6, "final_cleaned_factor": 0.6,
                    "rf_accept": True,
                },
                {
                    "entry_index": 4, "adj_mm": -0.2, "adj_t": 0.6,
                    "delta_value": -5.0, "allcuts": True, "nommcuts": True,
                    "cleaned_factor": 0.2, "final_cleaned_factor": 0.2,
                    "rf_accept": True,
                },
                {
                    "entry_index": 5, "adj_mm": 0.4, "adj_t": 2.5,
                    "delta_value": -5.0, "allcuts": True, "nommcuts": True,
                    "cleaned_factor": 0.5, "final_cleaned_factor": 0.5,
                    "rf_accept": True,
                },
            ],
            "rand": [
                {
                    "entry_index": 1, "adj_mm": 1.5, "adj_t": 0.7,
                    "delta_value": 3.0, "allcuts": False, "nommcuts": True,
                    "cleaned_factor": 1.0, "final_cleaned_factor": 1.0,
                    "rf_accept": True,
                },
            ],
            "dummy": [
                {
                    "entry_index": 2, "adj_mm": 0.8, "adj_t": 1.3,
                    "delta_value": 10.0, "allcuts": True, "nommcuts": True,
                    "cleaned_factor": 0.25, "final_cleaned_factor": 0.25,
                    "rf_accept": True,
                },
                {
                    "entry_index": 6, "adj_mm": 2.2, "adj_t": 1.6,
                    "delta_value": 4.0, "allcuts": True, "nommcuts": True,
                    "cleaned_factor": 0.4, "final_cleaned_factor": 0.4,
                    "rf_accept": False,
                },
            ],
            "dummy_rand": [
                {
                    "entry_index": 3, "adj_mm": 1.8, "adj_t": 1.8,
                    "delta_value": 0.0, "allcuts": True, "nommcuts": True,
                    "cleaned_factor": 0.8, "final_cleaned_factor": 0.8,
                    "rf_accept": True,
                },
            ],
        }
        coefficients = {
            "prompt": 2.0, "rand": -0.5, "dummy": -1.0,
            "dummy_rand": 0.25,
        }
        if identity_host:
            for records in host_records.values():
                for row in records:
                    row["cleaned_factor"] = 1.0
                    row["final_cleaned_factor"] = 1.0
        prepared_sources = {}
        lookup = {}
        per_t_host = [[], []]
        for label, records in host_records.items():
            entries = {}
            for row in records:
                entries[row["entry_index"]] = {
                    key: value for key, value in row.items()
                    if key not in (
                        "entry_index", "cleaned_factor", "final_cleaned_factor", "rf_accept"
                    )
                }
                t_index = event_contract.find_canonical_bin(row["adj_t"], t_edges)
                delta_index = event_contract.find_canonical_bin(row["delta_value"], delta_edges)
                lookup["{}:{}".format(label, row["entry_index"])] = {
                    "t_index": int(t_index), "delta_index": int(delta_index),
                    "proton_weight": 1.0 - row["cleaned_factor"],
                    "applied_proton_probability": 1.0 - row["cleaned_factor"],
                    "cleaned_factor": row["cleaned_factor"],
                    "final_cleaned_factor": row["final_cleaned_factor"],
                    "rf_accept": row["rf_accept"],
                }
                physical = {
                    **row,
                    "coefficient": coefficients[label],
                }
                if t_index >= 0:
                    per_t_host[t_index].append(physical)
            prepared_sources[label] = {
                "coefficient": coefficients[label],
                "tree_name": "Cut_Kaon_{}_noRF".format(label),
                "entries": entries,
            }
        proton_products = []
        host_full = []
        host_cut = []
        for t_index, rows in enumerate(per_t_host):
            full = _fill_host("host_full_t{}".format(t_index + 1), rows, "nommcuts")
            cut = _fill_host("host_cut_t{}".format(t_index + 1), rows, "allcuts")
            host_full.append(full)
            host_cut.append(cut)
            proton_products.append({
                "t_index": t_index,
                "t_edges": t_edges[t_index:t_index + 2],
                "coordinate_fingerprint": coordinate,
                "final_targets": {"h_mm_nosub": full, "h_mm": cut},
                "final_output_fingerprint": (
                    event_contract.fingerprint_histogram_content_error(full)
                ),
            })
            parents[t_index]["H_proton_cleaned_final_rf"] = full
            parents[t_index]["proton_output_fingerprint"] = (
                event_contract.fingerprint_histogram_content_error(full)
            )
        proton_application = {
            "accepted": not identity_host,
            "coordinate_fingerprint": coordinate,
            "canonical_t_products": tuple(proton_products),
            "final_targets": {
                "h_mm_nosub": _sum_histograms("host_global_full", host_full),
                "h_mm": _sum_histograms("host_global_cut", host_cut),
            },
            "diagnostics": {
                "rf_applied": False,
                "event_weight_source": "setting_wide_immutable_prepared_lookup",
                "lambda_preservation_gate": {
                    "status": "bypassed" if identity_host else "passed",
                    "production_action": "bypass" if identity_host else "apply",
                    "proton_cleaning_committed": not identity_host,
                },
            },
        }
        if identity_host:
            proton_application.update({
                "host_state": "identity_no_proton_cleaning",
                "rf_restoration_applied": False,
                "_prepared_event_weight_lookup": lookup,
            })
        proton_result = {
            "accepted": not identity_host,
            "method": "timing_t_event_weight",
            "coordinate_fingerprint": coordinate,
            "diagnostics": {
                "event_weight_source": "setting_wide_immutable_prepared_lookup",
                "rf_applied": False,
                "lambda_preservation_gate": {
                    "status": "bypassed" if identity_host else "passed",
                    "production_action": "bypass" if identity_host else "apply",
                    "proton_cleaning_committed": not identity_host,
                },
            },
        }
        if not identity_host:
            proton_result["_prepared_event_weight_lookup"] = lookup
        canonical_global = {
            "H_MM_estimated_contamination": _sum_histograms(
                "pion_global_full", full_templates
            )
        }
        canonical_binning = {
            "t_edges": t_edges,
            "phi_edges": [-180.0, 0.0, 180.0],
            "requested_num_phi_bins": 2,
            "actual_num_phi_bins": 2,
            "canonical_interval_pair_id": pair_id,
            "canonical_interval_pair_hash": pair_hash,
        }
        return {
            "pion_control_cache": pion_cache,
            "pion_parents": tuple(parents),
            "canonical_t_global": canonical_global,
            "proton_source_bundle": {"prepared_sources": prepared_sources},
            "proton_cleaning_result": proton_result,
            "proton_cleaning_application": proton_application,
            "inp_dict": {},
            "canonical_binning": canonical_binning,
            "delta_edge_source": "proton_cleaning_result.delta_edges",
        }

    def _build(self, fixture, actions=None):
        actions = actions or {0: "component_weight", 1: "component_weight"}

        def policy(parent, _inp):
            action = actions[int(parent["t_bin_index"])]
            return {
                "action": action,
                "fallback_mode": None if action == "component_weight" else action,
                "fit_accepted": action == "component_weight",
                "parent_id": parent["pion_parent_id"],
            }

        def weights(fit_result, **_kwargs):
            return fit_result["_phase_a_weight_payload"]

        with mock.patch.object(
            event_contract, "resolve_frozen_parent_application_policy", side_effect=policy
        ), mock.patch.object(
            event_contract, "build_simc_shape_pion_control_weights", side_effect=weights
        ), mock.patch.object(
            event_contract, "resolve_particle_subtraction_weight_clip_bounds",
            return_value=(0.0, None),
        ), mock.patch.object(
            event_contract, "resolve_particle_subtraction_weight_denominator_floor",
            return_value=1.0e-12,
        ):
            return event_contract.build_pion_hgcer_event_contract(**fixture)

    def _accepted_application_fixture(self, *, production_action):
        fixture = self._fixture(identity_host=True)
        application = copy.deepcopy(fixture["proton_cleaning_application"])
        result = copy.deepcopy(fixture["proton_cleaning_result"])
        upstream = copy.deepcopy(application["final_targets"])
        application["accepted"] = True
        application.pop("host_state", None)
        application["rf_restoration_applied"] = False
        application["raw_targets"] = copy.deepcopy(upstream)
        application["proton_targets"] = copy.deepcopy(upstream)
        application["cleaned_targets_pre_rf"] = copy.deepcopy(upstream)
        for product in application["canonical_t_products"]:
            final_targets = product["final_targets"]
            product["raw_targets"] = copy.deepcopy(final_targets)
            product["proton_targets"] = copy.deepcopy(final_targets)
            product["cleaned_targets_pre_rf"] = copy.deepcopy(final_targets)
        committed = production_action == "apply"
        gate = {
            "status": "pass" if committed else "fail",
            "production_action": production_action,
            "proton_cleaning_committed": committed,
        }
        application["diagnostics"]["lambda_preservation_gate"] = dict(gate)
        result["accepted"] = True
        result["diagnostics"]["lambda_preservation_gate"] = dict(gate)
        return result, application, upstream

    def test_exact_event_identity_weight_source_accounting_and_closure(self):
        fixture = self._fixture()
        parent_ids_before = [id(parent) for parent in fixture["pion_parents"]]
        hist_fingerprints_before = [
            event_contract.fingerprint_histogram_content_error(
                parent["final_diagnostic_application_result"][
                    "H_pion_subtraction_template_MM_nosub"
                ]
            )
            for parent in fixture["pion_parents"]
        ]
        host_fingerprints_before = [
            event_contract.fingerprint_histogram_content_error(
                parent["H_proton_cleaned_final_rf"]
            )
            for parent in fixture["pion_parents"]
        ]
        lookup_before = copy.deepcopy(
            fixture["proton_cleaning_result"]["_prepared_event_weight_lookup"]
        )
        contract = self._build(fixture)
        self.assertTrue(contract["available"], contract.get("reason"))
        self.assertEqual(
            contract["fingerprint_schema_version"],
            event_contract.EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION,
        )
        self.assertEqual(
            contract["fingerprint_inputs"]["fingerprint_schema_version"],
            event_contract.EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION,
        )
        self.assertTrue(contract["pion_closure"]["passed"])
        self.assertTrue(contract["host_closure"]["passed"])
        self.assertEqual(len(contract["pion_records"]), 4)
        prompt = next(
            record for record in contract["pion_records"]
            if record["source_label"] == "prompt"
        )
        self.assertEqual(prompt["canonical_t_index"], 0)
        self.assertEqual(prompt["raw_t"], 1.4)
        self.assertEqual(prompt["baseline_pion_weight_w0"], 2.0)
        self.assertEqual(prompt["signed_baseline_event_contribution"], 4.0)
        self.assertEqual(prompt["source_coefficient"], 2.0)
        self.assertEqual(prompt["signed_source_coefficient"], 2.0)
        self.assertIsNone(prompt["proton_cleaning_factor"])
        self.assertEqual(prompt["P_hgcer_xAtCer"], 1.25)
        self.assertTrue(all(
            record["P_hgcer_npeSum"] > 2.0
            for record in contract["pion_records"]
        ))
        self.assertTrue(all(
            provenance["authoritative_weight_source"]
            == "frozen_final_diagnostic_application_result"
            for provenance in contract["baseline_weight_provenance"]
        ))
        negative = next(
            record for record in contract["pion_records"]
            if record["source_label"] == "rand"
        )
        self.assertEqual(negative["signed_baseline_event_contribution"], -1.5)
        self.assertEqual(negative["refinement_geometry_status"], "delta_overflow")
        final_edge = next(
            record for record in contract["pion_records"]
            if record["source_label"] == "dummy"
        )
        self.assertEqual(final_edge["delta_index"], 1)
        metrics = contract["pion_closure"]["per_t"][0]["source_accounting"]["rand"]
        self.assertEqual(metrics["signed_weighted_sum"], -1.5)
        self.assertEqual(metrics["sumw2"], 2.25)
        self.assertEqual(metrics["effective_entries"], 1.0)
        self.assertTrue(
            contract["pion_closure"]["source_accounting_closure"]["passed"]
        )
        self.assertEqual(parent_ids_before, [id(parent) for parent in fixture["pion_parents"]])
        self.assertEqual(
            hist_fingerprints_before,
            [
                event_contract.fingerprint_histogram_content_error(
                    parent["final_diagnostic_application_result"][
                        "H_pion_subtraction_template_MM_nosub"
                    ]
                )
                for parent in fixture["pion_parents"]
            ],
        )
        self.assertEqual(
            host_fingerprints_before,
            [
                event_contract.fingerprint_histogram_content_error(
                    parent["H_proton_cleaned_final_rf"]
                )
                for parent in fixture["pion_parents"]
            ],
        )
        self.assertEqual(
            lookup_before,
            fixture["proton_cleaning_result"]["_prepared_event_weight_lookup"],
        )
        json.dumps(contract, allow_nan=False)

    def test_host_uses_existing_final_factor_and_no_pion_refinement(self):
        contract = self._build(self._fixture())
        prompt = next(
            record for record in contract["kaon_host_records"]
            if record["source_label"] == "prompt" and record["analysis_MM"] == 0.5
        )
        self.assertEqual(prompt["proton_cleaning_factor"], 0.6)
        self.assertEqual(prompt["final_cleaned_factor"], 0.6)
        self.assertEqual(prompt["signed_host_event_contribution"], 1.2)
        self.assertIsNone(prompt["pion_refinement_factor"])
        self.assertEqual(prompt["source_target_state"], "post_proton_noRF")
        self.assertEqual(prompt["host_state"], "proton_cleaned")
        self.assertFalse(prompt["rf_restoration_applied"])
        self.assertNotIn("proton_rf_accept", prompt)
        ignored_legacy_rf = next(
            record for record in contract["kaon_host_records"]
            if record["source_label"] == "dummy" and record["analysis_MM"] == 2.2
        )
        self.assertEqual(ignored_legacy_rf["signed_host_event_contribution"], -0.4)
        self.assertEqual(len(contract["host_records_outside_geometry"]), 1)
        global_snapshot = contract["host_closure"]["global_full"]["reconstructed"]
        self.assertNotEqual(global_snapshot["underflow"]["content"], 0.0)
        self.assertNotEqual(global_snapshot["overflow"]["content"], 0.0)

    def test_frozen_final_weight_payload_is_authoritative_and_rebuild_is_audit(self):
        fixture = self._fixture()
        frozen_references = [
            parent["final_diagnostic_application_result"]["H_pion_control_model"]
            for parent in fixture["pion_parents"]
        ]
        rebuilt_references = []
        for parent in fixture["pion_parents"]:
            final_payload = parent["final_diagnostic_application_result"]
            rebuilt = _Histogram(
                "rebuilt_audit_t{}".format(parent["t_bin_index"] + 1),
                final_payload["H_pion_control_model"].axis.edges,
            )
            rebuilt_references.append(rebuilt)
            parent["fit_result"]["_phase_a_weight_payload"] = {
                "H_pion_control_model": rebuilt,
                "weights": list(final_payload["weights"]),
                "diagnostics": {"source": "rebuilt_audit"},
            }
        evaluated_references = []
        real_evaluator = event_contract.simc_shape_pion_weight_from_value

        def evaluator(value, reference, weights):
            evaluated_references.append(reference)
            return real_evaluator(value, reference, weights)

        with mock.patch.object(
            event_contract,
            "simc_shape_pion_weight_from_value",
            side_effect=evaluator,
        ):
            contract = self._build(fixture)
        self.assertTrue(contract["available"], contract.get("reason"))
        self.assertTrue(evaluated_references)
        self.assertTrue(all(
            reference in frozen_references for reference in evaluated_references
        ))
        self.assertTrue(all(
            reference not in rebuilt_references for reference in evaluated_references
        ))
        self.assertTrue(all(
            provenance["component_rebuild_audit"]["passed"]
            for provenance in contract["baseline_weight_provenance"]
        ))

    def test_fingerprint_excludes_phi_binning(self):
        first_fixture = self._fixture()
        second_fixture = self._fixture()
        second_fixture["canonical_binning"].update({
            "phi_edges": [-180.0, -60.0, 30.0, 180.0],
            "requested_num_phi_bins": 3,
            "actual_num_phi_bins": 3,
        })
        first = self._build(first_fixture)
        second = self._build(second_fixture)
        self.assertEqual(first["contract_fingerprint"], second["contract_fingerprint"])
        fingerprint_json = json.dumps(first["fingerprint_inputs"]).lower()
        self.assertNotIn("phi_edges", fingerprint_json)
        self.assertNotIn("requested_num_phi_bins", fingerprint_json)
        self.assertNotIn("actual_num_phi_bins", fingerprint_json)

    def test_pair_id_is_serialized_and_validated_but_excluded_from_fingerprints(self):
        first_fixture = self._fixture()
        second_fixture = self._fixture()
        renamed_pair_id = "pair-phase-a-renamed"
        second_fixture["canonical_binning"]["canonical_interval_pair_id"] = (
            renamed_pair_id
        )
        second_fixture["pion_parents"] = tuple(
            dict(parent, canonical_interval_pair_id=renamed_pair_id)
            for parent in second_fixture["pion_parents"]
        )

        first = self._build(first_fixture)
        second = self._build(second_fixture)
        self.assertTrue(first["available"], first.get("reason"))
        self.assertTrue(second["available"], second.get("reason"))
        self.assertEqual(
            first["pion_event_population_fingerprint"],
            second["pion_event_population_fingerprint"],
        )
        self.assertEqual(first["contract_fingerprint"], second["contract_fingerprint"])
        self.assertEqual(
            [record["pion_parent_id"] for record in first["pion_records"]],
            [record["pion_parent_id"] for record in second["pion_records"]],
        )
        self.assertTrue(all(
            record["canonical_interval_pair_id"] == renamed_pair_id
            and record["canonical_interval_pair_hash"] == "pair-hash-phase-a"
            for record in second["pion_records"]
        ))
        fingerprint_json = json.dumps(second["fingerprint_inputs"])
        self.assertNotIn("canonical_interval_pair_id", fingerprint_json)
        self.assertIn("canonical_interval_pair_hash", fingerprint_json)

        for key in (
            "canonical_interval_pair_id", "canonical_interval_pair_hash",
        ):
            invalid_fixture = self._fixture()
            invalid_parents = list(invalid_fixture["pion_parents"])
            invalid_parents[0] = dict(
                invalid_parents[0], **{key: "wrong-{}".format(key)}
            )
            invalid_fixture["pion_parents"] = tuple(invalid_parents)
            unavailable = self._build(invalid_fixture)
            self.assertFalse(unavailable["available"])
            self.assertIn("pion_parent_{}_mismatch:t1".format(key), unavailable["reason"])

    def test_pair_hash_and_scientific_records_change_fingerprints(self):
        baseline = self._build(self._fixture())

        hash_fixture = self._fixture()
        changed_pair_hash = "pair-hash-phase-a-changed"
        hash_fixture["canonical_binning"]["canonical_interval_pair_hash"] = (
            changed_pair_hash
        )
        hash_fixture["pion_parents"] = tuple(
            dict(parent, canonical_interval_pair_hash=changed_pair_hash)
            for parent in hash_fixture["pion_parents"]
        )
        hash_changed = self._build(hash_fixture)
        self.assertNotEqual(
            baseline["pion_event_population_fingerprint"],
            hash_changed["pion_event_population_fingerprint"],
        )
        self.assertNotEqual(
            baseline["contract_fingerprint"], hash_changed["contract_fingerprint"])

        science_fixture = self._fixture()
        by_t = list(science_fixture["pion_control_cache"]["by_t"])
        records = list(by_t[0]["records"])
        records[0] = dict(records[0], P_hgcer_npeSum=99.0)
        by_t[0] = dict(by_t[0], records=tuple(records))
        science_fixture["pion_control_cache"]["by_t"] = by_t
        science_changed = self._build(science_fixture)
        self.assertTrue(science_changed["available"], science_changed.get("reason"))
        self.assertNotEqual(
            baseline["pion_event_population_fingerprint"],
            science_changed["pion_event_population_fingerprint"],
        )
        self.assertNotEqual(
            baseline["contract_fingerprint"], science_changed["contract_fingerprint"])

    def test_identity_closure_from_existing_application_passes_exact_identity(self):
        _result, application, upstream = self._accepted_application_fixture(
            production_action="bypass"
        )
        before = {
            key: event_contract.fingerprint_histogram_content_error(histogram)
            for key, histogram in upstream.items()
        }
        closure = event_contract.build_identity_host_closure_from_application(
            application, upstream
        )
        self.assertTrue(closure["passed"])
        self.assertTrue(closure["identity_transform_closure"]["passed"])
        self.assertTrue(closure["upstream_noRF_closure"]["passed"])
        self.assertTrue(closure["global_constructed_strictly_from_per_t"])
        self.assertTrue(all(
            entry["passed"]
            for entry in closure["identity_transform_closure"]["per_t"]
        ))
        self.assertTrue(
            closure["upstream_noRF_closure"]["upstream_references_unchanged"]
        )
        self.assertEqual(before, {
            key: event_contract.fingerprint_histogram_content_error(histogram)
            for key, histogram in upstream.items()
        })

    def test_identity_closure_from_existing_application_detects_nonidentity(self):
        _result, application, upstream = self._accepted_application_fixture(
            production_action="bypass"
        )
        application["final_targets"]["h_mm_nosub"].contents[1] += 1.0
        closure = event_contract.build_identity_host_closure_from_application(
            application, upstream
        )
        self.assertFalse(closure["passed"])
        self.assertFalse(closure["identity_transform_closure"]["passed"])
        self.assertFalse(
            closure["identity_transform_closure"]["global_full"]["passed"]
        )

    def test_accepted_timing_model_lambda_bypass_builds_identity_closure(self):
        result, application, upstream = self._accepted_application_fixture(
            production_action="bypass"
        )
        gate = result["diagnostics"]["lambda_preservation_gate"]
        self.assertTrue(result["accepted"])
        self.assertTrue(application["accepted"])
        self.assertEqual(gate["status"], "fail")
        self.assertEqual(gate["production_action"], "bypass")
        self.assertFalse(gate["proton_cleaning_committed"])
        self.assertNotIn("identity_host_closure", application)
        self.assertNotIn("identity_host_closure", application["diagnostics"])
        with mock.patch.object(
            event_contract,
            "_build_identity_no_proton_cleaning_application",
            side_effect=AssertionError("accepted-bypass must not use identity builder"),
        ):
            committed_host = event_contract.finalize_committed_host_application(
                result, application, upstream
            )
        closure = application["diagnostics"]["identity_host_closure"]
        self.assertTrue(application["accepted"])
        self.assertEqual(committed_host["host_state"], "identity_no_proton_cleaning")
        self.assertTrue(closure["passed"])
        self.assertTrue(closure["identity_transform_closure"]["passed"])
        self.assertTrue(closure["upstream_noRF_closure"]["passed"])

        qa = refinement_plots.proton_main_qa_payload(
            result, application, committed_host
        )
        lines = refinement_plots.proton_main_summary_lines(qa)
        for expected in (
            "overall: PASS",
            "identity-transform: PASS",
            "upstream noRF: PASS",
            "No proton subtraction was applied.",
        ):
            self.assertIn(expected, lines)
        for forbidden in (
            "overall: not recorded",
            "identity-transform: not recorded",
            "upstream noRF: not recorded",
        ):
            self.assertNotIn(forbidden, lines)

    def test_accepted_timing_model_lambda_pass_stays_proton_cleaned(self):
        result, application, upstream = self._accepted_application_fixture(
            production_action="apply"
        )
        committed_host = event_contract.finalize_committed_host_application(
            result, application, upstream
        )
        self.assertTrue(application["accepted"])
        self.assertEqual(committed_host["host_state"], "proton_cleaned")
        self.assertNotIn("identity_host_closure", application["diagnostics"])
        qa = refinement_plots.proton_main_qa_payload(
            result, application, committed_host
        )
        self.assertEqual(qa["closure_mode"], "proton_cleaned")

    def test_committed_host_finalization_is_keyed_from_final_host_state(self):
        event_source = (
            REPO_ROOT / "src" / "cuts" / "pion_hgcer_event_contract.py"
        ).read_text(encoding="utf-8")
        start = event_source.index("def finalize_committed_host_application(")
        end = event_source.index("\ndef _resolve_host_state(", start)
        finalization = event_source[start:end]
        self.assertIn(
            'if committed_host["host_state"] == "identity_no_proton_cleaning":',
            finalization,
        )
        self.assertNotIn("if not proton_cleaning_result", finalization)
        self.assertNotIn("if not proton_cleaning_application", finalization)

        rand_source = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("finalize_committed_host_application(", rand_source)
        self.assertIn("component_targets,", rand_source)

    def test_identity_no_proton_cleaning_host_closes_with_unity_factors(self):
        contract = self._build(self._fixture(identity_host=True))
        self.assertTrue(contract["available"], contract.get("reason"))
        self.assertEqual(contract["host_state"], "identity_no_proton_cleaning")
        self.assertEqual(contract["source_target_state"], "post_proton_noRF")
        self.assertFalse(contract["rf_restoration_applied"])
        self.assertTrue(contract["host_closure"]["passed"])
        self.assertTrue(all(
            record["proton_cleaning_factor"] == 1.0
            and record["final_cleaned_factor"] == 1.0
            for record in contract["kaon_host_records"]
        ))

    def test_rejected_result_builds_runtime_identity_application_and_parents(self):
        fixture = self._fixture(identity_host=True)
        # These setting-wide histograms were filled independently by the
        # fixture's ordinary pre-proton host path, before the identity builder.
        templates = fixture["proton_cleaning_application"]["final_targets"]
        template_fingerprints = {
            key: event_contract.fingerprint_histogram_content_error(histogram)
            for key, histogram in templates.items()
        }
        with mock.patch.object(
            event_contract,
            "build_identity_host_closure_from_application",
            wraps=event_contract.build_identity_host_closure_from_application,
        ) as closure_builder:
            application = event_contract._build_identity_no_proton_cleaning_application(
                proton_source_bundle=fixture["proton_source_bundle"],
                target_templates=templates,
                t_edges=fixture["canonical_binning"]["t_edges"],
                delta_edges=fixture["pion_control_cache"]["delta_edges"],
                coordinate_fingerprint="coordinate-phase-a",
                proton_cleaning_result=fixture["proton_cleaning_result"],
            )
        closure_builder.assert_called_once()
        self.assertTrue(application["accepted"])
        self.assertEqual(
            application["host_state"], "identity_no_proton_cleaning"
        )
        self.assertEqual(application["source_target_state"], "post_proton_noRF")
        self.assertFalse(application["rf_restoration_applied"])
        self.assertEqual(len(application["canonical_t_products"]), 2)
        self.assertTrue(
            application["diagnostics"]["identity_host_closure"]["passed"]
        )
        closure = application["diagnostics"]["identity_host_closure"]
        self.assertTrue(closure["identity_transform_closure"]["passed"])
        self.assertTrue(closure["upstream_noRF_closure"]["passed"])
        for selection in ("full", "cut"):
            self.assertTrue(
                closure["upstream_noRF_closure"][selection][
                    "raw_vs_upstream"
                ]["passed"]
            )
            self.assertTrue(
                closure["upstream_noRF_closure"][selection][
                    "final_vs_upstream"
                ]["passed"]
            )
        self.assertEqual(len(application["_prepared_event_weight_lookup"]), 7)
        for payload in application["_prepared_event_weight_lookup"].values():
            self.assertEqual(payload["proton_weight"], 0.0)
            self.assertEqual(payload["applied_proton_probability"], 0.0)
            self.assertEqual(payload["cleaned_factor"], 1.0)
            self.assertEqual(payload["final_cleaned_factor"], 1.0)
        strict_sum = _sum_histograms(
            "identity_runtime_strict_sum",
            [
                product["final_targets"]["h_mm_nosub"]
                for product in application["canonical_t_products"]
            ],
        )
        self.assertTrue(event_contract._histogram_closure(
            strict_sum, application["final_targets"]["h_mm_nosub"], 1.0e-10
        )["passed"])
        global_snapshot = event_contract._hist_snapshot(
            application["final_targets"]["h_mm_nosub"]
        )
        self.assertNotEqual(global_snapshot["underflow"]["content"], 0.0)
        self.assertNotEqual(global_snapshot["overflow"]["content"], 0.0)
        self.assertEqual(
            template_fingerprints,
            {
                key: event_contract.fingerprint_histogram_content_error(histogram)
                for key, histogram in templates.items()
            },
        )
        with_legacy_rf = copy.deepcopy(fixture["proton_source_bundle"])
        for source in with_legacy_rf["prepared_sources"].values():
            for entry in source["entries"].values():
                entry["rf_accept"] = False
        legacy_rf_application = (
            event_contract._build_identity_no_proton_cleaning_application(
                proton_source_bundle=with_legacy_rf,
                target_templates=templates,
                t_edges=fixture["canonical_binning"]["t_edges"],
                delta_edges=fixture["pion_control_cache"]["delta_edges"],
                coordinate_fingerprint="coordinate-phase-a",
                proton_cleaning_result=fixture["proton_cleaning_result"],
            )
        )
        self.assertEqual(
            event_contract.fingerprint_histogram_content_error(
                application["final_targets"]["h_mm_nosub"]
            ),
            event_contract.fingerprint_histogram_content_error(
                legacy_rf_application["final_targets"]["h_mm_nosub"]
            ),
        )

        fixture["proton_cleaning_application"] = application
        contract = self._build(fixture)
        self.assertTrue(contract["available"], contract.get("reason"))
        self.assertTrue(contract["host_closure"]["passed"])
        self.assertEqual(contract["host_state"], "identity_no_proton_cleaning")

    def test_identity_upstream_closure_detects_independent_corruption(self):
        fixture = self._fixture(identity_host=True)
        templates = fixture["proton_cleaning_application"]["final_targets"]

        def build(source_bundle, references):
            return event_contract._build_identity_no_proton_cleaning_application(
                proton_source_bundle=source_bundle,
                target_templates=references,
                t_edges=fixture["canonical_binning"]["t_edges"],
                delta_edges=fixture["pion_control_cache"]["delta_edges"],
                coordinate_fingerprint="coordinate-phase-a",
                proton_cleaning_result=fixture["proton_cleaning_result"],
            )

        missing_prepared_event = copy.deepcopy(fixture["proton_source_bundle"])
        missing_prepared_event["prepared_sources"]["prompt"]["entries"].pop(0)
        with self.assertRaises(event_contract.EventContractUnavailable) as caught:
            build(missing_prepared_event, templates)
        reason = str(caught.exception)
        self.assertIn("identity_host_upstream_noRF_closure_failed", reason)
        for comparison in (
            "full/raw_vs_upstream", "full/final_vs_upstream",
            "cut/raw_vs_upstream", "cut/final_vs_upstream",
        ):
            self.assertIn(comparison, reason)

        corruption_cases = (
            ("regular_content", "h_mm_nosub", 1, "content"),
            ("error_only", "h_mm_nosub", 1, "error"),
            ("underflow", "h_mm_nosub", 0, "content"),
            ("overflow", "h_mm", -1, "content"),
        )
        for label, target_key, bin_index, corruption_kind in corruption_cases:
            with self.subTest(label=label):
                corrupted = {
                    key: histogram.Clone("{}_{}".format(key, label))
                    for key, histogram in templates.items()
                }
                if corruption_kind == "error":
                    corrupted[target_key].sumw2[bin_index] += 1.0
                else:
                    corrupted[target_key].contents[bin_index] += 1.0
                before = {
                    key: event_contract.fingerprint_histogram_content_error(histogram)
                    for key, histogram in corrupted.items()
                }
                with self.assertRaisesRegex(
                    event_contract.EventContractUnavailable,
                    "identity_host_upstream_noRF_closure_failed",
                ):
                    build(fixture["proton_source_bundle"], corrupted)
                self.assertEqual(before, {
                    key: event_contract.fingerprint_histogram_content_error(histogram)
                    for key, histogram in corrupted.items()
                })

    def test_lambda_gate_committed_state_controls_host_classification(self):
        bypass = self._fixture(identity_host=True)
        bypass["proton_cleaning_result"]["accepted"] = True
        bypass["proton_cleaning_application"]["accepted"] = True
        bypass["proton_cleaning_application"].pop("host_state")
        bypass_contract = self._build(bypass)
        self.assertTrue(bypass_contract["available"], bypass_contract.get("reason"))
        self.assertEqual(
            bypass_contract["host_state"], "identity_no_proton_cleaning"
        )
        self.assertFalse(bypass_contract["proton_cleaning_committed"])
        self.assertEqual(
            bypass_contract["lambda_gate_production_action"], "bypass"
        )

        committed = self._build(self._fixture())
        self.assertTrue(committed["available"], committed.get("reason"))
        self.assertEqual(committed["host_state"], "proton_cleaned")
        self.assertTrue(committed["proton_cleaning_committed"])
        self.assertEqual(committed["lambda_gate_production_action"], "apply")
        self.assertTrue(any(
            record["final_cleaned_factor"] not in (0.0, 1.0)
            for record in committed["kaon_host_records"]
        ))

        contradictory = self._fixture()
        contradictory["proton_cleaning_application"]["host_state"] = (
            "identity_no_proton_cleaning"
        )
        unavailable = self._build(contradictory)
        self.assertFalse(unavailable["available"])
        self.assertIn("committed_state_contradiction", unavailable["reason"])

    def test_rejected_cleaning_without_explicit_identity_provenance_is_unavailable(self):
        fixture = self._fixture(identity_host=True)
        fixture["proton_cleaning_application"].pop("host_state")
        fixture["proton_cleaning_application"].pop(
            "_prepared_event_weight_lookup"
        )
        contract = self._build(fixture)
        self.assertFalse(contract["available"])
        self.assertEqual(
            contract["reason"],
            "identity_no_proton_cleaning_host_provenance_unavailable",
        )

    def test_noRF_sources_are_required_and_rf_accept_is_not_contract_semantics(self):
        fixture = self._fixture()
        contract = self._build(fixture)
        self.assertTrue(contract["available"], contract.get("reason"))
        self.assertNotIn("rf_accept", json.dumps(contract))

        invalid = self._fixture()
        invalid["proton_source_bundle"]["prepared_sources"]["prompt"][
            "tree_name"
        ] = "Cut_Kaon_prompt_RF"
        unavailable = self._build(invalid)
        self.assertFalse(unavailable["available"])
        self.assertIn("proton_host_source_not_noRF", unavailable["reason"])

        rf_restored = self._fixture()
        rf_restored["proton_cleaning_application"]["diagnostics"][
            "rf_applied"
        ] = True
        rf_restored["proton_cleaning_result"]["diagnostics"][
            "rf_applied"
        ] = True
        unavailable = self._build(rf_restored)
        self.assertFalse(unavailable["available"])
        self.assertIn("rf_restoration_not_explicitly_disabled", unavailable["reason"])

        contradictory_application = self._fixture()
        contradictory_application["proton_cleaning_application"][
            "rf_restoration_applied"
        ] = True
        unavailable = self._build(contradictory_application)
        self.assertFalse(unavailable["available"])
        self.assertIn("rf_restoration_not_explicitly_disabled", unavailable["reason"])

        mismatched_alias = self._fixture()
        mismatched_alias["proton_cleaning_result"][
            "_prepared_event_weight_lookup"
        ]["prompt:0"]["final_cleaned_factor"] = 0.5
        unavailable = self._build(mismatched_alias)
        self.assertFalse(unavailable["available"])
        self.assertIn("noRF_factor_alias_mismatch", unavailable["reason"])

    def test_absolute_support_effective_entries_preserves_signed_sum(self):
        metrics = event_contract._new_metrics()
        for contribution in (2.0, -1.0):
            event_contract._record_metrics(metrics, {
                "allcuts": True,
                "nommcuts": True,
                "signed_source_coefficient": 1.0,
                "signed_baseline_event_contribution": contribution,
            })
        result = event_contract._finalize_metrics(metrics)
        self.assertEqual(result["signed_weighted_sum"], 1.0)
        self.assertEqual(result["absolute_weighted_support"], 3.0)
        self.assertEqual(result["sumw2"], 5.0)
        self.assertEqual(result["effective_entries"], 9.0 / 5.0)

    def test_signed_prompt_random_dummy_and_dummy_random_algebra(self):
        contract = self._build(self._fixture())
        contributions = {
            record["source_label"]: record["signed_baseline_event_contribution"]
            for record in contract["pion_records"]
        }
        self.assertEqual(contributions, {
            "prompt": 4.0,
            "rand": -1.5,
            "dummy": -4.0,
            "dummy_rand": 1.25,
        })
        identity = self._build(self._fixture(identity_host=True))
        host_contributions = {}
        for record in identity["kaon_host_records"]:
            host_contributions.setdefault(record["source_label"], 0.0)
            host_contributions[record["source_label"]] += record[
                "signed_host_event_contribution"
            ]
        self.assertEqual(host_contributions, {
            "dummy": -2.0,
            "dummy_rand": 0.25,
            "prompt": 6.0,
            "rand": -0.5,
        })

    def test_phase_a2_has_no_method_a_or_refinement_fields(self):
        contract = self._build(self._fixture(identity_host=True))
        keys = set()

        def collect(value):
            if isinstance(value, dict):
                keys.update(value)
                for child in value.values():
                    collect(child)
            elif isinstance(value, (list, tuple)):
                for child in value:
                    collect(child)

        collect(contract)
        self.assertTrue({
            "L_A", "L_B", "L", "C", "C_final", "refined_pion_weight"
        }.isdisjoint(keys))
        self.assertFalse(contract["refinement_applied"])

    def test_exact_host_t_and_delta_geometry_assignment_matrix(self):
        mutations = (
            ("prompt:0", "t_index", 1, "proton_event_t_assignment_mismatch"),
            ("prompt:0", "t_index", None, "proton_event_t_assignment_mismatch"),
            ("prompt:5", "t_index", 0, "proton_event_t_assignment_mismatch"),
            ("prompt:0", "delta_index", 1, "proton_event_delta_assignment_mismatch"),
            ("prompt:0", "delta_index", None, "proton_event_delta_assignment_mismatch"),
            ("rand:1", "delta_index", 0, "proton_event_delta_assignment_mismatch"),
        )
        for signature, key, value, reason in mutations:
            with self.subTest(signature=signature, key=key, value=value):
                fixture = self._fixture()
                fixture["proton_cleaning_result"][
                    "_prepared_event_weight_lookup"
                ][signature][key] = value
                contract = self._build(fixture)
                self.assertFalse(contract["available"])
                self.assertIn(reason, contract["reason"])
        unchanged = self._build(self._fixture())
        self.assertTrue(unchanged["available"], unchanged.get("reason"))

    def test_exact_pion_cached_t_assignment_rejects_none_for_valid_event(self):
        fixture = self._fixture()
        fixture["pion_control_cache"]["by_t"][0]["records"][0][
            "t_index"
        ] = None
        contract = self._build(fixture)
        self.assertFalse(contract["available"])
        self.assertIn("pion_event_cached_t_assignment_mismatch", contract["reason"])

    def test_existing_zero_payload_is_reused_but_skip_is_unavailable(self):
        zero_fixture = self._fixture()
        for parent in zero_fixture["pion_parents"]:
            payload = parent["final_diagnostic_application_result"]
            reference = payload["H_pion_control_model"]
            payload["weights"] = [0.0] * (reference.GetNbinsX() + 2)
            payload["H_pion_subtraction_template_MM"] = _Histogram(
                "zero_cut_t{}".format(parent["t_bin_index"] + 1)
            )
            payload["H_pion_subtraction_template_MM_nosub"] = _Histogram(
                "zero_full_t{}".format(parent["t_bin_index"] + 1)
            )
            parent["final_diagnostic_application_status"]["final_status"] = "zero"
        zero = self._build(zero_fixture, actions={0: "zero", 1: "zero"})
        self.assertTrue(zero["available"])
        self.assertIsNone(zero["pion_closure"]["named_canonical_global"])
        skipped = self._build(
            self._fixture(), actions={0: "component_weight", 1: "skip_bin"}
        )
        self.assertFalse(skipped["available"])
        self.assertIn("baseline_parent_policy_unavailable", skipped["reason"])

    def test_existing_single_scale_payload_is_reused_without_refitting(self):
        fixture = self._fixture()
        for t_index, parent in enumerate(fixture["pion_parents"]):
            payload = parent["final_diagnostic_application_result"]
            reference = payload["H_pion_control_model"]
            scale = 1.5 + t_index
            scalar_weights = [scale] * (reference.GetNbinsX() + 2)
            records = fixture["pion_control_cache"]["by_t"][t_index]["records"]
            payload["weights"] = scalar_weights
            payload["particle_subtraction_effective_scale"] = scale
            payload["H_pion_subtraction_template_MM"] = _fill_template(
                "single_cut_t{}".format(t_index + 1), records, reference,
                scalar_weights, "allcuts",
            )
            payload["H_pion_subtraction_template_MM_nosub"] = _fill_template(
                "single_full_t{}".format(t_index + 1), records, reference,
                scalar_weights, "nommcuts",
            )
            parent["final_diagnostic_application_status"][
                "final_status"
            ] = "applied_fallback"
        contract = self._build(
            fixture, actions={0: "single_scale", 1: "single_scale"}
        )
        self.assertTrue(contract["available"], contract.get("reason"))
        self.assertEqual(
            [entry["application_action"] for entry in contract["pion_closure"]["per_t"]],
            ["single_scale", "single_scale"],
        )
        self.assertIsNone(contract["pion_closure"]["named_canonical_global"])

    def test_provenance_mismatch_is_unavailable_and_summary_is_lightweight(self):
        fixture = self._fixture()
        parents = list(fixture["pion_parents"])
        parents[0] = dict(parents[0], coordinate_fingerprint="wrong-coordinate")
        fixture["pion_parents"] = tuple(parents)
        contract = self._build(fixture)
        self.assertFalse(contract["available"])
        self.assertIn("coordinate_fingerprint_mismatch", contract["reason"])
        summary = event_contract.summarize_pion_hgcer_event_contract(contract)
        self.assertFalse(summary["available"])
        self.assertEqual(
            summary["fingerprint_schema_version"],
            event_contract.EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION,
        )
        self.assertNotIn("pion_records", summary)
        json.dumps(summary, allow_nan=False)

    def test_rand_sub_hook_is_after_parent_freeze_and_phase_a_has_no_renderer(self):
        rand_source = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(
            encoding="utf-8"
        )
        parent_call = rand_source.index("build_setting_t_bin_pion_parents(")
        contract_call = rand_source.index(
            "build_pion_hgcer_event_contract(", parent_call
        )
        identity_call = rand_source.index(
            "_build_identity_no_proton_cleaning_application("
        )
        cache_call = rand_source.index(
            "_build_authoritative_pion_control_source_cache(", identity_call
        )
        self.assertLess(identity_call, cache_call)
        self.assertGreater(contract_call, parent_call)
        self.assertIn('"P_hgcer_xAtCer"', rand_source)
        self.assertIn('"P_hgcer_yAtCer"', rand_source)
        module_source = (
            REPO_ROOT / "src" / "cuts" / "pion_hgcer_event_contract.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("TCanvas", module_source)
        self.assertNotIn(".Print(", module_source)
        self.assertNotIn("write_", module_source)

    def test_phase_a_diff_allowlist_and_rand_sub_patch_are_small(self):
        baseline = "af4fe473e"
        phase_a_head = "57e366806"
        changed = subprocess.run(
            [
                "git", "-c", "core.safecrlf=false", "diff", "--name-only",
                baseline, phase_a_head, "--",
            ],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
        self.assertEqual(changed, [
            "src/cuts/pion_hgcer_event_contract.py",
            "src/cuts/rand_sub.py",
            "testing/test_pion_hgcer_event_contract.py",
        ])
        numstat = subprocess.run(
            [
                "git", "-c", "core.safecrlf=false", "diff", "--numstat",
                baseline, phase_a_head, "--", "src/cuts/rand_sub.py",
            ],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip().split()
        self.assertEqual(numstat[-1], "src/cuts/rand_sub.py")
        self.assertLess(int(numstat[0]), 200)
        self.assertLess(int(numstat[1]), 30)
        phase_a3_changed = subprocess.run(
            [
                "git", "-c", "core.safecrlf=false", "diff", "--name-only",
                "bd8d9e9444bd46ddc3734346f92105347e658cdb", phase_a_head, "--",
            ],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.splitlines()
        self.assertEqual(phase_a3_changed, [
            "src/cuts/pion_hgcer_event_contract.py",
            "testing/test_pion_hgcer_event_contract.py",
        ])


@unittest.skipUnless(ROOT is not None, "PyROOT is not available")
class PhaseAEventContractRootTests(unittest.TestCase):
    def test_real_th1_sumw2_underflow_overflow_closure(self):
        authoritative = ROOT.TH1D("phase_a_root_authoritative", "", 2, 0.0, 2.0)
        authoritative.SetDirectory(0)
        authoritative.Sumw2()
        for value, weight in ((-0.2, 1.5), (0.5, -0.25), (2.2, 0.75)):
            authoritative.Fill(value, weight)
        reconstructed = authoritative.Clone("phase_a_root_reconstructed")
        reconstructed.SetDirectory(0)
        closure = event_contract._histogram_closure(
            reconstructed, authoritative, 1.0e-10
        )
        self.assertTrue(closure["passed"])
        self.assertEqual(
            closure["reconstructed"]["underflow"]["content"], 1.5
        )
        self.assertEqual(
            closure["reconstructed"]["overflow"]["content"], 0.75
        )

    def test_identity_builder_closes_against_independent_root_upstream(self):
        templates = {}
        for key in ("h_mm", "h_mm_nosub"):
            histogram = ROOT.TH1D(
                "phase_a3_root_upstream_{}".format(key), "", 2, 0.0, 2.0
            )
            histogram.SetDirectory(0)
            histogram.Sumw2()
            templates[key] = histogram
        rows = (
            ("prompt", 0, -0.2, 0.3, 2.0),
            ("rand", 1, 0.5, 0.7, -0.5),
            ("dummy", 2, 2.2, 1.4, -1.0),
            ("dummy_rand", 3, 1.5, 1.8, 0.25),
        )
        prepared_sources = {}
        for label, entry_index, mm_value, t_value, coefficient in rows:
            prepared_sources[label] = {
                "coefficient": coefficient,
                "tree_name": "Cut_Kaon_{}_noRF".format(label),
                "entries": {
                    entry_index: {
                        "adj_mm": mm_value,
                        "adj_t": t_value,
                        "delta_value": 0.0,
                        "allcuts": True,
                        "nommcuts": True,
                    }
                },
            }
            templates["h_mm"].Fill(mm_value, coefficient)
            templates["h_mm_nosub"].Fill(mm_value, coefficient)
        before = {
            key: event_contract.fingerprint_histogram_content_error(histogram)
            for key, histogram in templates.items()
        }
        application = event_contract._build_identity_no_proton_cleaning_application(
            proton_source_bundle={"prepared_sources": prepared_sources},
            target_templates=templates,
            t_edges=[0.0, 1.0, 2.0],
            delta_edges=[-10.0, 0.0, 10.0],
            coordinate_fingerprint="phase-a3-root-coordinate",
            proton_cleaning_result={
                "accepted": False,
                "diagnostics": {"rf_applied": False},
            },
        )
        upstream_closure = application["diagnostics"]["identity_host_closure"][
            "upstream_noRF_closure"
        ]
        self.assertTrue(upstream_closure["passed"])
        self.assertTrue(upstream_closure["upstream_references_unchanged"])
        self.assertEqual(before, {
            key: event_contract.fingerprint_histogram_content_error(histogram)
            for key, histogram in templates.items()
        })


if __name__ == "__main__":
    unittest.main()
