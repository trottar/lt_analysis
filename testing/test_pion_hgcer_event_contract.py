"""Focused Phase-A event-contract tests with optional PyROOT closure."""

from __future__ import annotations

import copy
import json
import math
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
        if record[selection] and record["rf_accept"]:
            contribution = record["coefficient"] * record["final_cleaned_factor"]
            if contribution:
                histogram.Fill(record["adj_mm"], contribution)
    return histogram


class PhaseAEventContractTests(unittest.TestCase):
    def _fixture(self):
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
                    "delta_index": 0, "P_hgcer_npeSum": 0.5,
                    "P_hgcer_xAtCer": 1.25, "P_hgcer_yAtCer": -2.5,
                    "allcuts": True, "nommcuts": True,
                },
                {
                    "source_label": "rand", "entry_index": 1,
                    "coefficient": -0.5, "source_tree_name": "Cut_Pion_rand_noRF",
                    "rf_state": "noRF", "raw_MM": 1.25, "adj_MM": 1.25,
                    "raw_t": 0.35, "adj_t": 0.35, "ssdelta": 12.0,
                    "delta_index": None, "P_hgcer_npeSum": 1.0,
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
                    "delta_index": 1, "P_hgcer_npeSum": 0.25,
                    "P_hgcer_xAtCer": -3.0, "P_hgcer_yAtCer": 4.0,
                    "allcuts": True, "nommcuts": True,
                },
                {
                    "source_label": "dummy_rand", "entry_index": 3,
                    "coefficient": 0.25, "source_tree_name": "Cut_Pion_rand_noRF",
                    "rf_state": "noRF", "raw_MM": 1.75, "adj_MM": 1.75,
                    "raw_t": 1.7, "adj_t": 1.7, "ssdelta": 0.0,
                    "delta_index": 1, "P_hgcer_npeSum": 1.5,
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
                    "cleaned_factor": 0.6, "final_cleaned_factor": 0.5,
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
                    "cleaned_factor": 0.4, "final_cleaned_factor": 0.0,
                    "rf_accept": False,
                },
            ],
        }
        coefficients = {"prompt": 2.0, "rand": -0.5, "dummy": -1.0}
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
            })
        proton_application = {
            "coordinate_fingerprint": coordinate,
            "canonical_t_products": tuple(proton_products),
            "final_targets": {
                "h_mm_nosub": _sum_histograms("host_global_full", host_full),
                "h_mm": _sum_histograms("host_global_cut", host_cut),
            },
        }
        proton_result = {
            "accepted": True,
            "method": "timing_t_event_weight",
            "coordinate_fingerprint": coordinate,
            "_prepared_event_weight_lookup": lookup,
            "diagnostics": {
                "event_weight_source": "setting_wide_immutable_prepared_lookup"
            },
        }
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
        contract = self._build(fixture)
        self.assertTrue(contract["available"], contract.get("reason"))
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
        json.dumps(contract, allow_nan=False)

    def test_host_uses_existing_final_factor_and_no_pion_refinement(self):
        contract = self._build(self._fixture())
        prompt = next(
            record for record in contract["kaon_host_records"]
            if record["source_label"] == "prompt" and record["analysis_MM"] == 0.5
        )
        self.assertEqual(prompt["proton_cleaning_factor"], 0.6)
        self.assertEqual(prompt["final_cleaned_factor"], 0.5)
        self.assertEqual(prompt["signed_host_event_contribution"], 1.0)
        self.assertIsNone(prompt["pion_refinement_factor"])
        self.assertEqual(len(contract["host_records_outside_geometry"]), 1)
        global_snapshot = contract["host_closure"]["global_full"]["reconstructed"]
        self.assertNotEqual(global_snapshot["underflow"]["content"], 0.0)
        self.assertEqual(global_snapshot["overflow"]["content"], 0.0)

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
        self.assertGreater(contract_call, parent_call)
        self.assertIn('"P_hgcer_xAtCer"', rand_source)
        self.assertIn('"P_hgcer_yAtCer"', rand_source)
        module_source = (
            REPO_ROOT / "src" / "cuts" / "pion_hgcer_event_contract.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("TCanvas", module_source)
        self.assertNotIn(".Print(", module_source)
        self.assertNotIn("write_", module_source)


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


if __name__ == "__main__":
    unittest.main()
