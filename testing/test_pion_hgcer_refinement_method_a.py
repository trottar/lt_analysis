"""Pure-Python contracts for Phase-B Method-A HGCer diagnostics."""

from __future__ import annotations

from copy import deepcopy
import json
import math
from pathlib import Path
import sys
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

import pion_hgcer_refinement_method_a as method_a


COORDINATE_FINGERPRINT = "coordinate-fingerprint"
_DEFAULT_COORDINATE = object()


def _phase_a(t_edges=(0.0, 1.0), delta_edges=(0.0, 10.0)):
    return {
        "schema_version": "pion_hgcer_event_contract/v1",
        "status": "available",
        "available": True,
        "contract_fingerprint": "phase-a-fingerprint",
        "coordinate_fingerprint": COORDINATE_FINGERPRINT,
        "physical_pion_control_mask_fingerprint": "physical-control-mask",
        "canonical_t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "pion_records": [{"baseline_weight": 0.75}],
        "kaon_host_records": [{"final_cleaned_factor": 1.0}],
    }


def _provenance():
    specifications = {
        "prompt": ("prompt", 1.0),
        "rand": ("random", -1.0),
        "dummy": ("dummy_prompt", -1.0),
        "dummy_rand": ("dummy_random", 1.0),
    }
    return {
        "pion": {
            label: {
                "source_label": label,
                "source_role": role,
                "tree_name": "Cut_Pion_Events_{}_noRF".format(label),
                "rf_state": "noRF",
                "coefficient": coefficient,
                "proton_factor_scope": "none",
            }
            for label, (role, coefficient) in specifications.items()
        }
    }


def _record(
    entry_index,
    npe,
    *,
    source="prompt",
    weight=None,
    analysis_t=0.5,
    t_index=0,
    delta=5.0,
    delta_index=0,
    allcuts=True,
    nommcuts=True,
    x=_DEFAULT_COORDINATE,
    y=_DEFAULT_COORDINATE,
):
    default_coefficient = {"prompt": 1.0, "rand": -1.0, "dummy": -1.0, "dummy_rand": 1.0}[source]
    coefficient = default_coefficient if weight is None else float(weight)
    return {
        "side": "pion",
        "source_label": source,
        "entry_index": int(entry_index),
        "coefficient": coefficient,
        "proton_cleaning_factor": None,
        "diagnostic_weight": coefficient,
        "analysis_t": float(analysis_t),
        "coordinate_fingerprint": COORDINATE_FINGERPRINT,
        "canonical_t_index": t_index,
        "ssdelta": float(delta),
        "delta_index": delta_index,
        "P_hgcer_npeSum": float(npe),
        "P_hgcer_xAtCer": (
            float(entry_index) if x is _DEFAULT_COORDINATE else x
        ),
        "P_hgcer_yAtCer": (
            float(-entry_index) if y is _DEFAULT_COORDINATE else y
        ),
        "allcuts": bool(allcuts),
        "nommcuts": bool(nommcuts),
        "rf_applied_to_diagnostic": False,
        "phi": 123.0,
    }


def _diagnostic(records, t_edges=(0.0, 1.0), delta_edges=(0.0, 10.0)):
    provenance = _provenance()
    for record in records:
        provenance["pion"][record["source_label"]]["coefficient"] = record[
            "diagnostic_weight"
        ]
    return {
        "status": "available",
        "rf_restoration_applied": False,
        "coordinate_fingerprint": COORDINATE_FINGERPRINT,
        "t_edges": list(t_edges),
        "delta_edges": list(delta_edges),
        "source_provenance": provenance,
        "records": {"pion": tuple(records), "kaon": ()},
    }


def _supported_records(*, low=10, control=20, delta=5.0, delta_index=0, t=0.5, t_index=0):
    records = []
    for index in range(low):
        records.append(
            _record(
                index,
                0.25 + 1.75 * index / max(1, low - 1),
                analysis_t=t,
                t_index=t_index,
                delta=delta,
                delta_index=delta_index,
            )
        )
    for offset in range(control):
        records.append(
            _record(
                low + offset,
                2.01 + offset / 10.0,
                analysis_t=t,
                t_index=t_index,
                delta=delta,
                delta_index=delta_index,
            )
        )
    return records


class PionHGCerMethodATests(unittest.TestCase):
    def test_population_boundary_partition_observables_and_position_statistics(self):
        records = _supported_records(low=10, control=20)
        records.extend(
            [
                _record(100, 0.0),
                _record(101, -0.1),
                _record(102, 0.01, nommcuts=False, allcuts=False),
            ]
        )
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(records), _phase_a()
        )
        self.assertTrue(result["available"])
        cell = result["cells"][0]
        self.assertEqual(cell["support_class"], "supported")
        self.assertEqual(cell["prompt_positive_count"], 30)
        self.assertEqual(cell["prompt_low_count"], 10)
        self.assertEqual(cell["prompt_control_count"], 20)
        self.assertTrue(cell["partition_closure_passed"])
        self.assertAlmostEqual(cell["f_low"], 1.0 / 3.0)
        self.assertAlmostEqual(cell["R_low_control"], 0.5)
        self.assertLessEqual(cell["minimum_npe"], 0.25)
        self.assertGreater(cell["maximum_npe"], 2.0)
        self.assertIsNotNone(cell["hgcer_x_mean"])
        self.assertIsNotNone(cell["hgcer_y_rms"])
        self.assertEqual(result["summary"]["records_nonpositive_response"], 2)
        self.assertEqual(result["summary"]["records_not_nommcuts_or_allcuts"], 1)

    def test_wilson_interval_and_type7_quantiles_are_bounded_and_deterministic(self):
        low, high = method_a._wilson_interval(1, 2)
        self.assertAlmostEqual(low, 0.09453120573423074)
        self.assertAlmostEqual(high, 0.9054687942657693)
        self.assertLessEqual(0.0, low)
        self.assertLessEqual(low, 0.5)
        self.assertLessEqual(0.5, high)
        self.assertLessEqual(high, 1.0)
        self.assertAlmostEqual(method_a._type7_quantile([0.0, 10.0, 20.0, 30.0], 0.25), 7.5)
        self.assertAlmostEqual(method_a._type7_quantile([0.0, 10.0, 20.0, 30.0], 0.50), 15.0)
        self.assertLess(
            method_a._ratio_from_fraction(low), method_a._ratio_from_fraction(high)
        )
        for successes in (0, 1, 25):
            interval_low, interval_high = method_a._wilson_interval(successes, 25)
            fraction = float(successes) / 25.0
            self.assertLessEqual(0.0, interval_low)
            self.assertLessEqual(interval_low, fraction)
            self.assertLessEqual(fraction, interval_high)
            self.assertLessEqual(interval_high, 1.0)
        self.assertGreater(method_a._wilson_interval(0, 25)[1], 0.0)

    def test_support_thresholds_zero_low_and_zero_control(self):
        zero_low = [_record(index, 3.0) for index in range(25)]
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(zero_low), _phase_a()
        )
        cell = result["cells"][0]
        self.assertEqual(cell["support_class"], "marginal")
        self.assertEqual(cell["f_low"], 0.0)
        self.assertEqual(cell["R_low_control"], 0.0)

        zero_control = [_record(index, 1.0) for index in range(30)]
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(zero_control), _phase_a()
        )
        cell = result["cells"][0]
        self.assertEqual(cell["support_class"], "unsupported")
        self.assertEqual(cell["method_A_status"], "unavailable")
        self.assertIsNone(cell["f_low"])
        self.assertIsNone(cell["R_low_control"])
        self.assertIsNone(cell["R_low_control_high"])

    def test_signed_support_uses_absolute_neff_and_rejects_cancelling_denominator(self):
        records = _supported_records(low=5, control=25)
        for index in range(25):
            records.append(
                _record(100 + index, 3.0, source="rand", weight=-1.0)
            )
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(records), _phase_a()
        )
        cell = result["cells"][0]
        expected_neff = cell["signed_positive_abs_support"] ** 2 / cell[
            "signed_positive_sumw2"
        ]
        self.assertAlmostEqual(cell["signed_positive_neff"], expected_neff)
        self.assertEqual(cell["signed_control_yield"], 0.0)
        self.assertEqual(cell["signed_control_abs_support"], 50.0)
        self.assertIsNone(cell["signed_R_low_control"])
        self.assertEqual(cell["prompt_vs_signed_status"], "signed_unavailable")

    def test_signed_ratio_requires_low_neff_and_non_cancelling_positive_yield(self):
        sparse_low = _supported_records(low=5, control=25)
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(sparse_low), _phase_a()
        )
        cell = result["cells"][0]
        self.assertEqual(cell["support_class"], "supported")
        self.assertEqual(cell["signed_low_neff"], 5.0)
        self.assertEqual(cell["signed_low_yield"], 5.0)
        self.assertIsNone(cell["signed_R_low_control"])
        self.assertEqual(cell["prompt_vs_signed_status"], "signed_unavailable")

        cancelling_low = _supported_records(low=10, control=25)
        cancelling_low.extend(
            _record(100 + index, 1.0, source="rand", weight=-1.0)
            for index in range(9)
        )
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(cancelling_low), _phase_a()
        )
        cell = result["cells"][0]
        self.assertGreaterEqual(cell["signed_low_neff"], 10.0)
        self.assertEqual(cell["signed_low_yield"], 1.0)
        self.assertEqual(cell["signed_low_abs_support"], 19.0)
        self.assertIsNone(cell["signed_R_low_control"])
        self.assertEqual(cell["prompt_vs_signed_status"], "signed_unavailable")

    def test_prompt_random_dummy_and_dummy_random_algebra_is_applied_once(self):
        records = _supported_records(low=10, control=20)
        records.extend(
            [
                _record(100, 1.0, source="rand", weight=-0.5),
                _record(101, 3.0, source="dummy", weight=-0.25),
                _record(102, 3.0, source="dummy_rand", weight=0.125),
            ]
        )
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(records), _phase_a()
        )
        cell = result["cells"][0]
        self.assertAlmostEqual(cell["signed_positive_yield"], 29.375)
        self.assertAlmostEqual(cell["signed_low_yield"], 9.5)
        self.assertAlmostEqual(cell["signed_control_yield"], 19.875)
        self.assertAlmostEqual(cell["signed_positive_abs_support"], 30.875)
        self.assertAlmostEqual(cell["signed_positive_sumw2"], 30.328125)
        self.assertIsNotNone(cell["signed_R_low_control"])

    def test_prompt_vs_signed_statuses_are_deterministic(self):
        records = [
            {"npe": 1.0, "x": 0.0, "y": 0.0} for _ in range(10)
        ] + [{"npe": 3.0, "x": 0.0, "y": 0.0} for _ in range(30)]
        prompt = method_a._prompt_metrics(records, method_a._resolved_config(None))
        interval_width = prompt["ratio_high"] - prompt["ratio"]
        self.assertEqual(
            method_a._prompt_vs_signed_status(
                prompt, prompt["ratio"], 0.0, method_a._resolved_config(None)
            ),
            "consistent",
        )
        self.assertEqual(
            method_a._prompt_vs_signed_status(
                prompt,
                prompt["ratio_high"] + interval_width,
                1.5 * interval_width,
                method_a._resolved_config(None),
            ),
            "marginal",
        )
        self.assertEqual(
            method_a._prompt_vs_signed_status(
                prompt,
                prompt["ratio_high"] + 2.0 * interval_width,
                0.5 * interval_width,
                method_a._resolved_config(None),
            ),
            "inconsistent",
        )
        self.assertEqual(
            method_a._prompt_vs_signed_status(
                prompt, None, None, method_a._resolved_config(None)
            ),
            "signed_unavailable",
        )

    def test_zero_low_prompt_uses_wilson_interval_for_signed_comparison(self):
        records = [_record(index, 3.0) for index in range(25)]
        records.extend(
            _record(100 + index, 1.0, source="dummy_rand", weight=0.1)
            for index in range(10)
        )
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(records), _phase_a()
        )
        cell = result["cells"][0]
        self.assertEqual(cell["support_class"], "marginal")
        self.assertEqual(cell["method_A_status"], "available")
        self.assertEqual(cell["f_low"], 0.0)
        self.assertGreater(cell["f_low_high"], 0.0)
        self.assertAlmostEqual(cell["signed_R_low_control"], 0.04)
        self.assertLessEqual(
            cell["signed_R_low_control"], cell["R_low_control_high"]
        )
        self.assertEqual(cell["prompt_vs_signed_status"], "consistent")

    def test_missing_xy_preserves_core_response_and_splits_geometry_summaries(self):
        finite_records = _supported_records(low=10, control=20)
        missing_records = deepcopy(finite_records)
        missing_records[0]["P_hgcer_xAtCer"] = None
        missing_records[1]["P_hgcer_yAtCer"] = float("nan")

        finite = method_a.build_pion_hgcer_method_a(
            _diagnostic(finite_records), _phase_a()
        )
        missing = method_a.build_pion_hgcer_method_a(
            _diagnostic(missing_records), _phase_a()
        )
        self.assertTrue(missing["available"])
        finite_cell = finite["cells"][0]
        missing_cell = missing["cells"][0]
        core_fields = (
            "prompt_positive_count",
            "prompt_low_count",
            "prompt_control_count",
            "f_low",
            "R_low_control",
            "support_class",
            "method_A_status",
        )
        for field in core_fields:
            self.assertEqual(finite_cell[field], missing_cell[field])
        self.assertEqual(missing_cell["hgcer_x_valid_count"], 29)
        self.assertEqual(missing_cell["hgcer_y_valid_count"], 29)
        self.assertEqual(missing_cell["hgcer_xy_valid_count"], 28)
        self.assertEqual(missing_cell["hgcer_xy_missing_count"], 2)
        self.assertEqual(missing_cell["low_hgcer_xy_valid_count"], 8)
        self.assertEqual(missing_cell["control_hgcer_xy_valid_count"], 20)

        low_x = method_a._distribution(range(1, 10))
        low_y = method_a._distribution([0] + list(range(-2, -10, -1)))
        control_x = method_a._distribution(range(10, 30))
        control_y = method_a._distribution(range(-10, -30, -1))
        for prefix, expected_x, expected_y in (
            ("low", low_x, low_y),
            ("control", control_x, control_y),
        ):
            self.assertAlmostEqual(
                missing_cell["{}_hgcer_x_mean".format(prefix)], expected_x["mean"]
            )
            self.assertAlmostEqual(
                missing_cell["{}_hgcer_x_rms".format(prefix)], expected_x["rms"]
            )
            self.assertAlmostEqual(
                missing_cell["{}_hgcer_x_median".format(prefix)], expected_x["q50"]
            )
            self.assertAlmostEqual(
                missing_cell["{}_hgcer_y_mean".format(prefix)], expected_y["mean"]
            )
            self.assertAlmostEqual(
                missing_cell["{}_hgcer_y_rms".format(prefix)], expected_y["rms"]
            )
            self.assertAlmostEqual(
                missing_cell["{}_hgcer_y_median".format(prefix)], expected_y["q50"]
            )

        identity_fields = (
            "source_label",
            "entry_index",
            "canonical_t_index",
            "delta_index",
            "P_hgcer_npeSum",
            "allcuts",
            "nommcuts",
            "coefficient",
            "diagnostic_weight",
        )
        self.assertEqual(
            [tuple(record[field] for field in identity_fields) for record in finite_records],
            [tuple(record[field] for field in identity_fields) for record in missing_records],
        )

    def test_allcuts_is_secondary_and_detects_biased_response(self):
        identical = _supported_records(low=10, control=30)
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(identical), _phase_a()
        )
        cell = result["cells"][0]
        self.assertEqual(cell["nommcuts_vs_allcuts_status"], "consistent")
        self.assertEqual(cell["f_low"], cell["f_low_allcuts"])

        biased = []
        for index in range(50):
            biased.append(_record(index, 1.0, allcuts=index < 5))
        for index in range(50):
            biased.append(_record(100 + index, 3.0, allcuts=index < 45))
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(biased), _phase_a()
        )
        cell = result["cells"][0]
        self.assertAlmostEqual(cell["f_low"], 0.5)
        self.assertAlmostEqual(cell["f_low_allcuts"], 0.1)
        self.assertEqual(cell["nommcuts_vs_allcuts_status"], "inconsistent")

    def test_final_edges_outside_records_and_assignment_mismatches(self):
        t_edges = (0.0, 1.0, 2.0)
        delta_edges = (0.0, 10.0, 20.0)
        records = _supported_records(
            low=5, control=20, t=2.0, t_index=1, delta=20.0, delta_index=1
        )
        records.extend(
            [
                _record(100, 1.0, analysis_t=-0.1, t_index=None),
                _record(101, 1.0, delta=21.0, delta_index=None),
            ]
        )
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(records, t_edges, delta_edges),
            _phase_a(t_edges, delta_edges),
        )
        self.assertTrue(result["available"])
        final_cell = next(
            cell for cell in result["cells"]
            if cell["t_index"] == 1 and cell["delta_index"] == 1
        )
        self.assertEqual(final_cell["prompt_positive_count"], 25)
        self.assertEqual(result["summary"]["records_outside_t"], 1)
        self.assertEqual(result["summary"]["records_outside_delta"], 1)

        mismatch = deepcopy(records)
        mismatch[0]["canonical_t_index"] = 0
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(mismatch, t_edges, delta_edges),
            _phase_a(t_edges, delta_edges),
        )
        self.assertFalse(result["available"])
        self.assertEqual(result["reason"], "canonical_t_assignment_mismatch")

    def test_no_cross_t_pooling_and_no_interpolation(self):
        t_edges = (0.0, 1.0, 2.0)
        delta_edges = (0.0, 10.0, 20.0, 30.0)
        records = []
        records.extend(_supported_records(delta=5.0, delta_index=0, t=0.5, t_index=0))
        records.extend(_supported_records(delta=25.0, delta_index=2, t=0.5, t_index=0))
        records.extend(_supported_records(delta=5.0, delta_index=0, t=1.5, t_index=1))
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(records, t_edges, delta_edges),
            _phase_a(t_edges, delta_edges),
        )
        middle = next(
            cell for cell in result["cells"]
            if cell["t_index"] == 0 and cell["delta_index"] == 1
        )
        other_t = next(
            cell for cell in result["cells"]
            if cell["t_index"] == 1 and cell["delta_index"] == 0
        )
        self.assertEqual(middle["method_A_status"], "unavailable")
        self.assertIsNone(middle["f_low"])
        self.assertEqual(other_t["prompt_positive_count"], 30)

    def test_provenance_failures_are_unavailable_and_no_proton_factor_is_allowed(self):
        records = _supported_records()
        diagnostic = _diagnostic(records)
        diagnostic["source_provenance"]["pion"]["prompt"]["rf_state"] = "RF_or_unknown"
        result = method_a.build_pion_hgcer_method_a(diagnostic, _phase_a())
        self.assertFalse(result["available"])
        self.assertEqual(result["reason"], "response_not_noRF")

        records[0]["proton_cleaning_factor"] = 0.5
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(records), _phase_a()
        )
        self.assertFalse(result["available"])
        self.assertEqual(result["reason"], "pion_record_has_proton_factor")

        phase = _phase_a()
        phase["available"] = False
        phase["status"] = "unavailable"
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(_supported_records()), phase
        )
        self.assertEqual(result["reason"], "phase_a_contract_unavailable")

    def test_fingerprint_is_phi_independent_and_inputs_remain_immutable(self):
        records = _supported_records()
        diagnostic = _diagnostic(records)
        phase = _phase_a()
        before_diagnostic = deepcopy(diagnostic)
        before_phase = deepcopy(phase)
        first = method_a.build_pion_hgcer_method_a(diagnostic, phase)
        self.assertEqual(diagnostic, before_diagnostic)
        self.assertEqual(phase, before_phase)
        for record in records:
            record["phi"] = -999.0
        second = method_a.build_pion_hgcer_method_a(diagnostic, phase)
        self.assertEqual(first["fingerprint"], second["fingerprint"])
        self.assertEqual(first["event_population_fingerprint"], second["event_population_fingerprint"])
        for original, changed in zip(before_diagnostic["records"]["pion"], records):
            original["phi"] = changed["phi"]
        self.assertEqual(diagnostic, before_diagnostic)
        self.assertEqual(phase, before_phase)
        first["t_edges"][0] = -10.0
        self.assertEqual(phase["canonical_t_edges"][0], 0.0)

        rejected = method_a.build_pion_hgcer_method_a(
            diagnostic, phase, config={"support": {"phi_bins": 8}}
        )
        self.assertFalse(rejected["available"])
        self.assertEqual(rejected["reason"], "phi_dependent_configuration_forbidden")

    def test_json_safety_summary_and_forbidden_output_contract(self):
        result = method_a.build_pion_hgcer_method_a(
            _diagnostic(_supported_records()), _phase_a()
        )
        encoded = json.dumps(result, allow_nan=False)
        summary = method_a.summarize_pion_hgcer_method_a(result)
        json.dumps(summary, allow_nan=False)
        self.assertNotIn("cells", summary)
        self.assertIn('"schema_version": "pion_hgcer_method_a/v1"', encoded)

        forbidden = {"C_A", "C", "C_final", "refined_pion_weight"}

        def walk(value):
            if isinstance(value, dict):
                self.assertTrue(forbidden.isdisjoint(value))
                for child in value.values():
                    walk(child)
            elif isinstance(value, list):
                for child in value:
                    walk(child)

        walk(result)
        source = (REPO_ROOT / "src/cuts/pion_hgcer_refinement_method_a.py").read_text(
            encoding="utf-8"
        )
        for token in ("TCanvas", ".Print(", "write_", "pion_hgcer_transfer"):
            self.assertNotIn(token, source)

    def test_static_handoff_and_runtime_order(self):
        diagnostic_source = (REPO_ROOT / "src/cuts/pion_hgcer_diagnostics.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('"P_hgcer_xAtCer": _finite(', diagnostic_source)
        self.assertIn('"P_hgcer_yAtCer": _finite(', diagnostic_source)
        runtime_source = (REPO_ROOT / "src/cuts/rand_sub.py").read_text(encoding="utf-8")
        event_position = runtime_source.index("build_pion_hgcer_event_contract(")
        method_position = runtime_source.index("build_pion_hgcer_method_a(")
        self.assertLess(event_position, method_position)
        self.assertIn('histDict["_pion_hgcer_method_a"]', runtime_source)
        self.assertNotIn("pion_hgcer_method_a_json", runtime_source)


if __name__ == "__main__":
    unittest.main()
