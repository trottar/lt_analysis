"""Current-lineage numerical/provenance identities; deterministic and ROOT-free."""

from copy import deepcopy
import unittest
from unittest.mock import patch

from testing import test_e8_4_production_impact_audit as fixtures


calc = fixtures._audit_producer
plots = fixtures.plots


def _seal(record):
    record.pop("fingerprint", None)
    record["fingerprint"] = calc._f6_3_audit_digest(record)


def _audit(source, final=None):
    return calc._build_f6_3_identity_audit(
        source, final or {(row["t_index"], row["phi_index"]): row["MM_0"] for row in source["children"]},
        [0.0, 1.0], [-180.0, 0.0, 180.0],
    )


class FlowHistogram(fixtures._Histogram):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.underflow = 91.0
        self.overflow = -37.0

    def GetBinContent(self, index):
        if index == 0:
            return self.underflow
        if index == self.GetNbinsX() + 1:
            return self.overflow
        return super().GetBinContent(index)

    def Clone(self, name):
        clone = super().Clone(name)
        clone.underflow, clone.overflow = self.underflow, self.overflow
        return clone


class CurrentLineageIdentityTests(unittest.TestCase):
    def test_exact_algebra_yields_signed_cancellation_and_current_aggregate(self):
        source = fixtures._f6_3_source(y0=(1.0, 0.0), ya=(0.5, 0.0))
        child = source["children"][0]
        child["MM_0"].contents[:] = [10.0, -9.0]
        child["MM_A"].contents[:] = [8.0, -7.5]
        child["B_pi_0"].contents[:] = [2.0, -1.0]
        child["B_pi_A"].contents[:] = [4.0, -2.5]
        before = fixtures._source_histogram_snapshot(source)
        record = _audit(source)
        self.assertTrue(record["identity_passed"])
        row = record["children"][0]
        self.assertEqual(row["maximum_absolute_bin_residual"], 0.0)
        self.assertEqual(row["failing_bin_count"], 0)
        self.assertEqual((row["MM_0_integral"], row["MM_A_integral"], row["delta_MM_integral"]), (1.0, 0.5, -0.5))
        self.assertEqual(row["signed_support"]["MM_0"], {
            "positive_support": 10.0, "negative_support": -9.0,
            "signed_integral": 1.0, "absolute_support": 19.0})
        self.assertEqual(row["signed_support"]["delta_MM"], {
            "positive_support": 1.5, "negative_support": -2.0,
            "signed_integral": -0.5, "absolute_support": 3.5})
        aggregate = record["aggregates"][0]
        self.assertEqual(aggregate["child_inventory"], [[0, 0], [0, 1]])
        self.assertEqual(aggregate["B_pi_0_current"], [5.0, 2.5])
        self.assertEqual(aggregate["B_pi_A_current"], [7.0, 1.0])
        self.assertFalse(aggregate["candidate_f4_parent_sum_comparable"])
        self.assertIn("nommcut", aggregate["candidate_f4_parent_sum_reason"])
        self.assertEqual(fixtures._source_histogram_snapshot(source), before)

    def test_one_bin_mm_or_pion_perturbation_fails(self):
        for name in ("MM_A", "B_pi_A"):
            with self.subTest(histogram=name):
                source = fixtures._f6_3_source()
                source["children"][0][name].contents[0] += 0.01
                record = _audit(source)
                self.assertFalse(record["identity_passed"])
                self.assertEqual(record["children"][0]["failing_bin_count"], 1)

    def test_empty_current_child_remains_in_inventory_without_normalization(self):
        source = fixtures._f6_3_source()
        child = source["children"][1]
        for name in ("B_pi_0", "B_pi_A", "MM_0", "MM_A"):
            child[name].contents[:] = [0.0, 0.0]
        child.update(Y0=0.0, YA=0.0)
        audit = _audit(source)
        self.assertTrue(audit["identity_passed"])
        self.assertFalse(audit["children"][1]["populated"])
        self.assertEqual(audit["aggregates"][0]["child_inventory"], [[0, 0], [0, 1]])
        self.assertEqual(audit["aggregates"][0]["populated_child_inventory"], [[0, 0]])

    def test_y0_and_ya_mismatch_fail_without_replacing_scalar(self):
        for name in ("Y0", "YA"):
            source = fixtures._f6_3_source()
            source["children"][0][name] += 1.0
            before = source["children"][0][name]
            record = _audit(source)
            self.assertFalse(record["identity_passed"])
            self.assertFalse(record["children"][0]["yield_identity_passed"])
            self.assertEqual(source["children"][0][name], before)

    def test_geometry_and_final_baseline_content_or_error_mismatch_fail(self):
        for defect in ("axis", "content", "error"):
            source = fixtures._f6_3_source()
            final = {(row["t_index"], row["phi_index"]): row["MM_0"].Clone("final") for row in source["children"]}
            if defect == "axis":
                source["children"][0]["B_pi_A"].edges[1] += 0.001
            elif defect == "content":
                final[(0, 0)].contents[0] += 1.0
            else:
                final[(0, 0)].errors[0] += 1.0
            with self.subTest(defect=defect), self.assertRaisesRegex(RuntimeError, "binning|fingerprint"):
                _audit(source, final)

    def test_duplicate_or_missing_child_fails(self):
        for duplicate in (False, True):
            source = fixtures._f6_3_source()
            source["children"] = [source["children"][0]] * (2 if duplicate else 1)
            with self.assertRaisesRegex(RuntimeError, "inventory"):
                _audit(source)

    def test_flow_bins_excluded_from_arithmetic_but_fingerprinted(self):
        source = fixtures._f6_3_source()
        for row in source["children"]:
            for name in ("B_pi_0", "B_pi_A", "MM_0", "MM_A"):
                old = row[name]
                row[name] = FlowHistogram(old.contents, edges=old.edges, errors=old.errors)
        record = _audit(source)
        self.assertTrue(record["identity_passed"])
        self.assertEqual(record["children"][0]["MM_0_integral"], source["children"][0]["Y0"])
        final = {(row["t_index"], row["phi_index"]): row["MM_0"].Clone("final") for row in source["children"]}
        final[(0, 0)].overflow += 1.0
        with self.assertRaisesRegex(RuntimeError, "fingerprint"):
            _audit(source, final)

    def test_strict_scaled_tolerance_and_nonfinite_rejection(self):
        source = fixtures._f6_3_source()
        source["children"][0]["MM_A"].contents[0] += 1.0e-14
        self.assertTrue(_audit(source)["identity_passed"])
        source["children"][0]["MM_A"].contents[0] += 1.0e-7
        self.assertFalse(_audit(source)["identity_passed"])
        source["children"][0]["B_pi_A"].contents[0] = float("nan")
        with self.assertRaises(RuntimeError):
            _audit(source)

    def test_runtime_attachment_uses_exact_measured_final_and_fails_closed(self):
        source = fixtures._f6_3_source()
        measurements = {(row["t_index"], row["phi_index"]): {"final_histogram": row["MM_0"].Clone("measured")} for row in source["children"]}
        calc._attach_f6_3_identity_audit(source, measurements, [0.0, 1.0], [-180.0, 0.0, 180.0])
        self.assertTrue(source["available"])
        measurements[(0, 0)]["final_histogram"].contents[0] += 0.01
        calc._attach_f6_3_identity_audit(source, measurements, [0.0, 1.0], [-180.0, 0.0, 180.0])
        self.assertFalse(source["available"])
        self.assertIn("final_baseline_fingerprint", source["reason"])


class ConsumerProvenanceTests(unittest.TestCase):
    def test_unproven_simc_retains_all_three_t_parents_and_six_blocked_page_ids(self):
        source, baseline, support, closure = fixtures.E84ShareablePagesTests().fixtures()
        audit = support["normalization_audit"]
        audit.update(absolute_comparison_available=False,
                     histogram_units="sum_iter_weight_times_normfac_divided_by_Ncontribute")
        audit.pop("absolute_units_source_authority")
        _seal(audit)
        payload = plots.build_full_background_subtraction_e8_4_payload(
            source, baseline, simc_support=support, parent_closure=closure)
        self.assertTrue(payload["available"])
        manifest, failures = [], []
        with patch.object(plots, "_e8_4_draw_histogram", wraps=plots._e8_4_draw_histogram) as draw:
            plots._render_full_background_subtraction_e8_4_pages(fixtures._FakeRoot, "fixture.pdf", payload, manifest, failures)
        self.assertEqual(failures, [])
        self.assertFalse(any("share_" in call.args[2] for call in draw.call_args_list))
        blocked = [page for page in manifest if "simc" in page["page_id"]]
        self.assertEqual({page["page_id"] for page in blocked}, {
            "full_background.e8_4.{}.t{}".format(family, t)
            for family in ("method_a_vs_simc", "baseline_method_a_simc") for t in (1, 2, 3)})
        self.assertTrue(all(page["available"] is False and page["reason"] == audit["absolute_comparison_reason"] for page in blocked))
        self.assertEqual(len(payload["current_lineage_aggregate_audit"]), 3)

    def build(self, source, support=None):
        baseline = fixtures._e8_2_payload()
        return plots.build_full_background_subtraction_e8_4_payload(
            source, baseline, simc_support=support or fixtures._support(source, baseline),
            parent_closure=fixtures._closure())

    def test_comparable_synthetic_simc_accepted_without_scaling_or_mutation(self):
        source = fixtures._f6_3_source()
        support = fixtures._support(source, fixtures._e8_2_payload())
        before = fixtures._source_histogram_snapshot(source)
        simc_before = deepcopy(support["mm"][0][0].contents)
        payload = self.build(source, support)
        self.assertTrue(payload["available"])
        self.assertTrue(payload["simc_absolute_comparison_available"])
        self.assertTrue(payload["current_lineage_identity_audit"]["identity_passed"])
        self.assertEqual(payload["simc_normalization_audit"]["children"][0]["SIMC_integral_in_the_same_lambda_window"], 12.5)
        self.assertFalse(payload["method_b_numerical_dependency"])
        self.assertFalse(payload["canonical_child_renormalization_performed"])
        self.assertEqual(fixtures._source_histogram_snapshot(source), before)
        self.assertEqual(support["mm"][0][0].contents, simc_before)
        payload["current_lineage_aggregate_audit"][0]["B_pi_0_current"][0] = 999.0
        self.assertNotEqual(source["identity_audit"]["aggregates"][0]["B_pi_0_current"][0], 999.0)

    def test_missing_stale_historical_or_inconsistent_current_audit_rejected(self):
        mutations = (
            lambda source: source.pop("identity_audit"),
            lambda source: source["identity_audit"].update(authority={"historical": True}),
            lambda source: source["identity_audit"].update(lineage="historical_f6_1"),
            lambda source: source["identity_audit"]["aggregates"][0].update(lineage="historical_f6_1"),
            lambda source: source["identity_audit"]["aggregates"][0]["B_pi_0_current"].__setitem__(0, 999.0),
            lambda source: source["identity_audit"]["children"][0]["histogram_fingerprints"].update(final_baseline="wrong"),
            lambda source: source["children"][0]["MM_A"].contents.__setitem__(0, 99.0),
        )
        for mutation in mutations:
            source = fixtures._f6_3_source()
            mutation(source)
            if "identity_audit" in source:
                _seal(source["identity_audit"])
            with self.subTest(mutation=mutation):
                self.assertFalse(self.build(source)["available"])

    def test_simc_missing_ambiguous_units_or_display_normalization_rejected(self):
        for defect in ("missing", "units", "fallback", "object", "factor", "missing_unit_authority"):
            source = fixtures._f6_3_source()
            support = fixtures._support(source, fixtures._e8_2_payload())
            audit = support["normalization_audit"]
            if defect == "missing":
                support.pop("normalization_audit")
            elif defect == "units":
                audit["histogram_units"] = "ambiguous"
            elif defect == "fallback":
                audit["display_normalization_applied"] = True
            elif defect == "object":
                audit["children"][0]["histogram_fingerprint"] = "wrong"
            elif defect == "factor":
                audit["normalization_factor_applied"] *= 2.0
            else:
                audit.pop("absolute_units_source_authority")
            _seal(audit)
            with self.subTest(defect=defect):
                self.assertFalse(self.build(source, support)["available"])

    def test_source_unproven_simc_blocks_only_absolute_pages_and_keeps_current_audit(self):
        source = fixtures._f6_3_source()
        support = fixtures._support(source, fixtures._e8_2_payload())
        audit = support["normalization_audit"]
        audit.update(absolute_comparison_available=False,
                     histogram_units="sum_iter_weight_times_normfac_divided_by_Ncontribute")
        audit.pop("absolute_units_source_authority")
        _seal(audit)
        before = fixtures._source_histogram_snapshot(source)
        simc_before = deepcopy(support["mm"][0][0].contents)
        payload = self.build(source, support)
        self.assertTrue(payload["available"])
        self.assertIsNone(payload["reason"])
        self.assertFalse(payload["simc_absolute_comparison_available"])
        self.assertEqual(payload["simc_absolute_comparison_reason"], "SIMC_normfac_luminosity_and_charge_units_not_source_proven")
        self.assertEqual(payload["simc_normalization_audit"], audit)
        self.assertTrue(payload["current_lineage_identity_audit"]["identity_passed"])
        self.assertEqual(len(payload["current_lineage_aggregate_audit"]), 1)
        fixtures._FakeText.lines[:] = []
        fixtures._Histogram.draw_calls[:] = []
        manifest, failures = [], []
        with patch.object(plots, "_e8_4_draw_histogram", wraps=plots._e8_4_draw_histogram) as draw:
            plots._render_full_background_subtraction_e8_4_pages(fixtures._FakeRoot, "fixture.pdf", payload, manifest, failures)
        self.assertEqual(failures, [])
        self.assertTrue(draw.called)  # Data-only panels were retained.
        ids = {row["page_id"] for row in manifest}
        for suffix in ("authority", "parent_closure", "setting_summary", "pion_consequence", "final_mm", "signed_difference", "yield_impact", "yield_summary.t1"):
            self.assertIn("full_background.e8_4." + suffix, ids)
        blocked = [row for row in manifest if "simc" in row["page_id"]]
        self.assertEqual(len(blocked), 2)
        for row in blocked:
            self.assertFalse(row["available"])
            self.assertFalse(row["simc_absolute_comparison_available"])
            self.assertEqual(row["reason"], payload["simc_absolute_comparison_reason"])
            self.assertEqual(row["represented_phi_inventory"][0]["phi_index"], 0)
        self.assertTrue(any(payload["simc_absolute_comparison_reason"] in text for text in fixtures._FakeText.lines))
        self.assertFalse(any("share_" in name for name, option in fixtures._Histogram.draw_calls))
        self.assertEqual(fixtures._source_histogram_snapshot(source), before)
        self.assertEqual(support["mm"][0][0].contents, simc_before)


if __name__ == "__main__":
    unittest.main()
