"""Presentation ownership and visibility checks; no ROOT/farm claim."""
from copy import deepcopy
import tempfile
import unittest
from unittest.mock import patch

from testing import test_e8_4_production_impact_audit as fixtures
from testing import test_e8_3_detached_method_a_reweighting_audit as historical

plots = fixtures.plots


def snapshot(payload):
    def freeze(value):
        if isinstance(value, fixtures._Histogram):
            return (tuple(value.contents), tuple(value.errors), tuple(value.edges), value.directory)
        if isinstance(value, dict):
            return {key: freeze(item) for key, item in value.items()}
        if isinstance(value, (tuple, list)):
            return tuple(freeze(item) for item in value)
        return deepcopy(value)
    return freeze(payload)


class VisualizationTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.E84ShareablePagesTests()
        self.inputs = self.fixture.fixtures()
        self.payload = self.fixture.build(*self.inputs)
        self.assertTrue(self.payload["available"], self.payload.get("reason"))
        fixtures._Histogram.draw_calls[:] = []
        fixtures._FakeText.lines[:] = []

    def capture(self, renderer, payload=None):
        payload = self.payload if payload is None else payload
        drawn = []
        original = plots._e8_4_draw_histogram
        def capture(*args, **kwargs):
            result = original(*args, **kwargs)
            drawn.append(result)
            return result
        with patch.object(plots, "_e8_4_draw_histogram", side_effect=capture):
            self.assertTrue(renderer(fixtures._FakeRoot, "ignored.pdf", payload, payload["per_t"][0]))
        return drawn

    def assert_comparison(self, method_a, baseline):
        self.assertEqual((method_a.color, method_a.line_style, method_a.width), (6, 2, 4))
        self.assertEqual((baseline.color, baseline.line_style, baseline.width), (4, 1, 2))
        self.assertEqual((baseline.marker_style, baseline.marker_color, baseline.marker_size), (24, 4, 0.6))
        self.assertIn("p", baseline.draw_option.split())
        self.assertEqual((method_a.minimum, method_a.maximum), (baseline.minimum, baseline.maximum))

    def test_final_mm_identical_curves_baseline_last_and_redundant_style(self):
        source = fixtures._f6_3_source(y0=(3.0, 4.0), ya=(3.0, 4.0))
        payload = fixtures._build(source, fixtures._e8_2_payload())
        before = snapshot(payload)
        drawn = self.capture(plots._e8_4_render_final_mm_page, payload)
        for method_a, baseline in zip(drawn[::2], drawn[1::2]):
            self.assertIn("mma_", method_a.clone_name)
            self.assertIn("mm0_", baseline.clone_name)
            self.assertEqual(method_a.contents, baseline.contents)
            self.assert_comparison(method_a, baseline)
        self.assertEqual(snapshot(payload), before)

    def test_pion_input_then_method_a_then_baseline(self):
        drawn = self.capture(plots._e8_4_render_pion_consequence_page)
        for common, method_a, baseline in zip(drawn[::3], drawn[1::3], drawn[2::3]):
            self.assertIn("pion_input", common.clone_name)
            self.assertIn("bpia_", method_a.clone_name)
            self.assertIn("bpi0_", baseline.clone_name)
            self.assertEqual((common.color, common.line_style, common.width), (1, 3, 1))
            self.assert_comparison(method_a, baseline)

    def test_authorized_simc_first_method_a_then_baseline(self):
        drawn = self.capture(plots._e8_4_render_baseline_simc_page)
        for simc, method_a, baseline in zip(drawn[::3], drawn[1::3], drawn[2::3]):
            self.assertIn("SIMC", simc.clone_name)
            self.assertEqual((simc.color, simc.line_style, simc.width), (1, 1, 2))
            self.assert_comparison(method_a, baseline)
        drawn = self.capture(plots._e8_4_render_method_a_simc_page)
        self.assertTrue(all("SIMC" in hist.clone_name for hist in drawn[::2]))
        self.assertTrue(all("MM_A" in hist.clone_name for hist in drawn[1::2]))

    def test_all_pages_preserve_exact_payload_and_page_inventory(self):
        before = snapshot(self.payload)
        source_before = fixtures._source_histogram_snapshot(self.inputs[0])
        manifest, failures = [], []
        plots._render_full_background_subtraction_e8_4_pages(fixtures._FakeRoot, "ignored.pdf", self.payload, manifest, failures)
        from testing.run_e8_4_fix5_left_lowe_plot_gate import NEW_PAGE_IDS, OLD_PAGES
        expected = set(OLD_PAGES) | {(page_id, "t" + page_id.rsplit("t", 1)[1]) for page_id in NEW_PAGE_IDS if ".t" in page_id}
        expected.add(("full_background.e8_4.parent_closure", "setting"))
        self.assertEqual({(row["page_id"], row["scope"]) for row in manifest}, expected)
        self.assertEqual(len(manifest), 24)
        self.assertEqual(failures, [])
        self.assertEqual(snapshot(self.payload), before)
        self.assertEqual(fixtures._source_histogram_snapshot(self.inputs[0]), source_before)

    def test_blocked_simc_six_placeholders_no_simc_draws(self):
        source, baseline, support, closure = self.fixture.fixtures()
        audit = support["normalization_audit"]
        audit.update(absolute_comparison_available=False, absolute_units_source_authority=None,
                     reason="SIMC_normfac_luminosity_and_charge_units_not_source_proven")
        audit.pop("fingerprint")
        audit["fingerprint"] = fixtures._audit_producer._f6_3_audit_digest(audit)
        payload = self.fixture.build(source, baseline, support, closure)
        self.assertTrue(payload["available"], payload.get("reason"))
        before = snapshot(payload)
        manifest, failures = [], []
        plots._render_full_background_subtraction_e8_4_pages(fixtures._FakeRoot, "ignored.pdf", payload, manifest, failures)
        self.assertEqual(failures, [])
        self.assertEqual(len(manifest), 24)
        blocked = [row for row in manifest if "simc" in row["page_id"]]
        self.assertEqual(len(blocked), 6)
        for row in blocked:
            self.assertIs(row["available"], False)
            self.assertIs(row["simc_absolute_comparison_available"], False)
            self.assertEqual(row["reason"], payload["simc_absolute_comparison_reason"])
        self.assertTrue(all(row.get("available", True) for row in manifest if "simc" not in row["page_id"]))
        self.assertFalse(any("SIMC" in name for name, _ in fixtures._Histogram.draw_calls))
        self.assertTrue(any(payload["simc_absolute_comparison_reason"] in line for line in fixtures._FakeText.lines))
        self.assertEqual(snapshot(payload), before)

    def test_signed_support_and_aggregate_display_follow_stored_records_only(self):
        # Deliberately distinct stored values test presentation ownership only;
        # this modified fixture is never passed back to scientific validators.
        group = self.payload["per_t"][0]
        for child in group["children"]:
            child["identity_audit"]["signed_support"] = {
                "MM_0": dict(positive_support=1000.25, negative_support=-1000, signed_integral=0.25, absolute_support=2000.25),
                "MM_A": dict(positive_support=900, negative_support=-899.5, signed_integral=0.5, absolute_support=1799.5),
                "delta_MM": dict(positive_support=100, negative_support=-99.75, signed_integral=0.25, absolute_support=199.75),
            }
        self.payload["current_lineage_aggregate_audit"][0].update(
            aggregate_B_pi_0_integral=12345.25, aggregate_B_pi_A_integral=23456.5, aggregate_delta_integral=11111.25)
        before = snapshot(self.payload)
        with patch.object(fixtures._Histogram, "GetBinContent", side_effect=AssertionError("no re-sum")), patch.object(fixtures._Histogram, "GetBinError", side_effect=AssertionError("no re-integrate")):
            self.assertTrue(plots._e8_4_render_yield_summary_page(fixtures._FakeRoot, "ignored.pdf", self.payload, group))
        lines = fixtures._FakeText.lines
        for child in group["children"]:
            self.assertTrue(any(line.startswith("phi {} ".format(child["phi_index"] + 1)) and
                                "MM_0 S=0.25 Abs=2000.25; MM_A S=0.5 Abs=1799.5; delta S=0.25 Abs=199.75" in line for line in lines))
        self.assertIn("Stored aggregate B_pi_0=12345.25; B_pi_A=23456.5; delta=11111.25", lines)
        self.assertIn("current F.6.3 candidate lineage; Lambda/allcut MM-template aggregate", lines)
        self.assertTrue(any("Not the broader F.4" in line for line in lines))
        self.assertEqual(snapshot(self.payload), before)

    def test_missing_stored_support_fails_closed_without_manifest_success(self):
        self.payload["per_t"][0]["children"][0]["identity_audit"].pop("signed_support")
        manifest, failures = [], []
        plots._render_full_background_subtraction_e8_4_pages(fixtures._FakeRoot, "ignored.pdf", self.payload, manifest, failures)
        self.assertNotIn("full_background.e8_4.yield_summary.t1", [row["page_id"] for row in manifest])
        self.assertIn("E.8.4 full_background.e8_4.yield_summary.t1 page unavailable", failures)

    def test_historical_visible_labels_in_all_page_families(self):
        historical.DetachedMethodAReweightingAuditTests.setUpClass()
        historical_fixture = historical.DetachedMethodAReweightingAuditTests()
        with tempfile.TemporaryDirectory() as directory:
            payload = historical_fixture._payload(directory)
        self.assertTrue(payload["available"], payload.get("reason"))
        captured = []
        def text_page(*args, **kwargs):
            captured.extend(args[4])
            return True
        with patch.object(plots, "_e8_text_page", side_effect=text_page), patch.object(plots, "_e8_3_signed_histogram", return_value=historical._FakeHistogram()):
            historical._FakePaveText.instances[:] = []
            manifest, failures = [], []
            plots._render_full_background_subtraction_e8_3_pages(historical._FakeRoot, "ignored.pdf", payload, manifest, failures)
        self.assertEqual(failures, [])
        labels = captured + [line for text in historical._FakePaveText.instances for line in text.lines]
        self.assertGreaterEqual(sum("Historical accepted F.6.1" in line or "historical accepted F.6.1" in line for line in labels), 8)
        self.assertGreaterEqual(sum("not the current F.6.3/E.8.4 candidate lineage" in line for line in labels), 8)
        self.assertTrue(plots._e8_4_render_setting_summary_page(fixtures._FakeRoot, "ignored.pdf", self.payload))
        self.assertTrue(any("stored Fix.5.4 current-lineage audit" in line for line in fixtures._FakeText.lines))


if __name__ == "__main__":
    unittest.main()
