"""Presentation ownership and visibility checks; no ROOT/farm claim."""
from copy import deepcopy
import ast
import inspect
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
        ordered = [("full_background.e8_4.authority", "setting")]
        for index in range(1, 4):
            scope = "t{}".format(index)
            ordered.extend(("full_background.e8_4." + suffix, scope) for suffix in
                           ("pion_consequence", "final_mm", "signed_difference", "yield_impact"))
            ordered.extend(("full_background.e8_4.{}.{}".format(suffix, scope), scope)
                           for suffix in ("method_a_vs_simc", "baseline_method_a_simc", "yield_summary"))
        ordered.extend(("full_background.e8_4." + suffix, "setting")
                       for suffix in ("parent_closure", "setting_summary"))
        self.assertEqual([(row["page_id"], row["scope"]) for row in manifest], ordered)
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
        for index, child in enumerate(group["children"]):
            child["identity_audit"]["signed_support"] = {
                "MM_0": dict(positive_support=1000.25, negative_support=-1000, signed_integral=0.25 + index, absolute_support=2000.25 + index),
                "MM_A": dict(positive_support=900, negative_support=-899.5, signed_integral=0.5 + index, absolute_support=1799.5 + index),
                "delta_MM": dict(positive_support=100, negative_support=-99.75, signed_integral=0.75 + index, absolute_support=199.75 + index),
            }
        self.payload["current_lineage_aggregate_audit"][0].update(
            aggregate_B_pi_0_integral=12345.25, aggregate_B_pi_A_integral=23456.5, aggregate_delta_integral=11111.25)
        before = snapshot(self.payload)
        with patch.object(fixtures._Histogram, "GetBinContent", side_effect=AssertionError("no re-sum")), patch.object(fixtures._Histogram, "GetBinError", side_effect=AssertionError("no re-integrate")):
            self.assertTrue(plots._e8_4_render_yield_summary_page(fixtures._FakeRoot, "ignored.pdf", self.payload, group))
        lines = fixtures._FakeText.lines
        support_lines = plots._e8_4_stored_support_lines(self.payload, group)
        self.assertTrue(all(len(line) <= 96 for line in support_lines), support_lines)
        self.assertEqual(len(support_lines), 23)
        for child in group["children"]:
            row = next(i for i, line in enumerate(support_lines)
                       if line.startswith("phi {} ".format(child["phi_index"] + 1)))
            self.assertIn("[{:.0f}, {:.0f})".format(child["phi_low"], child["phi_high"]), support_lines[row])
            for key, label, offset in (("MM_0", "MM_0", 0), ("MM_A", "MM_A", 1), ("delta_MM", "delta", 1)):
                stored = child["identity_audit"]["signed_support"][key]
                self.assertIn("{} S={:.8g} Abs={:.8g}".format(
                    label, stored["signed_integral"], stored["absolute_support"]), support_lines[row + offset])
            self.assertIn(support_lines[row], lines)
            self.assertIn(support_lines[row + 1], lines)
        self.assertIn("Stored aggregate B_pi_0=12345.25; B_pi_A=23456.5; delta=11111.25", lines)
        self.assertIn("current F.6.3 candidate lineage; Lambda/allcut MM-template aggregate", lines)
        self.assertTrue(any("Not the broader F.4" in line for line in lines))
        self.assertEqual(snapshot(self.payload), before)

    def test_stored_support_exponent_tokens_remain_bounded_at_original_precision(self):
        group = self.payload["per_t"][0]
        for child in group["children"]:
            for stored in child["identity_audit"]["signed_support"].values():
                stored.update(signed_integral=-1.2345678e-307, absolute_support=1.2345678e307)
        self.payload["current_lineage_aggregate_audit"][0].update(
            aggregate_B_pi_0_integral=-1.2345678e-307,
            aggregate_B_pi_A_integral=1.2345678e307, aggregate_delta_integral=-1.2345678e307)
        before = snapshot(self.payload)
        lines = plots._e8_4_stored_support_lines(self.payload, group)
        self.assertTrue(all(len(line) <= 96 for line in lines), lines)
        self.assertIn("Stored aggregate B_pi_0=-1.2345678e-307; B_pi_A=1.2345678e+307; delta=-1.2345678e+307", lines)
        self.assertEqual(snapshot(self.payload), before)

    def test_authority_flags_and_statements_are_bounded_and_immutable(self):
        flags = ("baseline_production_mutated", "production_promotion_performed",
                 "method_b_numerical_dependency", "empirical_residual_used",
                 "event_correction_persisted", "canonical_child_renormalization_performed",
                 "baseline_public_output_unchanged")
        # Distinct stored Boolean values prove this is formatting, not a default.
        for index, name in enumerate(flags):
            self.payload[name] = bool(index % 2)
        before = snapshot(self.payload)
        with patch.object(plots, "_e8_text_page", return_value=True) as page:
            self.assertTrue(plots._e8_4_render_authority_page(fixtures._FakeRoot, "ignored.pdf", self.payload))
        lines = page.call_args.args[4]
        self.assertTrue(all(len(line) <= 96 for line in lines), lines)
        self.assertLessEqual(len(page.call_args.args[3]), 96)
        self.assertIn("Non-production flags:", lines)
        for name in flags:
            self.assertEqual(sum(name in line for line in lines), 1)
            self.assertIn("  {}={}".format(name, self.payload[name]), lines)
        joined = " ".join(lines)
        self.assertIn("historical E.8.3/F.6.1 is a separate lineage.", joined)
        self.assertIn("Method B is numerically absent.", joined)
        self.assertIn("not automatically a systematic uncertainty.", joined)
        self.assertEqual(snapshot(self.payload), before)

    def test_yield_support_region_fits_nine_children_without_panel_overlap(self):
        before = snapshot(self.payload)
        with patch.object(plots, "_e8_add_text", wraps=plots._e8_add_text) as text, patch.object(
                fixtures._FakeCanvas, "SetPad") as pads:
            self.assertTrue(plots._e8_4_render_yield_summary_page(
                fixtures._FakeRoot, "ignored.pdf", self.payload, self.payload["per_t"][0]))
        support = next(call for call in text.call_args_list if call.args[1] == (0.04, 0.03, 0.96, 0.54))
        self.assertEqual(support.kwargs["size"], 0.021)
        self.assertEqual(len(support.args[2]), 23)
        self.assertGreater(support.args[1][3] - support.args[1][1], len(support.args[2]) * support.kwargs["size"])
        self.assertEqual(len(pads.call_args_list), 3)
        for call in pads.call_args_list:
            self.assertGreater(call.args[1], support.args[1][3])
            self.assertEqual(call.args[3], 0.90)
        self.assertEqual(snapshot(self.payload), before)

    def test_target_titles_are_ascii_safe_and_retain_identity(self):
        historical.DetachedMethodAReweightingAuditTests.setUpClass()
        instance = historical.DetachedMethodAReweightingAuditTests()
        with tempfile.TemporaryDirectory() as directory:
            payload = instance._payload(directory)
        before = deepcopy(payload)
        historical._FakePaveText.instances[:] = []
        for parent in payload["parents"]:
            self.assertTrue(plots._e8_3_render_tphi_page(historical._FakeRoot, "ignored.pdf", payload, parent))
            title = historical._FakePaveText.instances[-1].lines[0]
            self.assertEqual(title, "E.8.3 persisted canonical (t,phi) pion redistribution - {} t{}".format(
                payload["setting_id"], parent["canonical_t_index"] + 1))
            self.assertNotIn("\u2014", title)
            self.assertNotIn("\u00e2\u20ac\u201d", title)
        self.assertEqual(payload, before)
        before = snapshot(self.payload)
        fixtures._FakeText.lines[:] = []
        self.assertTrue(plots._e8_4_render_setting_summary_page(fixtures._FakeRoot, "ignored.pdf", self.payload))
        title = fixtures._FakeText.lines[0]
        self.assertEqual(title, "E.8.4 stored canonical t-by-phi impact summary - {}".format(self.payload["setting_id"]))
        self.assertNotIn("\u2014", title)
        self.assertNotIn("\u00e2\u20ac\u201d", title)
        self.assertEqual(snapshot(self.payload), before)

    def test_affected_renderers_have_no_numerical_producer_or_integration_calls(self):
        forbidden = {"GetBinContent", "GetBinError", "Integral", "Scale", "Normalize",
                     "Rebin", "Fill", "Add", "calculate_yield_data", "calculate_yield_simc",
                     "fill_simc_shape_pion_subtraction_templates", "bin_data", "bin_simc", "bg_fit"}
        for function in (plots._e8_3_render_tphi_page, plots._e8_4_render_setting_summary_page,
                         plots._e8_4_render_authority_page, plots._e8_4_stored_support_lines,
                         plots._e8_4_render_yield_summary_page):
            tree = ast.parse(inspect.getsource(function))
            calls = {node.func.attr if isinstance(node.func, ast.Attribute) else node.func.id
                     for node in ast.walk(tree) if isinstance(node, ast.Call)
                     and isinstance(node.func, (ast.Name, ast.Attribute))}
            self.assertFalse(calls & forbidden, (function.__name__, calls & forbidden))
            if function is plots._e8_4_stored_support_lines:
                self.assertNotIn("sum", calls)

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
