"""Focused deterministic consumer-contract coverage for E.8.4."""

from __future__ import annotations

from copy import deepcopy
import inspect
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO_ROOT / "src" / "utility"), str(REPO_ROOT / "src" / "cuts")]

import full_background_subtraction_plots as plots
from root_histogram_ownership import fingerprint_histogram_content_error
from testing.test_e8_2_baseline_stage_audit import _load_calculate_yield_module


_audit_producer = _load_calculate_yield_module()
_audit_producer.fingerprint_histogram_content_error = fingerprint_histogram_content_error


def _refresh_identity(source):
    children = source["children"]
    source["identity_audit"] = _audit_producer._build_f6_3_identity_audit(
        source, {(row["t_index"], row["phi_index"]): row["MM_0"].Clone("final_yield_histogram") for row in children},
        sorted(set([row["t_low"] for row in children] + [row["t_high"] for row in children])),
        sorted(set([row["phi_low"] for row in children] + [row["phi_high"] for row in children])),
    )
    return source


class _Axis:
    def __init__(self, edges):
        self.edges = list(edges)

    def GetBinLowEdge(self, index):
        return self.edges[int(index) - 1]

    def GetXmin(self):
        return self.edges[0]

    def GetXmax(self):
        return self.edges[-1]


class _Histogram:
    draw_calls = []

    def __init__(self, contents, edges=(1.0, 1.1, 1.2), errors=None):
        self.contents = [float(value) for value in contents]
        self.edges = [float(value) for value in edges]
        self.errors = list(errors or [0.0 for _unused in contents])
        self.directory = object()
        self.title = ""
        self.clone_name = ""

    def Clone(self, name):
        clone = type(self)(self.contents, self.edges, self.errors)
        clone.clone_name = str(name)
        return clone

    def SetDirectory(self, directory):
        self.directory = directory

    def GetNbinsX(self):
        return len(self.contents)

    def GetBinContent(self, index):
        if int(index) in (0, len(self.contents) + 1):
            return 0.0
        return self.contents[int(index) - 1]

    def GetBinError(self, index):
        if int(index) in (0, len(self.contents) + 1):
            return 0.0
        return self.errors[int(index) - 1]

    def GetXaxis(self):
        return _Axis(self.edges)

    def Add(self, other, factor=1.0):
        self.contents = [
            left + float(factor) * right
            for left, right in zip(self.contents, other.contents)
        ]

    def GetMinimum(self):
        return min(self.contents)

    def GetMaximum(self):
        return max(self.contents)

    def SetMinimum(self, value):
        self.minimum = value

    def SetMaximum(self, value):
        self.maximum = value

    def SetTitle(self, value):
        self.title = str(value)

    def GetTitle(self):
        return self.title

    def SetLineColor(self, value):
        self.color = value

    def SetLineWidth(self, value):
        self.width = value

    def SetStats(self, value):
        self.stats = value

    def Draw(self, option):
        self.draw_option = option
        type(self).draw_calls.append((self.clone_name, option))


class _FakeCanvas:
    printed = []

    def __init__(self, *args):
        self.args = args

    def Divide(self, *_args):
        return None

    def SetPad(self, *args):
        return None

    def SetTopMargin(self, value):
        return None

    def SetBottomMargin(self, value):
        return None

    def DrawFrame(self, *args):
        return _FakeGraph()

    def cd(self, *_args):
        return self

    def Print(self, name):
        self.printed.append(name)

    def Close(self):
        return None


class _FakeText:
    lines = []

    def __init__(self, *_args):
        self.local_lines = []

    def SetFillStyle(self, _value):
        return None

    def SetBorderSize(self, _value):
        return None

    def SetTextAlign(self, _value):
        return None

    def SetTextSize(self, _value):
        return None

    def AddText(self, text):
        self.local_lines.append(str(text))
        self.lines.append(str(text))

    def Draw(self):
        return None


class _FakeLine:
    def __init__(self, *_args):
        pass

    def SetLineStyle(self, _value):
        return None

    def SetLineColor(self, _value):
        return None

    def Draw(self):
        return None


class _FakeGraph:
    options = []
    points = []

    def __init__(self, *args):
        self.args = args
        if len(args) >= 3:
            self.points.append((list(args[1]), list(args[2])))

    def SetMarkerStyle(self, value):
        return None

    def SetMarkerColor(self, value):
        return None

    def SetTitle(self, value):
        self.title = value

    def SetLineColor(self, _value):
        return None

    def SetLineWidth(self, _value):
        return None

    def SetStats(self, _value):
        return None

    def Draw(self, option):
        self.options.append(option)


class _FakeLegend(_FakeText):
    entries = []
    def AddEntry(self, hist, label, option):
        self.entries.append(label)


class _FakeRoot:
    TCanvas = _FakeCanvas
    TPaveText = _FakeText
    TLine = _FakeLine
    TGraph = _FakeGraph
    TGraphErrors = _FakeGraph
    TLegend = _FakeLegend
    kBlack = 1
    kBlue = 4
    kMagenta = 6


def _histogram(value, *, errors=(0.1, 0.2), edges=(1.10, 1.13, 1.16)):
    return _Histogram((value, value + 0.5), edges=edges, errors=errors)


def _e8_2_payload(
    *, y0=(3.0, 4.0), setting="Left", epsilon="low",
    mm_edges=(0.70, 1.00, 1.30), t_edges=(0.0, 1.0), phi_edges=(-180.0, 0.0, 180.0),
):
    per_t = []
    for t_index in range(len(t_edges) - 1):
        children = []
        for phi_index in range(len(phi_edges) - 1):
            value = y0[(t_index * (len(phi_edges) - 1) + phi_index) % len(y0)]
            children.append({
                "t_index": t_index, "t_low": t_edges[t_index], "t_high": t_edges[t_index + 1],
                "phi_index": phi_index, "phi_low": phi_edges[phi_index],
                "phi_high": phi_edges[phi_index + 1], "valid": True,
                "final_yield": float(value), "statistical_error": 0.2,
                "total_error": 0.3,
            })
        per_t.append({
            "t_index": t_index, "t_low": t_edges[t_index], "t_high": t_edges[t_index + 1],
            "children": tuple(children),
        })
    return {
        "schema_version": plots.E8_2_PRESENTATION_SCHEMA_VERSION,
        "available": True, "reason": None, "setting": setting, "epsilon": epsilon,
        "t_edges": list(t_edges), "phi_edges": list(phi_edges),
        "mm_edges": list(mm_edges), "lambda_window": [1.08, 1.18],
        "per_t": tuple(per_t),
    }


def _f6_3_source(
    *, y0=(3.0, 4.0), ya=(3.6, 4.0), reverse=False,
    selected_setting_id="Left-lowe", mm_edges=(1.10, 1.13, 1.16),
    t_edges=(0.0, 1.0), phi_edges=(-180.0, 0.0, 180.0),
):
    children = []
    for t_index in range(len(t_edges) - 1):
        for phi_index in range(len(phi_edges) - 1):
            value_index = t_index * (len(phi_edges) - 1) + phi_index
            baseline = y0[value_index % len(y0)]
            adjusted = ya[value_index % len(ya)]
            child_offset = value_index
            children.append({
                "t_index": t_index, "t_low": t_edges[t_index], "t_high": t_edges[t_index + 1],
                "phi_index": phi_index, "phi_low": phi_edges[phi_index],
                "phi_high": phi_edges[phi_index + 1], "valid": True, "status": "available",
                "pion_input": _Histogram((baseline + 2.0 + child_offset, 2.5 + child_offset), edges=mm_edges, errors=(0.1, 0.2)),
                "B_pi_0": _histogram(2.0 + child_offset, edges=mm_edges),
                "B_pi_A": _Histogram((2.0 + child_offset - (adjusted - baseline), 2.5 + child_offset), edges=mm_edges, errors=(0.1, 0.2)),
                "MM_0": _Histogram((baseline, 0.0), edges=mm_edges, errors=(0.1, 0.2)),
                "MM_A": _Histogram((adjusted, 0.0), edges=mm_edges, errors=(0.1, 0.2)),
                "Y0": float(baseline), "Y0_statistical_error": 0.2,
                "Y0_total_error": 0.3, "YA": float(adjusted),
                "YA_statistical_error": 0.25, "YA_total_error": 0.35,
            })
    if reverse:
        children.reverse()
    return _refresh_identity({
        "schema_version": "f6_3_parallel_method_a_source/v1", "available": True,
        "reason": None, "selected_setting_id": selected_setting_id,
        "branch_role": "parallel_nonproduction_method_a_full_analysis",
        "baseline_production_mutated": False, "production_promotion_performed": False,
        "method_b_numerical_dependency": False, "empirical_residual_used": False,
        "event_correction_persisted": False,
        "canonical_child_renormalization_performed": False,
        "baseline_public_output_unchanged": True,
        "authority": {"live_cache_parity_passed": True, "accepted": "fixture"},
        "lambda_integration_window": [1.08, 1.18], "children": children,
    })


def _render_payloads():
    return {
        key: {"available": False, "reason": "fixture"}
        for key in ("d6", "d7", "d8", "d9", "d10", "d11", "e2", "e3", "e4", "e6", "e7", "f1", "e72", "e8")
    }


def _source_histogram_snapshot(source):
    return json.dumps([
        {
            name: {
                "contents": list(child[name].contents),
                "errors": list(child[name].errors),
                "edges": list(child[name].edges),
            }
            for name in ("pion_input", "B_pi_0", "B_pi_A", "MM_0", "MM_A")
        }
        for child in source["children"]
    ], sort_keys=True)


def _closure(count=1):
    return [{"canonical_t_index": index, "closure_passed": True,
             "baseline_parent_sum": 12.0 + index, "adjusted_parent_sum": 12.0 + index,
             "closure_residual": 0.0, "closure_relative_scale": 1.0e-10}
            for index in range(count)]


def _support(source, baseline):
    edges = (1.10, 1.13, 1.16)
    if isinstance(source, dict) and isinstance(source.get("children"), (list, tuple)) and source["children"]:
        edges = source["children"][0]["pion_input"].edges
    support = {"mm": [[_histogram(6.0, edges=edges) for _ in range(len(baseline["phi_edges"]) - 1)]
                      for _ in range(len(baseline["t_edges"]) - 1)]}
    record = _audit_producer._build_f6_3_simc_provenance(
        {"phi_setting": baseline["setting"], "InFile_SIMC": type("RootFile", (), {"GetName": lambda self: "synthetic.root"})(),
         "normfac_simc": 0.5, "simc_normfactor": 5.0, "simc_nevents": 10},
        support, baseline["t_edges"], baseline["phi_edges"],
        {"EPSSET": baseline["epsilon"], "mm_min": 1.08, "mm_max": 1.18},
    )
    # This is explicit synthetic unit authority, never a farm/SIMC claim.
    record.update(absolute_comparison_available=True, histogram_units="yield_per_effective_data_charge",
                  absolute_units_source_authority="synthetic_deterministic_fixture")
    record.pop("fingerprint")
    record["fingerprint"] = _audit_producer._f6_3_audit_digest(record)
    support["normalization_audit"] = record
    return support


def _build(source, baseline, **kwargs):
    return plots.build_full_background_subtraction_e8_4_payload(
        source, baseline, simc_support=_support(source, baseline),
        parent_closure=_closure(len(baseline["t_edges"]) - 1), **kwargs)


class E84ProductionImpactAuditTests(unittest.TestCase):
    def setUp(self):
        patcher = patch.object(plots, "_e8_4_read_candidate_parent_closure", return_value=_closure())
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_valid_consumer_retains_canonical_geometry_deltas_and_detaches_histograms(self):
        source = _f6_3_source(reverse=True)
        before = deepcopy(source["children"][0]["MM_A"].contents)
        payload = _build(source, _e8_2_payload())

        self.assertTrue(payload["available"])
        self.assertEqual(payload["setting_id"], "Left-lowe")
        self.assertEqual(payload["epsilon"], "low")
        self.assertEqual(payload["lambda_window"], [1.08, 1.18])
        self.assertEqual(payload["mm_edges"], [1.10, 1.13, 1.16])
        self.assertNotEqual(payload["mm_edges"], _e8_2_payload()["mm_edges"])
        children = payload["per_t"][0]["children"]
        self.assertEqual([child["phi_index"] for child in children], [0, 1])
        self.assertAlmostEqual(children[0]["delta_y"], 0.6)
        self.assertAlmostEqual(children[0]["delta_y_over_y0"], 0.2)
        self.assertEqual(children[1]["delta_y"], 0.0)
        self.assertEqual(children[1]["delta_y_over_y0"], 0.0)
        self.assertIsNot(children[0]["histograms"]["MM_A"], source["children"][1]["MM_A"])
        self.assertEqual(source["children"][0]["MM_A"].contents, before)

    def test_semantic_epsilon_maps_to_f6_3_canonical_setting_identity(self):
        for setting, epsilon, selected_setting_id in (
            ("Left", "low", "Left-lowe"),
            ("Center", "high", "Center-highe"),
        ):
            with self.subTest(selected_setting_id=selected_setting_id):
                payload = _build(
                    _f6_3_source(selected_setting_id=selected_setting_id),
                    _e8_2_payload(setting=setting, epsilon=epsilon),
                )
                self.assertTrue(payload["available"])
                self.assertEqual(payload["setting"], setting)
                self.assertEqual(payload["epsilon"], epsilon)
                self.assertEqual(payload["setting_id"], selected_setting_id)

    def test_all_one_branch_and_display_differences_are_zero_without_source_mutation(self):
        source = _f6_3_source(y0=(3.0, 4.0), ya=(3.0, 4.0))
        for child in source["children"]:
            child["B_pi_A"].contents[:] = child["B_pi_0"].contents
            child["MM_A"].contents[:] = child["MM_0"].contents
        source_contents = list(source["children"][0]["B_pi_A"].contents)
        payload = _build(source, _e8_2_payload())
        child = payload["per_t"][0]["children"][0]
        difference = plots._e8_4_signed_difference(
            _FakeRoot, child["histograms"]["B_pi_A"], child["histograms"]["B_pi_0"],
            "fixture", "fixture", _FakeRoot.kBlue,
        )
        self.assertEqual(difference.contents, [0.0, 0.0])
        self.assertEqual(child["delta_y"], 0.0)
        self.assertEqual(child["delta_y_over_y0"], 0.0)
        self.assertEqual(source["children"][0]["B_pi_A"].contents, source_contents)

    def test_fail_closed_source_identity_geometry_baseline_and_scalar_contracts(self):
        cases = []
        cases.append((None, _e8_2_payload(), "f6_3_source_schema_invalid"))
        unavailable = _f6_3_source(); unavailable["available"] = False; unavailable["reason"] = "literal branch reason"
        cases.append((unavailable, _e8_2_payload(), "literal branch reason"))
        wrong_schema = _f6_3_source(); wrong_schema["schema_version"] = "bad"
        cases.append((wrong_schema, _e8_2_payload(), "f6_3_source_schema_invalid"))
        wrong_role = _f6_3_source(); wrong_role["branch_role"] = "public"
        cases.append((wrong_role, _e8_2_payload(), "e8_4_branch_role_invalid"))
        wrong_flag = _f6_3_source(); wrong_flag["production_promotion_performed"] = True
        cases.append((wrong_flag, _e8_2_payload(), "e8_4_nonproduction_flags_invalid"))
        no_parity = _f6_3_source(); no_parity["authority"]["live_cache_parity_passed"] = False
        cases.append((no_parity, _e8_2_payload(), "e8_4_live_cache_parity_not_passed"))
        unsupported_epsilon = _e8_2_payload(epsilon="lowe")
        cases.append((_f6_3_source(), unsupported_epsilon, "e8_4_baseline_epsilon_identity_invalid"))
        wrong_setting = _f6_3_source(); wrong_setting["selected_setting_id"] = "Right-highe"
        cases.append((wrong_setting, _e8_2_payload(), "e8_4_setting_identity_mismatch"))
        wrong_window = _f6_3_source(); wrong_window["lambda_integration_window"] = [1.07, 1.18]
        cases.append((wrong_window, _e8_2_payload(), "e8_4_lambda_window_mismatch"))
        malformed_window_integer = _f6_3_source(); malformed_window_integer["lambda_integration_window"] = 1
        cases.append((malformed_window_integer, _e8_2_payload(), "e8_4_lambda_window_invalid"))
        malformed_window_mapping = _f6_3_source(); malformed_window_mapping["lambda_integration_window"] = {"low": 1.08, "high": 1.18}
        cases.append((malformed_window_mapping, _e8_2_payload(), "e8_4_lambda_window_invalid"))
        missing_child = _f6_3_source(); missing_child["children"].pop()
        cases.append((missing_child, _e8_2_payload(), "e8_4_child_inventory_invalid"))
        malformed_children_integer = _f6_3_source(); malformed_children_integer["children"] = 1
        cases.append((malformed_children_integer, _e8_2_payload(), "e8_4_child_inventory_invalid"))
        malformed_children_mapping = _f6_3_source(); malformed_children_mapping["children"] = {"t0phi0": {}}
        cases.append((malformed_children_mapping, _e8_2_payload(), "e8_4_child_inventory_invalid"))
        duplicate = _f6_3_source(); duplicate["children"][1]["phi_index"] = 0
        cases.append((duplicate, _e8_2_payload(), "e8_4_child_inventory_invalid"))
        bad_geometry = _f6_3_source(); bad_geometry["children"][0]["phi_high"] = 1.0
        cases.append((bad_geometry, _e8_2_payload(), "e8_4_child_geometry_invalid"))
        mismatch = _f6_3_source(); mismatch["children"][0]["Y0"] = 9.0
        cases.append((mismatch, _e8_2_payload(), "e8_4_baseline_y0_identity_mismatch"))
        stat_mismatch = _f6_3_source(); stat_mismatch["children"][0]["Y0_statistical_error"] = 0.21
        cases.append((stat_mismatch, _e8_2_payload(), "e8_4_baseline_y0_identity_mismatch"))
        total_mismatch = _f6_3_source(); total_mismatch["children"][0]["Y0_total_error"] = 0.31
        cases.append((total_mismatch, _e8_2_payload(), "e8_4_baseline_y0_identity_mismatch"))
        missing_hist = _f6_3_source(); missing_hist["children"][0].pop("MM_A")
        cases.append((missing_hist, _e8_2_payload(), "e8_4_required_histogram_clone_failed"))
        nonfinite_hist = _f6_3_source(); nonfinite_hist["children"][0]["MM_A"].contents[0] = float("nan")
        cases.append((nonfinite_hist, _e8_2_payload(), "e8_4_required_histogram_nonfinite"))
        nonfinite = _f6_3_source(); nonfinite["children"][0]["YA_total_error"] = float("nan")
        cases.append((nonfinite, _e8_2_payload(), "e8_4_YA_total_error_nonfinite"))
        for source, baseline, reason in cases:
            with self.subTest(reason=reason):
                payload = _build(source, baseline)
                self.assertFalse(payload["available"])
                self.assertEqual(payload["reason"], reason)

    def test_internal_f6_3_analysis_geometry_is_required_and_retained(self):
        narrow_edges = (1.10, 1.13, 1.16)
        wide_edges = (0.70, 1.00, 1.30)
        source = _f6_3_source(t_edges=(0.0, 1.0, 2.0), mm_edges=narrow_edges)
        payload = _build(
            source,
            _e8_2_payload(t_edges=(0.0, 1.0, 2.0), mm_edges=wide_edges),
        )
        self.assertTrue(payload["available"])
        self.assertEqual(payload["mm_edges"], list(narrow_edges))
        self.assertEqual(len(payload["per_t"]), 2)
        self.assertEqual([child["phi_index"] for child in payload["per_t"][1]["children"]], [0, 1])

        within_child = _f6_3_source()
        within_child["children"][0]["B_pi_A"] = _histogram(2.1, edges=(1.10, 1.14, 1.16))
        failed_within_child = _build(
            within_child, _e8_2_payload(),
        )
        self.assertFalse(failed_within_child["available"])
        self.assertEqual(failed_within_child["reason"], "e8_4_histogram_binning_mismatch")

        later_child = _f6_3_source(t_edges=(0.0, 1.0, 2.0))
        for name in ("pion_input", "B_pi_0", "B_pi_A", "MM_0", "MM_A"):
            later_child["children"][-1][name] = _histogram(8.0, edges=(1.11, 1.14, 1.16))
        failed_later_child = _build(
            later_child, _e8_2_payload(t_edges=(0.0, 1.0, 2.0)),
        )
        self.assertFalse(failed_later_child["available"])
        self.assertEqual(failed_later_child["reason"], "e8_4_histogram_binning_mismatch")

    def test_zero_y0_has_an_explicit_undefined_fraction_without_epsilon_fallback(self):
        payload = _build(
            _f6_3_source(y0=(0.0, 4.0), ya=(1.0, 4.0)), _e8_2_payload(y0=(0.0, 4.0)),
        )
        child = payload["per_t"][0]["children"][0]
        self.assertTrue(payload["available"])
        self.assertEqual(child["delta_y"], 1.0)
        self.assertIsNone(child["delta_y_over_y0"])
        self.assertEqual(child["delta_y_over_y0_reason"], "y0_zero_fraction_undefined")

    def test_renderer_pages_have_stable_ids_points_only_fraction_and_no_source_mutation(self):
        source = _f6_3_source(y0=(0.0, 4.0), ya=(1.0, 4.0))
        payload = _build(source, _e8_2_payload(y0=(0.0, 4.0)))
        before = _source_histogram_snapshot(source)
        manifest, failures = [], []
        _FakeCanvas.printed[:] = []; _FakeGraph.options[:] = []; _FakeText.lines[:] = []
        _Histogram.draw_calls[:] = []
        plots._render_full_background_subtraction_e8_4_pages(
            _FakeRoot, "fixture.pdf", payload, manifest, failures,
        )
        self.assertEqual(failures, [])
        self.assertEqual(
            [page["page_id"] for page in manifest],
            [
                "full_background.e8_4.authority",
                "full_background.e8_4.pion_consequence",
                "full_background.e8_4.final_mm",
                "full_background.e8_4.signed_difference",
                "full_background.e8_4.yield_impact",
                "full_background.e8_4.method_a_vs_simc.t1",
                "full_background.e8_4.baseline_method_a_simc.t1",
                "full_background.e8_4.yield_summary.t1",
                "full_background.e8_4.parent_closure",
                "full_background.e8_4.setting_summary",
            ],
        )
        self.assertTrue(_FakeGraph.options)
        self.assertTrue(all(option in ("AP", "P same") for option in _FakeGraph.options))
        pion_draws = {
            name: option for name, option in _Histogram.draw_calls
            if name in {
                "H_e8_4_pion_input_t1_phi1",
                "H_e8_4_bpi0_t1_phi1",
                "H_e8_4_bpia_t1_phi1",
            }
        }
        self.assertEqual(pion_draws, {
            "H_e8_4_pion_input_t1_phi1": "hist e",
            "H_e8_4_bpi0_t1_phi1": "hist e same",
            "H_e8_4_bpia_t1_phi1": "hist e same",
        })
        self.assertEqual(_source_histogram_snapshot(source), before)

    def test_procedure_order_places_e84_after_e83_and_before_handoff_and_renders_unavailable(self):
        order = []
        with patch.object(plots, "_import_root", return_value=object()), patch.object(
            plots, "_render_full_background_subtraction_e8_context_page", side_effect=lambda *_: order.append("e8") or True,
        ), patch.object(plots, "_render_full_background_subtraction_e8_parent_pages", side_effect=lambda *_: order.append("parents")), patch.object(
            plots, "_render_full_background_subtraction_e8_2_pages", side_effect=lambda *_: order.append("e8_2"),
        ), patch.object(plots, "_render_full_background_subtraction_e8_3_pages", side_effect=lambda *_: order.append("e8_3")), patch.object(
            plots, "_render_full_background_subtraction_e8_4_pages", side_effect=lambda *_: order.append("e8_4")), patch.object(
            plots, "_render_full_background_subtraction_e8_handoff_page", side_effect=lambda *_: order.append("handoff") or True,
        ):
            result = plots.render_full_background_subtraction_procedure_pages(
                "fixture.pdf", {"available": False}, {"available": False},
                e8_payload={"available": True, "parents": ({"t_index": 0},)},
                e8_2_payload={"available": True}, e8_3_payload={"available": True},
                e8_4_payload={"available": True},
            )
        self.assertEqual(order, ["e8", "parents", "e8_2", "e8_3", "e8_4", "handoff"])
        self.assertEqual(result["manifest"][-1]["page_id"], "full_background.e8.handoff")
        with patch.object(plots, "_import_root", return_value=object()), patch.object(
            plots, "_render_full_background_subtraction_e8_context_page", return_value=True,
        ), patch.object(plots, "_render_full_background_subtraction_e8_handoff_page", return_value=True), patch.object(
            plots, "_render_full_background_subtraction_e8_4_unavailable_page", return_value=True,
        ):
            unavailable = plots.render_full_background_subtraction_procedure_pages(
                "fixture.pdf", {"available": False}, {"available": False},
                e8_payload={"available": True, "parents": ()},
                e8_4_payload={"available": False, "reason": "literal unavailable"},
            )
        self.assertIn("full_background.e8_4.unavailable", [page["page_id"] for page in unavailable["manifest"]])

    def test_finalizer_passes_e84_payload_through_the_existing_single_transaction(self):
        e8_2 = _e8_2_payload()
        e8_4 = _build(_f6_3_source(), e8_2)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdf_path, manifest_path = root / "preliminary.pdf", root / "preliminary.json"
            pdf_path.write_text("preliminary", encoding="utf-8")
            manifest_path.write_text("preliminary", encoding="utf-8")
            hist = {
                "_xsect_support_simc": _support(_f6_3_source(), e8_2),
                "_e8_2_baseline_stage_source": {"fixture": "delegated"},
                "_f6_3_parallel_method_a_source": _f6_3_source(),
                "_e8_2_full_background_render_state": plots.capture_full_background_subtraction_e8_2_render_state(
                    pdf_path=str(pdf_path), page_manifest_path=str(manifest_path),
                    page_manifest_setting={
                        "kinematic_token": "Q4p4W2p74", "epsilon_filename_token": "lowe",
                        "phi_setting": "Left", "particle_type": "kaon",
                    }, payloads=_render_payloads(),
                ),
            }
            seen = []

            def open_pdf(path):
                Path(path).write_text("replacement", encoding="utf-8")
                return True

            def render(*_args, **kwargs):
                seen.append(kwargs["e8_4_payload"])
                kwargs["page_manifest"].append({"page_id": "fixture"})
                return {"failures": []}

            with patch.object(plots, "build_full_background_subtraction_e8_2_payload", return_value=e8_2), patch.object(
                plots, "open_full_background_subtraction_pdf", side_effect=open_pdf,
            ), patch.object(plots, "close_full_background_subtraction_pdf", return_value=True), patch.object(
                plots, "render_full_background_subtraction_procedure_pages", side_effect=render,
            ):
                status = plots.finalize_full_background_subtraction_e8_2(hist)
        self.assertEqual(len(seen), 1)
        self.assertTrue(seen[0]["available"])
        self.assertTrue(status["e8_4_available"])
        self.assertIsNone(status["e8_4_reason"])

    def test_unavailable_e84_keeps_successful_baseline_finalization_available(self):
        e8_2 = _e8_2_payload()
        unavailable = _f6_3_source()
        unavailable["available"] = False
        unavailable["reason"] = "literal sidecar unavailable"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdf_path, manifest_path = root / "preliminary.pdf", root / "preliminary.json"
            pdf_path.write_text("preliminary", encoding="utf-8")
            manifest_path.write_text("preliminary", encoding="utf-8")
            hist = {
                "_xsect_support_simc": _support(_f6_3_source(), e8_2),
                "_e8_2_baseline_stage_source": {"fixture": "delegated"},
                "_f6_3_parallel_method_a_source": unavailable,
                "_e8_2_full_background_render_state": plots.capture_full_background_subtraction_e8_2_render_state(
                    pdf_path=str(pdf_path), page_manifest_path=str(manifest_path),
                    page_manifest_setting={
                        "kinematic_token": "Q4p4W2p74", "epsilon_filename_token": "lowe",
                        "phi_setting": "Left", "particle_type": "kaon",
                    }, payloads=_render_payloads(),
                ),
            }

            def open_pdf(path):
                Path(path).write_text("replacement", encoding="utf-8")
                return True

            def renderer(*_args, **kwargs):
                self.assertFalse(kwargs["e8_4_payload"]["available"])
                self.assertEqual(kwargs["e8_4_payload"]["reason"], "literal sidecar unavailable")
                kwargs["page_manifest"].append({"page_id": "fixture"})
                return {"failures": []}

            with patch.object(plots, "build_full_background_subtraction_e8_2_payload", return_value=e8_2), patch.object(
                plots, "open_full_background_subtraction_pdf", side_effect=open_pdf,
            ), patch.object(plots, "close_full_background_subtraction_pdf", return_value=True), patch.object(
                plots, "render_full_background_subtraction_procedure_pages", side_effect=renderer,
            ):
                status = plots.finalize_full_background_subtraction_e8_2(hist)
            self.assertEqual(pdf_path.read_text(encoding="utf-8"), "replacement")
            self.assertTrue(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(status["status"], "available")
        self.assertFalse(status["e8_4_available"])
        self.assertEqual(status["e8_4_reason"], "literal sidecar unavailable")
        self.assertNotIn("_e8_2_full_background_render_state", hist)

    def test_malformed_e84_sidecar_keeps_successful_baseline_finalization_available(self):
        e8_2 = _e8_2_payload()
        malformed = _f6_3_source()
        lambda_container = {"low": 1.08, "high": 1.18}
        malformed["lambda_integration_window"] = lambda_container
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdf_path, manifest_path = root / "preliminary.pdf", root / "preliminary.json"
            pdf_path.write_text("preliminary", encoding="utf-8")
            manifest_path.write_text("preliminary", encoding="utf-8")
            hist = {
                "_xsect_support_simc": _support(_f6_3_source(), e8_2),
                "_e8_2_baseline_stage_source": {"fixture": "delegated"},
                "_f6_3_parallel_method_a_source": malformed,
                "_e8_2_full_background_render_state": plots.capture_full_background_subtraction_e8_2_render_state(
                    pdf_path=str(pdf_path), page_manifest_path=str(manifest_path),
                    page_manifest_setting={
                        "kinematic_token": "Q4p4W2p74", "epsilon_filename_token": "lowe",
                        "phi_setting": "Left", "particle_type": "kaon",
                    }, payloads=_render_payloads(),
                ),
            }

            def open_pdf(path):
                Path(path).write_text("replacement", encoding="utf-8")
                return True

            def renderer(*_args, **kwargs):
                self.assertIs(kwargs["e8_2_payload"], e8_2)
                self.assertFalse(kwargs["e8_4_payload"]["available"])
                self.assertEqual(kwargs["e8_4_payload"]["reason"], "e8_4_lambda_window_invalid")
                kwargs["page_manifest"].append({"page_id": "fixture"})
                return {"failures": []}

            with patch.object(plots, "build_full_background_subtraction_e8_2_payload", return_value=e8_2), patch.object(
                plots, "open_full_background_subtraction_pdf", side_effect=open_pdf,
            ), patch.object(plots, "close_full_background_subtraction_pdf", return_value=True), patch.object(
                plots, "render_full_background_subtraction_procedure_pages", side_effect=renderer,
            ):
                status = plots.finalize_full_background_subtraction_e8_2(hist)
            self.assertEqual(pdf_path.read_text(encoding="utf-8"), "replacement")
            self.assertTrue(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(status["status"], "available")
        self.assertFalse(status["e8_4_available"])
        self.assertEqual(status["e8_4_reason"], "e8_4_lambda_window_invalid")
        self.assertIs(malformed["lambda_integration_window"], lambda_container)
        self.assertNotIn("_e8_2_full_background_render_state", hist)

    def test_unavailable_e84_is_rendered_without_replacing_baseline_or_weakening_pair_recovery(self):
        e8_2 = _e8_2_payload()
        unavailable = _f6_3_source(); unavailable["available"] = False; unavailable["reason"] = "literal sidecar unavailable"
        preserved_pdf = None
        preserved_manifest = None
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdf_path, manifest_path = root / "preliminary.pdf", root / "preliminary.json"
            pdf_path.write_text("preliminary", encoding="utf-8")
            manifest_path.write_text("preliminary", encoding="utf-8")
            hist = {
                "_xsect_support_simc": _support(_f6_3_source(), e8_2),
                "_e8_2_baseline_stage_source": {"fixture": "delegated"},
                "_f6_3_parallel_method_a_source": unavailable,
                "_e8_2_full_background_render_state": plots.capture_full_background_subtraction_e8_2_render_state(
                    pdf_path=str(pdf_path), page_manifest_path=str(manifest_path),
                    page_manifest_setting={"phi_setting": "Left"}, payloads=_render_payloads(),
                ),
            }

            def open_pdf(path):
                Path(path).write_text("replacement", encoding="utf-8")
                return True

            def renderer(*_args, **kwargs):
                self.assertFalse(kwargs["e8_4_payload"]["available"])
                self.assertEqual(kwargs["e8_4_payload"]["reason"], "literal sidecar unavailable")
                return {"failures": ["fixture renderer failure"]}

            with patch.object(plots, "build_full_background_subtraction_e8_2_payload", return_value=e8_2), patch.object(
                plots, "open_full_background_subtraction_pdf", side_effect=open_pdf,
            ), patch.object(plots, "close_full_background_subtraction_pdf", return_value=True), patch.object(
                plots, "render_full_background_subtraction_procedure_pages", side_effect=renderer,
            ):
                status = plots.finalize_full_background_subtraction_e8_2(hist)
            preserved_pdf = pdf_path.read_text(encoding="utf-8")
            preserved_manifest = manifest_path.read_text(encoding="utf-8")
        self.assertEqual(status["status"], "unavailable")
        self.assertTrue(status["preliminary_artifact_preserved"])
        self.assertFalse(status["e8_4_available"])
        self.assertEqual(status["e8_4_reason"], "literal sidecar unavailable")
        self.assertEqual(preserved_pdf, "preliminary")
        self.assertEqual(preserved_manifest, "preliminary")
        self.assertIn("_e8_2_full_background_render_state", hist)

    def test_static_consumer_guards_exclude_producer_and_shift_uncertainty_work(self):
        source = inspect.getsource(plots.build_full_background_subtraction_e8_4_payload)
        renderer = inspect.getsource(plots._render_full_background_subtraction_e8_4_pages)
        for forbidden in (
            "calculate_yield_data", "bin_data", "fill_simc_shape",
            "bg_fit", "Normalize", ".Scale(", ".Rebin(",
        ):
            self.assertNotIn(forbidden, source + renderer)
        payload = _build(_f6_3_source(), _e8_2_payload())
        self.assertFalse(any("delta_y" in key and "error" in key for key in payload["per_t"][0]["children"][0]))


class E84ShareablePagesTests(unittest.TestCase):
    def fixtures(self):
        geometry = {"t_edges": (0.0, 0.25, 0.5, 0.75), "phi_edges": tuple(range(-180, 181, 40))}
        source = _f6_3_source(y0=(0.0, 4.0), ya=(1.0, 3.0), **geometry)
        baseline = _e8_2_payload(y0=(0.0, 4.0), **geometry)
        return source, baseline, _support(source, baseline), _closure(3)

    def build(self, source, baseline, support, closure):
        return plots.build_full_background_subtraction_e8_4_payload(
            source, baseline, simc_support=support, parent_closure=closure)

    def test_all_27_cells_overlays_stored_summary_closure_and_input_immutability(self):
        source, baseline, support, closure = self.fixtures()
        before = _source_histogram_snapshot(source)
        simc_before = [(list(hist.contents), list(hist.errors), list(hist.edges), hist.directory)
                       for row in support["mm"] for hist in row]
        payload = self.build(source, baseline, support, closure)
        self.assertTrue(payload["available"], payload.get("reason"))
        for group in payload["per_t"]:
            self.assertEqual(len(group["children"]), 9)
            for child in group["children"]:
                self.assertIsNot(child["histograms"]["SIMC"], support["mm"][child["t_index"]][child["phi_index"]])
        _FakeLegend.entries[:] = []; _FakeGraph.points[:] = []; _FakeText.lines[:] = []; _Histogram.draw_calls[:] = []
        manifest, failures = [], []
        plots._render_full_background_subtraction_e8_4_pages(_FakeRoot, "fixture.pdf", payload, manifest, failures)
        self.assertEqual(failures, [])
        from testing.run_e8_4_fix5_left_lowe_plot_gate import NEW_PAGE_IDS, OLD_PAGES
        for page_id in NEW_PAGE_IDS:
            self.assertEqual(sum(page["page_id"] == page_id for page in manifest), 1)
        for page in manifest:
            if "simc" in page["page_id"]:
                self.assertTrue(page["available"])
                self.assertTrue(page["simc_absolute_comparison_available"])
                self.assertIsNone(page["reason"])
        self.assertTrue(set(OLD_PAGES) <= {(page["page_id"], page["scope"]) for page in manifest})
        self.assertEqual(len(manifest), 24)
        self.assertEqual(_FakeLegend.entries.count("SIMC"), 54)
        self.assertEqual(_FakeLegend.entries.count("Method A data"), 54)
        self.assertEqual(_FakeLegend.entries.count("Baseline data"), 27)
        self.assertIn("parent-normalization", inspect.getsource(plots._e8_4_render_parent_closure_page).lower().replace(" ", "-"))
        self.assertTrue(any("undefined" in text for text in _FakeText.lines))
        self.assertTrue(any(-1.0 in ys for _, ys in _FakeGraph.points))
        self.assertEqual(_source_histogram_snapshot(source), before)
        self.assertEqual([(list(hist.contents), list(hist.errors), list(hist.edges), hist.directory)
                          for row in support["mm"] for hist in row], simc_before)
        self.assertEqual(payload["parent_closure"][0]["baseline_parent_sum"], closure[0]["baseline_parent_sum"])
        self.assertEqual(payload["lambda_window"], source["lambda_integration_window"])

    def test_simc_missing_dimension_axis_nonfinite_and_setting_wide_fallback_fail_closed(self):
        for mutation, reason in (
            (lambda support: support.clear(), "e8_4_simc_mm_support_missing"),
            (lambda support: support.update(mm=support["mm"][:2]), "e8_4_simc_child_geometry_mismatch"),
            (lambda support: support["mm"][1].pop(), "e8_4_simc_child_geometry_mismatch"),
            (lambda support: support["mm"][1].__setitem__(3, None), "e8_4_simc_child_missing"),
            (lambda support: setattr(support["mm"][0][0], "edges", [1.1, 1.14, 1.16]), "e8_4_simc_mm_binning_mismatch"),
            (lambda support: support["mm"][0][0].contents.__setitem__(0, float("nan")), "e8_4_simc_child_nonfinite"),
            (lambda support: support["mm"][0][0].errors.__setitem__(0, float("inf")), "e8_4_simc_child_nonfinite"),
            (lambda support: support.update(mm=_histogram(1.0)), "e8_4_simc_mm_support_missing"),
        ):
            with self.subTest(reason=reason):
                source, baseline, support, closure = self.fixtures()
                mutation(support)
                result = self.build(source, baseline, support, closure)
                self.assertFalse(result["available"])
                self.assertEqual(result["reason"], reason)

    def test_missing_yield_and_method_b_contamination_are_unavailable(self):
        source, baseline, support, closure = self.fixtures()
        source["children"][0].pop("YA")
        self.assertFalse(self.build(source, baseline, support, closure)["available"])
        source, baseline, support, closure = self.fixtures()
        source["method_b_numerical_dependency"] = True
        self.assertFalse(self.build(source, baseline, support, closure)["available"])

    def test_empty_pad_and_new_page_renderer_failure_are_explicit(self):
        source, baseline, support, closure = self.fixtures()
        source["children"][0]["MM_A"].contents[:] = [0.0, 0.0]
        source["children"][0]["MM_A"].errors[:] = [0.0, 0.0]
        source["children"][0]["YA"] = 0.0
        source["children"][0]["B_pi_A"].contents[:] = source["children"][0]["pion_input"].contents
        _refresh_identity(source)
        payload = self.build(source, baseline, support, closure)
        _FakeText.lines[:] = []
        self.assertTrue(plots._e8_4_render_method_a_simc_page(_FakeRoot, "fixture.pdf", payload, payload["per_t"][0]))
        self.assertIn("EMPTY", _FakeText.lines)
        manifest, failures = [], []
        with patch.object(plots, "_e8_4_render_method_a_simc_page", return_value=False):
            plots._render_full_background_subtraction_e8_4_pages(_FakeRoot, "fixture.pdf", payload, manifest, failures)
        self.assertEqual(len(failures), 3)
        self.assertTrue(all("method_a_vs_simc.t" in reason for reason in failures))

    def test_summary_never_integrates_or_rescales_and_uses_stored_values(self):
        functions = (plots._e8_4_render_simc_page, plots._e8_4_render_yield_summary_page,
                     plots.build_full_background_subtraction_e8_4_payload, plots._e8_4_render_parent_closure_page)
        code = "\n".join(inspect.getsource(function) for function in functions)
        for forbidden in (".Scale(", ".Integral(", ".Fill(", "TFile", ".Normalize(", "find_yield_simc"):
            self.assertNotIn(forbidden, code)
        source, baseline, support, closure = self.fixtures()
        payload = self.build(source, baseline, support, closure)
        child = payload["per_t"][0]["children"][1]
        child["delta_y"] = -123.5
        _FakeGraph.points[:] = []
        self.assertTrue(plots._e8_4_render_yield_summary_page(_FakeRoot, "fixture.pdf", payload, payload["per_t"][0]))
        self.assertTrue(any(-123.5 in values for _, values in _FakeGraph.points))

    def test_candidate_parent_closure_is_hash_pinned_stored_input(self):
        source, baseline, support, closure = self.fixtures()
        from pion_hgcer_method_a_parallel_full_procedure import accepted_f6_3_artifact_paths
        setting = {"phi_setting": "Left", "epsilon_filename_token": "lowe", "kinematic_token": "Q4p4W2p74", "particle_type": "kaon"}
        raw = json.dumps({"correction": {"fingerprint": "fixture", "parents": [dict(parent, setting=setting) for parent in closure]}}).encode()
        source["authority"].update(accepted_f4_source_file_sha256=hashlib.sha256(raw).hexdigest(), accepted_f4_correction_fingerprint="fixture")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(accepted_f6_3_artifact_paths(directory, "Q4p4W2p74")["f4"])
            path.write_bytes(raw)
            result = plots._e8_4_read_candidate_parent_closure(source, Path(directory) / "procedure.pdf")
            self.assertEqual(len(result), 3)
            self.assertEqual(result[0]["baseline_parent_sum"], closure[0]["baseline_parent_sum"])
            path.write_bytes(raw + b" ")
            with self.assertRaisesRegex(plots._E8PayloadError, "identity_mismatch"):
                plots._e8_4_read_candidate_parent_closure(source, Path(directory) / "procedure.pdf")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
