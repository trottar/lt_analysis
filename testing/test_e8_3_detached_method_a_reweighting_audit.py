"""Focused pure-Python E.8.3 frozen-authority reader coverage."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [
    str(REPO_ROOT / "src" / "utility"), str(REPO_ROOT / "src" / "cuts"),
    str(REPO_ROOT / "testing"),
]

import full_background_subtraction_plots as plots
import pion_hgcer_method_a_reweighting_validation as f6_1
import test_pion_hgcer_method_a_reweighting_validation as f6_1_fixtures


def _write_json(directory, name, value):
    path = Path(directory) / name
    raw = json.dumps(value, sort_keys=True, indent=2, allow_nan=False).encode("utf-8")
    path.write_bytes(raw)
    return path, hashlib.sha256(raw).hexdigest()


class _FakeCanvas:
    printed = []

    def __init__(self, *args):
        self.args = args

    def Divide(self, *args):
        return None

    def cd(self, *args):
        return self

    def Print(self, value):
        self.printed.append(value)

    def Close(self):
        return None


class _FakeGraph:
    draw_options = []

    def __init__(self, *args):
        self.args = args

    def SetTitle(self, value):
        self.title = value

    def SetLineColor(self, value):
        self.color = value

    def SetLineWidth(self, value):
        self.width = value

    def SetStats(self, value):
        self.stats = value

    def Draw(self, option):
        self.draw_options.append(option)


class _FakeHistogram:
    def Draw(self, option):
        self.draw_option = option


class _FakePaveText:
    instances = []

    def __init__(self, *args):
        self.args = args
        self.lines = []
        self.instances.append(self)

    def SetFillStyle(self, value):
        self.fill_style = value

    def SetBorderSize(self, value):
        self.border_size = value

    def SetTextAlign(self, value):
        self.text_align = value

    def SetTextSize(self, value):
        self.text_size = value

    def AddText(self, value):
        self.lines.append(str(value))

    def Draw(self):
        return None


class _FakeRoot:
    TCanvas = _FakeCanvas
    TGraph = _FakeGraph
    TPaveText = _FakePaveText
    kBlack = 1
    kBlue = 4
    kMagenta = 6


class DetachedMethodAReweightingAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (
            cls.artifacts, cls.hashes, cls.f3, cls.f4, cls.f5, cls.authority,
            cls.f3_authority, cls.f4_authority,
        ) = f6_1_fixtures._chain()
        cls.f6_1 = f6_1.build_pion_hgcer_method_a_reweighting_validation_artifact(
            cls.artifacts, cls.f3, cls.f4, cls.f5,
            f1_input_file_hashes=cls.hashes,
            f3_input_file_sha256=f6_1_fixtures.F3_SHA,
            f4_input_file_sha256=f6_1_fixtures.F4_SHA,
            f5_input_file_sha256=f6_1_fixtures.F5_SHA,
            input_paths={"fixture": "only"},
            accepted_runtime_authority_by_kinematic=cls.authority,
            accepted_f4_runtime_authority_by_kinematic=cls.f4_authority,
            accepted_f3_runtime_authority_by_kinematic=cls.f3_authority,
        )

    def _authority(self):
        return {
            "f4_source_file_sha256": self.f5["propagation"]["f4_source_file_sha256"],
            "f4_correction_fingerprint": self.f4["correction"]["fingerprint"],
            "f4_artifact_fingerprint": self.f4["artifact_fingerprint"],
            "f5_propagation_fingerprint": self.f5["propagation"]["fingerprint"],
            "f5_artifact_fingerprint": self.f5["artifact_fingerprint"],
            "f6_1_artifact_fingerprint": self.f6_1["artifact_fingerprint"],
            "f6_1_validation_fingerprint": self.f6_1["validation"]["fingerprint"],
        }

    @staticmethod
    def _f6_2_payload(setting_id="Left-lowe"):
        return {
            "available": True,
            "setting_id": setting_id,
            "json_path": "accepted-f6-2.json",
            "input_sha256": plots.E8_F6_2_INPUT_SHA256,
            "artifact_fingerprint": plots.E8_F6_2_ARTIFACT_FINGERPRINT,
            "validation_fingerprint": plots.E8_F6_2_VALIDATION_FINGERPRINT,
        }

    def _payload(self, directory, *, setting_id="Left-lowe", f4=None, f5=None, f6_1_artifact=None, **kwargs):
        f4_path, f4_sha = _write_json(directory, "f4.json", self.f4 if f4 is None else f4)
        f5_path, f5_sha = _write_json(directory, "f5.json", self.f5 if f5 is None else f5)
        f6_1_path, f6_1_sha = _write_json(directory, "f6_1.json", self.f6_1 if f6_1_artifact is None else f6_1_artifact)
        expected = {
            "expected_f4_sha256": f4_sha, "expected_f5_sha256": f5_sha,
            "expected_f6_1_sha256": f6_1_sha,
        }
        expected.update({key: kwargs.pop(key) for key in tuple(expected) if key in kwargs})
        return plots.build_full_background_subtraction_e8_3_payload(
            f4_path, f5_path, f6_1_path, self._f6_2_payload(setting_id),
            setting_id=setting_id, authority=self._authority(), **expected, **kwargs,
        )

    def _tphi_lines(self, parent, *, setting_id="Left-lowe"):
        _FakeCanvas.printed[:] = []
        _FakePaveText.instances[:] = []
        self.assertTrue(plots._e8_3_render_tphi_page(
            _FakeRoot, "ignored.pdf", {"setting_id": setting_id}, parent,
        ))
        self.assertEqual(_FakeCanvas.printed, ["ignored.pdf"])
        return _FakePaveText.instances[-1].lines

    def test_exact_selected_persisted_authority_passes_through_without_child_normalization(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory)
        self.assertTrue(payload["available"])
        self.assertEqual(payload["setting_id"], "Left-lowe")
        self.assertEqual(len(payload["parents"]), 3)
        source_setting = next(
            item for item in self.f5["propagation"]["setting_templates"]
            if item["setting_id"] == "Left-lowe"
        )
        for index, parent in enumerate(payload["parents"]):
            self.assertEqual(parent["f5_baseline_signed_contents"], source_setting["baseline_signed_contents"][index])
            self.assertEqual(parent["f5_method_a_signed_contents"], source_setting["adjusted_signed_contents"][index])
            self.assertEqual(parent["f5_signed_delta_contents"], source_setting["signed_delta_contents"][index])
            self.assertEqual(parent["f5_event_counts"], source_setting["event_counts"][index])
        self.assertFalse(payload["canonical_child_renormalization_performed"])
        self.assertFalse(payload["method_b_numerical_dependency"])
        self.assertFalse(payload["empirical_residual_dependency"])

    def test_wrong_raw_sha_and_wrong_fingerprint_fail_closed_before_any_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory, expected_f5_sha256="0" * 64)
            self.assertFalse(payload["available"])
            self.assertIn("f5_input_sha256_mismatch", payload["reason"])
            wrong_authority = self._authority()
            wrong_authority["f5_propagation_fingerprint"] = "0" * 64
            f4_path, f4_sha = _write_json(directory, "f4.json", self.f4)
            f5_path, f5_sha = _write_json(directory, "f5.json", self.f5)
            f6_1_path, f6_1_sha = _write_json(directory, "f6_1.json", self.f6_1)
            payload = plots.build_full_background_subtraction_e8_3_payload(
                f4_path, f5_path, f6_1_path, self._f6_2_payload(),
                expected_f4_sha256=f4_sha, expected_f5_sha256=f5_sha,
                expected_f6_1_sha256=f6_1_sha, authority=wrong_authority,
            )
            self.assertFalse(payload["available"])
            self.assertIn("f5_propagation_authority_invalid", payload["reason"])

    def test_f4_raw_sha_correction_and_artifact_authority_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory, expected_f4_sha256="0" * 64)
            self.assertFalse(payload["available"])
            self.assertIn("f4_input_sha256_mismatch", payload["reason"])
            wrong_correction = deepcopy(self.f4)
            wrong_correction["correction"]["fingerprint"] = "0" * 64
            payload = self._payload(directory, f4=wrong_correction)
            self.assertFalse(payload["available"])
            self.assertIn("f4_correction_authority_invalid", payload["reason"])
            wrong_artifact = deepcopy(self.f4)
            wrong_artifact["artifact_fingerprint"] = "0" * 64
            payload = self._payload(directory, f4=wrong_artifact)
            self.assertFalse(payload["available"])
            self.assertIn("f4_artifact_fingerprint_invalid", payload["reason"])

    def test_f5_raw_sha_propagation_and_artifact_authority_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory, expected_f5_sha256="0" * 64)
            self.assertFalse(payload["available"])
            self.assertIn("f5_input_sha256_mismatch", payload["reason"])
            wrong_propagation = deepcopy(self.f5)
            wrong_propagation["propagation"]["fingerprint"] = "0" * 64
            payload = self._payload(directory, f5=wrong_propagation)
            self.assertFalse(payload["available"])
            self.assertIn("f5_propagation_authority_invalid", payload["reason"])
            wrong_artifact = deepcopy(self.f5)
            wrong_artifact["artifact_fingerprint"] = "0" * 64
            payload = self._payload(directory, f5=wrong_artifact)
            self.assertFalse(payload["available"])
            self.assertIn("f5_artifact_fingerprint_invalid", payload["reason"])

    def test_f6_1_raw_sha_artifact_and_validation_authority_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory, expected_f6_1_sha256="0" * 64)
            self.assertFalse(payload["available"])
            self.assertIn("f6_1_input_sha256_mismatch", payload["reason"])
            wrong_artifact = deepcopy(self.f6_1)
            wrong_artifact["artifact_fingerprint"] = "0" * 64
            payload = self._payload(directory, f6_1_artifact=wrong_artifact)
            self.assertFalse(payload["available"])
            self.assertIn("f6_1_artifact_fingerprint_invalid", payload["reason"])
            wrong_validation = deepcopy(self.f6_1)
            wrong_validation["validation"]["fingerprint"] = "0" * 64
            payload = self._payload(directory, f6_1_artifact=wrong_validation)
            self.assertFalse(payload["available"])
            self.assertIn("f6_1_validation_authority_invalid", payload["reason"])

    def test_wrong_kinematic_token_fails_closed_before_reading_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory, kinematic_token="Q3p0W2p32")
        self.assertFalse(payload["available"])
        self.assertEqual(payload["reason"], "unsupported_e8_3_current_setting")

    def test_nonfinite_and_malformed_canonical_geometry_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            bad = deepcopy(self.f5)
            bad["propagation"]["setting_templates"][0]["baseline_signed_contents"][0][0] = float("nan")
            f4_path, f4_sha = _write_json(directory, "f4.json", self.f4)
            f6_1_path, f6_1_sha = _write_json(directory, "f6_1.json", self.f6_1)
            f5_path = Path(directory) / "f5-nonfinite.json"
            raw = json.dumps(bad, sort_keys=True, allow_nan=True).encode("utf-8")
            f5_path.write_bytes(raw)
            good = plots.build_full_background_subtraction_e8_3_payload(
                f4_path, f5_path, f6_1_path, self._f6_2_payload(),
                expected_f4_sha256=f4_sha,
                expected_f5_sha256=hashlib.sha256(raw).hexdigest(),
                expected_f6_1_sha256=f6_1_sha, authority=self._authority(),
            )
            self.assertFalse(good["available"])
            self.assertIn("f5_json_invalid", good["reason"])
            malformed = deepcopy(self.f5)
            malformed["propagation"]["setting_templates"][0]["baseline_signed_contents"][0].pop()
            payload = self._payload(directory, f5=malformed)
            self.assertFalse(payload["available"])
            self.assertIn("f5_baseline_signed_contents_shape_invalid", payload["reason"])

    def test_missing_and_duplicate_canonical_f5_settings_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = deepcopy(self.f5)
            missing["propagation"]["setting_templates"].pop()
            payload = self._payload(directory, f5=missing)
            self.assertFalse(payload["available"])
            self.assertIn("f5_setting_inventory_invalid", payload["reason"])
            duplicate = deepcopy(self.f5)
            duplicate["propagation"]["setting_templates"][1] = deepcopy(
                duplicate["propagation"]["setting_templates"][0]
            )
            payload = self._payload(directory, f5=duplicate)
            self.assertFalse(payload["available"])
            self.assertIn("f5_setting_identity_invalid", payload["reason"])

    def test_missing_and_duplicate_canonical_f4_parents_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = deepcopy(self.f4)
            missing["correction"]["parents"].pop()
            payload = self._payload(directory, f4=missing)
            self.assertFalse(payload["available"])
            self.assertIn("f4_parent_inventory_invalid", payload["reason"])
            duplicate = deepcopy(self.f4)
            duplicate["correction"]["parents"][1] = deepcopy(
                duplicate["correction"]["parents"][0]
            )
            payload = self._payload(directory, f4=duplicate)
            self.assertFalse(payload["available"])
            self.assertIn("f4_parent_geometry_invalid", payload["reason"])

    def test_missing_and_duplicate_canonical_f6_1_parents_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = deepcopy(self.f6_1)
            missing["validation"]["parents"].pop()
            payload = self._payload(directory, f6_1_artifact=missing)
            self.assertFalse(payload["available"])
            self.assertIn("f6_1_parent_inventory_invalid", payload["reason"])
            duplicate = deepcopy(self.f6_1)
            duplicate["validation"]["parents"][1] = deepcopy(
                duplicate["validation"]["parents"][0]
            )
            payload = self._payload(directory, f6_1_artifact=duplicate)
            self.assertFalse(payload["available"])
            self.assertIn("f6_1_parent_geometry_invalid", payload["reason"])

    def test_malformed_f6_1_mm_edges_and_content_lengths_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            malformed_edges = deepcopy(self.f6_1)
            edges = malformed_edges["validation"]["parents"][0]["signed_background"]["analysis_MM"]["edges"]
            edges[1] = edges[0]
            payload = self._payload(directory, f6_1_artifact=malformed_edges)
            self.assertFalse(payload["available"])
            self.assertIn("f6_1_mm_edges_invalid", payload["reason"])
            malformed_contents = deepcopy(self.f6_1)
            malformed_contents["validation"]["parents"][0]["signed_background"]["analysis_MM"]["baseline_signed_contents"].pop()
            payload = self._payload(directory, f6_1_artifact=malformed_contents)
            self.assertFalse(payload["available"])
            self.assertIn("f6_1_mm_content_shape_invalid", payload["reason"])

    def test_current_setting_is_exact_and_f6_1_mm_arrays_have_no_reconstruction(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory, setting_id="Center-highe")
        self.assertTrue(payload["available"])
        self.assertEqual({parent["canonical_t_index"] for parent in payload["parents"]}, {0, 1, 2})
        source = next(
            parent for parent in self.f6_1["validation"]["parents"]
            if parent["setting_id"] == "Center-highe" and parent["canonical_t_index"] == 1
        )["signed_background"]["analysis_MM"]
        persisted = payload["parents"][1]["f6_1_analysis_mm"]
        self.assertEqual(persisted["baseline_signed_contents"], source["baseline_signed_contents"])
        self.assertEqual(persisted["method_a_signed_contents"], source["method_a_signed_contents"])
        self.assertEqual(persisted["signed_delta_contents"], source["signed_delta_contents"])

    def test_ratio_mask_only_contains_finite_nonzero_persisted_baseline_bins(self):
        points = plots._e8_3_ratio_points([0.0, 1.0, 2.0, 3.0], [2.0, 0.0, float("nan")], [4.0, 1.0, 3.0])
        self.assertEqual(points, [(0.5, 2.0)])

    def test_ratio_renderer_is_gap_safe_and_never_joins_masked_bins(self):
        parent = {
            "canonical_t_index": 0,
            "f6_1_analysis_mm": {
                "edges": [0.0, 1.0, 2.0, 3.0],
                "baseline_signed_contents": [2.0, 0.0, 3.0],
                "method_a_signed_contents": [4.0, 7.0, 9.0],
                "signed_delta_contents": [2.0, 7.0, 6.0],
            },
        }
        _FakeCanvas.printed[:] = []
        _FakeGraph.draw_options[:] = []
        with patch.object(plots, "_e8_3_signed_histogram", side_effect=[_FakeHistogram(), _FakeHistogram(), _FakeHistogram()]), \
             patch.object(plots, "_e8_add_text", return_value=object()):
            self.assertTrue(plots._e8_3_render_mm_page(
                _FakeRoot, "ignored.pdf", {"setting_id": "Left-lowe"}, parent,
            ))
        self.assertEqual(_FakeGraph.draw_options, ["AP"])
        self.assertNotIn("L", _FakeGraph.draw_options[0])

    def test_ratio_values_are_unclipped_and_uninterpolated(self):
        points = plots._e8_3_ratio_points(
            [0.0, 1.0, 2.0, 3.0], [1.0, 0.0, -1.0], [1.0e30, 7.0, 3.0],
        )
        self.assertEqual(points, [(0.5, 1.0e30), (2.5, -3.0)])

    def test_persisted_event_counts_drive_empty_status_without_reinterpreting_values(self):
        f5 = deepcopy(self.f5)
        source = next(
            item for item in f5["propagation"]["setting_templates"]
            if item["setting_id"] == "Left-lowe"
        )
        source["event_counts"][0][0] = 0
        source["event_counts"][0][1] = 7
        for name in ("baseline_signed_contents", "adjusted_signed_contents", "signed_delta_contents"):
            source[name][0][1] = 0.0
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory, f5=f5)
        self.assertTrue(payload["available"])
        parent = payload["parents"][0]
        self.assertEqual(parent["f5_event_counts"], source["event_counts"][0])
        self.assertEqual(parent["f5_baseline_signed_contents"], source["baseline_signed_contents"][0])
        self.assertEqual(parent["f5_method_a_signed_contents"], source["adjusted_signed_contents"][0])
        self.assertEqual(parent["f5_signed_delta_contents"], source["signed_delta_contents"][0])
        lines = self._tphi_lines(parent)
        child_lines = [line for line in lines if line.startswith("phi ")]
        self.assertEqual(len(child_lines), 9)
        self.assertIn("EMPTY", child_lines[0])
        self.assertIn("POPULATED n=7", child_lines[1])
        self.assertNotIn("EMPTY", child_lines[1])

    def test_persisted_parent_closure_is_rendered_without_new_metric(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory)
        parent = payload["parents"][0]
        closure = parent["f5_parent_closure"]
        lines = self._tphi_lines(parent)
        rendered = next(line for line in lines if line.startswith("Persisted parent closure:"))
        self.assertIn("{:.8g}".format(closure["baseline_parent_sum"]), rendered)
        self.assertIn("{:.8g}".format(closure["adjusted_parent_sum"]), rendered)
        self.assertIn("{:.8g}".format(abs(closure["closure_residual"])), rendered)
        self.assertIn("No relative closure metric", lines[-1])

    def test_runtime_route_uses_only_outpath_frozen_inputs_and_existing_pdf_lifecycle(self):
        source = (REPO_ROOT / "src" / "cuts" / "rand_sub.py").read_text(encoding="utf-8")
        start = source.index("full_background_subtraction_e8_kinematic =")
        stop = source.index("# The E.8 unavailable page", start)
        route = source[start:stop]
        self.assertIn("build_full_background_subtraction_e8_3_payload(", route)
        self.assertIn("OUTPATH,", route)
        for basename in (
            "method-a-parent-preserving-correction.json",
            "method-a-tphi-propagation.json",
            "method-a-reweighting-validation.json",
        ):
            self.assertIn(basename, route)
        self.assertNotIn("pion_component_fits", route)
        self.assertNotIn("method_b", route.lower())
        self.assertNotIn("calculate_yield", route)
        self.assertNotIn("f6_3", route.lower())
        self.assertNotIn("production_application", route.lower())
        self.assertIn("e8_3_payload=full_background_subtraction_e8_3_payload", source)
        self.assertIn('"e8_3": full_background_subtraction_e8_3_payload', source)

    def test_e8_3_pages_append_after_prior_e8_e8_2_pages_and_unavailable_is_nonfatal(self):
        order = []

        def record(label):
            def renderer(*args):
                order.append(label)
                return True
            return renderer

        def e8_2_renderer(*args):
            order.append("e8_2")

        def e8_3_renderer(*args):
            order.append("e8_3")

        with patch.object(plots, "_import_root", return_value=object()), \
             patch.object(plots, "_render_full_background_subtraction_e8_context_page", record("e8")), \
             patch.object(plots, "_render_full_background_subtraction_e8_2_pages", e8_2_renderer), \
             patch.object(plots, "_render_full_background_subtraction_e8_3_pages", e8_3_renderer), \
             patch.object(plots, "_render_full_background_subtraction_e8_handoff_page", record("handoff")):
            result = plots.render_full_background_subtraction_procedure_pages(
                "ignored.pdf", None, None,
                e8_payload={"available": True, "parents": ()},
                e8_2_payload={"available": True},
                e8_3_payload={"available": True, "parents": ()},
            )
        self.assertEqual(order, ["e8", "e8_2", "e8_3", "handoff"])
        self.assertEqual(result["failures"], [])
        with patch.object(plots, "_import_root", return_value=object()), \
             patch.object(plots, "_render_full_background_subtraction_e8_unavailable_page", return_value=True), \
             patch.object(plots, "_render_full_background_subtraction_e8_3_unavailable_page", return_value=True):
            unavailable = plots.render_full_background_subtraction_procedure_pages(
                "ignored.pdf", None, None, e8_payload={"available": False, "reason": "frozen"},
                e8_3_payload={"available": False, "reason": "missing"},
            )
        self.assertEqual(unavailable["failures"], [])
        self.assertEqual(
            [page["page_id"] for page in unavailable["manifest"]],
            ["full_background.e8.unavailable", "full_background.e8_3.unavailable"],
        )

    def test_unavailable_payload_is_a_presentation_result_not_a_production_exception(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = Path(directory) / "missing.json"
            payload = plots.build_full_background_subtraction_e8_3_payload(
                missing, missing, missing, self._f6_2_payload(),
            )
        self.assertFalse(payload["available"])
        self.assertFalse(payload["production_objects_mutated"])
        self.assertFalse(payload["production_application_performed"])

    def test_e8_3_payload_and_rendering_leave_baseline_sentinels_untouched(self):
        with tempfile.TemporaryDirectory() as directory:
            payload = self._payload(directory)
        baseline_sentinel = {
            "baseline_mm_0": [1.0, -2.0], "yield_0": {"value": 3.0, "error": 0.2},
        }
        payload_before = deepcopy(payload)
        baseline_before = deepcopy(baseline_sentinel)
        _FakeCanvas.printed[:] = []
        _FakeGraph.draw_options[:] = []
        _FakePaveText.instances[:] = []
        with patch.object(plots, "_e8_3_signed_histogram", side_effect=[_FakeHistogram(), _FakeHistogram(), _FakeHistogram()]), \
             patch.object(plots, "_e8_add_text", return_value=object()):
            self.assertTrue(plots._e8_3_render_mm_page(
                _FakeRoot, "ignored.pdf", payload, payload["parents"][0],
            ))
        self._tphi_lines(payload["parents"][0])
        self.assertEqual(payload, payload_before)
        self.assertEqual(baseline_sentinel, baseline_before)
        self.assertFalse(payload["production_objects_mutated"])
        self.assertFalse(payload["production_application_performed"])
        self.assertFalse(payload["method_b_numerical_dependency"])
        self.assertFalse(payload["empirical_residual_dependency"])


if __name__ == "__main__":
    unittest.main()
