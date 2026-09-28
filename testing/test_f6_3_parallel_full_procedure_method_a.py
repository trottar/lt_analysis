"""Deterministic F.6.3 authority and multiplier regressions (no ROOT/farm)."""

from __future__ import annotations

import copy
from contextlib import ExitStack
import json
from pathlib import Path
import inspect
import math
import sys
import unittest
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src" / "cuts"), str(ROOT / "src" / "utility")]

import pion_component_subtraction as pion
import pion_hgcer_method_a_parallel_full_procedure as f63
import pion_hgcer_method_a_parent_preserving_correction as f4
import pion_hgcer_method_a_tphi_propagation as f5
from testing.test_e8_2_baseline_stage_audit import _Histogram as _YieldHistogram
from testing.test_e8_2_baseline_stage_audit import _load_calculate_yield_module
from testing.test_e8_2_baseline_stage_audit import _public_yield_fixture


class _Axis:
    def FindBin(self, value):
        return 1


class _Reference:
    def GetXaxis(self):
        return _Axis()


class _Histogram:
    def __init__(self):
        self.fills = []

    def Fill(self, *values):
        self.fills.append(tuple(float(value) for value in values))


def _row(identity=("prompt", 7), *, coefficient=2.0, weight=3.0):
    return {
        "source_label": identity[0], "entry_index": identity[1],
        "t_index": 0, "phi_index": 0, "analysis_MM": 1.12,
        "analysis_t": 0.5, "signed_source_coefficient": coefficient,
        "baseline_pion_weight_w0": weight,
        "signed_baseline_event_contribution": coefficient * weight,
    }


def _accepted_rows(*rows):
    return {
        (str(row["source_label"]), int(row["entry_index"])): {
            name: row[name]
            for name in (
                "t_index", "phi_index", "analysis_MM", "analysis_t",
                "signed_source_coefficient", "baseline_pion_weight_w0",
                "signed_baseline_event_contribution",
            )
        }
        for row in rows
    }


def _cache():
    return {
        "entry_index": np.asarray([7], dtype=np.int32),
        "adj_MM": np.asarray([1.12]), "adj_t": np.asarray([0.5]),
        "Q2": np.asarray([1.0]), "W": np.asarray([2.0]),
        "epsilon": np.asarray([0.6]), "theta_cm_deg": np.asarray([30.0]),
        "ssxptar": np.asarray([0.0]), "ssyptar": np.asarray([0.0]),
        "hsxptar": np.asarray([0.0]), "hsyptar": np.asarray([0.0]),
        "allcuts": np.asarray([True]), "nommcuts": np.asarray([True]),
        "t_index": np.asarray([0], dtype=np.int32), "phi_index": np.asarray([0], dtype=np.int32),
        "coefficient": np.asarray([2.0]),
        "allcut_bin_index": {(0, 0): np.asarray([0], dtype=np.int32)},
        "nommcut_bin_index": {(0, 0): np.asarray([0], dtype=np.int32)},
    }


class ParallelAuthorityTests(unittest.TestCase):
    def test_live_cache_parity_passes_without_retaining_factor_values(self):
        result = f63.validate_live_cache_parity(
            {("prompt", 7): 1.25}, [_row()], _accepted_rows(_row()),
        )
        self.assertTrue(result["live_cache_parity_passed"])
        self.assertNotIn("1.25", str(result))

    def test_live_cache_identity_and_baseline_mismatches_fail_closed(self):
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "inventory"):
            f63.validate_live_cache_parity(
                {("prompt", 8): 1.0}, [_row()], _accepted_rows(_row()),
            )
        changed = _row(); changed["baseline_pion_weight_w0"] = 4.0
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "baseline_identity"):
            f63.validate_live_cache_parity(
                {("prompt", 7): 1.0}, [changed], _accepted_rows(_row()),
            )

    def test_live_cache_duplicate_nonfinite_and_nonpositive_fail_closed(self):
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "duplicate"):
            f63.validate_live_cache_parity(
                {("prompt", 7): 1.0}, [_row(), _row()], _accepted_rows(_row()),
            )
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "nonfinite_or_nonpositive"):
            f63.validate_live_cache_parity(
                {("prompt", 7): 0.0}, [_row()], _accepted_rows(_row()),
            )

    def test_reconstruction_pairs_only_selected_setting_and_rejects_duplicate_identity(self):
        persisted = {"fingerprint": "f4"}
        rows = {("Left-lowe", 0): [_row()]}
        review = [{"setting_id": "Left-lowe", "canonical_t_index": 0, "correction_factors": np.asarray([1.5])}]
        parsed = [{"setting_id": "Left-lowe"}]
        with mock.patch.object(f63, "_runtime_dependencies", return_value=(f4, f5, None, None, None)), \
             mock.patch.object(f5, "_validate_f4_artifact", return_value=({}, persisted, {"accepted_authority_match": True})), \
             mock.patch.object(f4._f3, "_validate_f1_artifacts", return_value=parsed), \
             mock.patch.object(f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", return_value=(persisted, review)), \
             mock.patch.object(f4, "_raw_application_rows", return_value=rows):
            factors, provenance, accepted_rows = f63.reconstruct_transient_factor_map(
                [{}], {}, {"artifact_fingerprint": "artifact"},
                f1_input_file_hashes={"Left-lowe": "a" * 64}, f3_input_file_sha256="b" * 64,
                f4_input_file_sha256="c" * 64, setting_id="Left-lowe",
            )
        self.assertEqual(factors, {("prompt", 7): 1.5})
        self.assertEqual(accepted_rows, _accepted_rows(_row(weight=3.0, coefficient=2.0)))
        self.assertEqual(provenance["transient_factor_population_count"], 1)
        self.assertNotIn("factors", provenance)
        self.assertNotIn("transient_factor_population_fingerprint", provenance)
        self.assertNotIn("1.5", str(provenance))

    def test_reconstruction_rejects_f4_and_f1_f3_authority_mismatches(self):
        persisted = {"fingerprint": "f4"}
        rows = {("Left-lowe", 0): [_row()]}
        review = [{"setting_id": "Left-lowe", "canonical_t_index": 0, "correction_factors": np.asarray([1.5])}]
        parsed = [{"setting_id": "Left-lowe"}]
        common = {
            "f1_input_file_hashes": {"Left-lowe": "a" * 64},
            "f3_input_file_sha256": "b" * 64,
            "f4_input_file_sha256": "c" * 64,
            "setting_id": "Left-lowe",
        }
        with mock.patch.object(f63, "_runtime_dependencies", return_value=(f4, f5, None, None, None)), \
             mock.patch.object(f5, "_validate_f4_artifact", return_value=({}, persisted, {})), \
             mock.patch.object(f4._f3, "_validate_f1_artifacts", return_value=parsed), \
             mock.patch.object(f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", return_value=({"fingerprint": "other"}, review)), \
             mock.patch.object(f4, "_raw_application_rows", return_value=rows):
            with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "shared_reproduction_mismatch"):
                f63.reconstruct_transient_factor_map([{}], {}, {}, **common)
        with mock.patch.object(f63, "_runtime_dependencies", return_value=(f4, f5, None, None, None)), \
             mock.patch.object(f5, "_validate_f4_artifact", side_effect=RuntimeError("source_sha_mismatch")):
            with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "shared_reproduction_failed"):
                f63.reconstruct_transient_factor_map([{}], {}, {}, **common)
        with mock.patch.object(f63, "_runtime_dependencies", return_value=(f4, f5, None, None, None)), \
             mock.patch.object(f5, "_validate_f4_artifact", return_value=({}, persisted, {})), \
             mock.patch.object(f4._f3, "_validate_f1_artifacts", return_value=parsed):
            bad = dict(common); bad["f1_input_file_hashes"] = {"Right-highe": "a" * 64}
            with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "f1_hash_inventory"):
                f63.reconstruct_transient_factor_map([{}], {}, {}, **bad)

    def test_reconstruction_rejects_length_and_duplicate_application_rows(self):
        persisted = {"fingerprint": "f4"}
        parsed = [{"setting_id": "Left-lowe"}]
        common = {
            "f1_input_file_hashes": {"Left-lowe": "a" * 64},
            "f3_input_file_sha256": "b" * 64,
            "f4_input_file_sha256": "c" * 64,
            "setting_id": "Left-lowe",
        }
        with mock.patch.object(f63, "_runtime_dependencies", return_value=(f4, f5, None, None, None)), \
             mock.patch.object(f5, "_validate_f4_artifact", return_value=({}, persisted, {})), \
             mock.patch.object(f4._f3, "_validate_f1_artifacts", return_value=parsed), \
             mock.patch.object(f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", return_value=(persisted, [{"setting_id": "Left-lowe", "canonical_t_index": 0, "correction_factors": np.asarray([1.0, 2.0])}])), \
             mock.patch.object(f4, "_raw_application_rows", return_value={("Left-lowe", 0): [_row()]}):
            with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "row_alignment"):
                f63.reconstruct_transient_factor_map([{}], {}, {}, **common)
        duplicate_review = [
            {"setting_id": "Left-lowe", "canonical_t_index": 0, "correction_factors": np.asarray([1.0])},
            {"setting_id": "Left-lowe", "canonical_t_index": 1, "correction_factors": np.asarray([1.0])},
        ]
        duplicate_rows = {("Left-lowe", 0): [_row()], ("Left-lowe", 1): [_row()]}
        with mock.patch.object(f63, "_runtime_dependencies", return_value=(f4, f5, None, None, None)), \
             mock.patch.object(f5, "_validate_f4_artifact", return_value=({}, persisted, {})), \
             mock.patch.object(f4._f3, "_validate_f1_artifacts", return_value=parsed), \
             mock.patch.object(f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", return_value=(persisted, duplicate_review)), \
             mock.patch.object(f4, "_raw_application_rows", return_value=duplicate_rows):
            with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "identity_duplicate"):
                f63.reconstruct_transient_factor_map([{}], {}, {}, **common)

    def test_live_cache_requires_every_f1_baseline_field_to_match(self):
        reference = _row()
        for field, changed_value in (
            ("signed_source_coefficient", 3.0),
            ("baseline_pion_weight_w0", 4.0),
            ("signed_baseline_event_contribution", 7.0),
            ("t_index", 1), ("phi_index", 1),
            ("analysis_MM", 1.13), ("analysis_t", 0.51),
        ):
            with self.subTest(field=field):
                changed = _row(); changed[field] = changed_value
                if field == "signed_source_coefficient":
                    changed["signed_baseline_event_contribution"] = changed_value * changed["baseline_pion_weight_w0"]
                with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "(f1_.*_mismatch|baseline_identity_mismatch)"):
                    f63.validate_live_cache_parity(
                        {("prompt", 7): 1.0}, [changed], _accepted_rows(reference),
                    )

    def test_live_cache_uses_accepted_f4_scaled_tolerance(self):
        accepted = _row()
        accepted["analysis_MM"] = 1.0e6
        within = _row()
        within["analysis_MM"] = 1.0e6 + 5.0e-7
        f63.validate_live_cache_parity(
            {("prompt", 7): 1.0}, [within], _accepted_rows(accepted),
        )
        outside = _row()
        outside["analysis_MM"] = 1.0e6 + 2.0e-6
        with self.assertRaisesRegex(f63.MethodAParallelFullProcedureError, "analysis_MM"):
            f63.validate_live_cache_parity(
                {("prompt", 7): 1.0}, [outside], _accepted_rows(accepted),
            )


class TemplateMultiplierTests(unittest.TestCase):
    def _fill(self, multipliers=None):
        templates = {"mm": _Histogram(), "mm_nosub": _Histogram()}
        pion.fill_simc_shape_pion_subtraction_templates(
            templates,
            [{"label": "prompt", "cache_section": _cache(), "coefficient": 2.0, "base_coefficient": 2.0}],
            _Reference(), np.asarray([0.0, 3.0]), {"t_index": 0, "phi_index": 0}, "kaon", "unpolarized",
            method_a_event_multipliers=multipliers,
        )
        return templates

    def test_omitted_and_all_one_multiplier_preserve_baseline(self):
        baseline = self._fill()
        ones = self._fill({("prompt", 7): 1.0})
        self.assertEqual(baseline["mm"].fills, ones["mm"].fills)
        self.assertEqual(baseline["mm_nosub"].fills, ones["mm_nosub"].fills)

    def test_nontrivial_multiplier_scales_allcuts_and_nommcuts_only(self):
        templates = self._fill({("prompt", 7): 1.5})
        self.assertEqual(templates["mm"].fills[0][-1], 9.0)
        self.assertEqual(templates["mm_nosub"].fills[0][-1], 9.0)

    def test_multiplier_lookup_miss_and_invalid_value_are_errors_not_c_one(self):
        with self.assertRaisesRegex(ValueError, "missing"):
            self._fill({})
        with self.assertRaisesRegex(ValueError, "invalid"):
            self._fill({("prompt", 7): math.nan})

    def test_template_extension_neither_normalizes_children_nor_reads_method_b(self):
        text = inspect.getsource(pion.fill_simc_shape_pion_subtraction_templates)
        self.assertNotIn("Scale(", text)
        self.assertNotIn("method_b", text.lower())


class FullParallelBranchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.calculate_yield = _load_calculate_yield_module()

    @staticmethod
    def _clone(histogram, *_args, **kwargs):
        clone = histogram.Clone(kwargs.get("name", "parallel_clone"))
        if kwargs.get("reset"):
            clone.Reset()
        return clone

    def _branch_fixture(self, factor=1.0, *, source_coefficient=2.0, base_coefficient=2.0, accepted_coefficient=2.0):
        calc = self.calculate_yield
        pion_input = _YieldHistogram((20.0,), edges=(1.0, 1.3), errors=(0.0,))
        baseline_template = _YieldHistogram((6.0,), edges=(1.0, 1.3), errors=(0.0,))
        baseline_final = _YieldHistogram((14.0,), edges=(1.0, 1.3), errors=(0.0,))
        dummy = _YieldHistogram((5.0,), edges=(1.0, 1.3), errors=(0.0,))
        payload = {
            "accepted": True,
            "H_pion_control_model": object(),
            "weights": object(),
            "H_pion_subtraction_template_MM": baseline_template,
            "H_pion_subtraction_template_MM_nosub": baseline_template.Clone("baseline_nosub"),
            "H_MM_before_pion_subtraction": pion_input,
            "H_MM_after_pion_subtraction": baseline_final,
        }
        entry = {
            "child_valid": True,
            "particle_subtraction_component_payload": payload,
            "H_MM_DUMMY_NORM": dummy,
            "bg_fit1_frac_err": 0.0,
            "bg_fit2_frac_err": 0.0,
        }
        cache = _cache()
        spec = {
            "label": "prompt", "cache_section": cache,
            "coefficient": source_coefficient,
            "base_coefficient": base_coefficient,
        }
        accepted = _accepted_rows(_row(coefficient=accepted_coefficient))
        patches = (
            mock.patch.object(calc, "OUTPATH", "unused"),
            mock.patch.object(calc, "get_active_bg_profile_name", return_value="no_empirical_residual"),
            mock.patch.object(calc, "resolve_particle_subtraction_mode", return_value="simc_shape_components"),
            mock.patch.object(calc, "resolve_pion_subtraction_scope", return_value="t_bin"),
            mock.patch.object(calc, "get_particle_subtraction_setting_key", return_value="Q4p4W2p74"),
            mock.patch.object(calc, "accepted_f6_3_artifact_paths", return_value={}),
            mock.patch.object(calc, "load_accepted_f6_3_authority", return_value=([{}], {}, {}, {"Left-lowe": "a", "f3": "b", "f4": "c"})),
            mock.patch.object(calc, "reconstruct_transient_factor_map", return_value=({("prompt", 7): factor}, {"accepted": True}, accepted)),
            mock.patch.object(calc, "iter_component_control_source_specs", return_value=[spec]),
            mock.patch.object(calc, "simc_shape_pion_weight_from_value", return_value=3.0),
            mock.patch.object(
                calc,
                "_component_cache_event_coefficient",
                side_effect=pion._component_cache_event_coefficient,
            ),
            mock.patch.object(calc, "clone_root_histogram", side_effect=self._clone),
        )
        def fill(templates, _specs, _reference, _weights, _context, _particle, _pol, *, method_a_event_multipliers=None):
            multiplier = method_a_event_multipliers[("prompt", 7)]
            templates["mm"].Fill(1.12, 6.0 * multiplier)
            templates["mm_nosub"].Fill(1.12, 6.0 * multiplier)
        patches += (mock.patch.object(calc, "fill_simc_shape_pion_subtraction_templates", side_effect=fill),)
        return {
            "calc": calc, "entry": entry, "cache": {"prompt": cache}, "patches": patches,
            "hist": {"phi_setting": "Left"},
            "processed": {"t_bin1phi_bin1": entry},
            "inp": {
                "ParticleType": "kaon", "EPSSET": "lowe", "POL": "unpolarized",
                "mm_min": 1.08, "mm_max": 1.18,
                "data_charge_err_left": 0.0, "dummy_charge_err_left": 0.0,
            },
            "measurements": {(0, 0): {"yield": 14.0, "statistical_error": 0.0, "total_error": 0.0}},
            "baseline_objects": (pion_input, baseline_template, baseline_final),
        }

    def _build(self, factor=1.0, *, fixture_kwargs=None, extra_patches=()):
        fixture = self._branch_fixture(factor, **(fixture_kwargs or {}))
        with ExitStack() as stack:
            for patcher in (*fixture["patches"], *extra_patches):
                stack.enter_context(patcher)
            source = fixture["calc"]._build_f6_3_parallel_method_a_source(
                fixture["hist"], fixture["processed"], fixture["cache"], [0.0, 1.0], [-180.0, 180.0],
                fixture["inp"], fixture["measurements"], 1.0, 1.0, 4,
            )
        return source, fixture

    def test_all_one_parallel_branch_reproduces_baseline_without_mutating_public_objects(self):
        source, fixture = self._build(1.0)
        self.assertTrue(source["available"])
        child = source["children"][0]
        self.assertEqual(child["B_pi_A"].contents, child["B_pi_0"].contents)
        self.assertEqual(child["MM_A"].contents, child["MM_0"].contents)
        self.assertEqual(child["YA"], child["Y0"])
        self.assertEqual(child["YA_statistical_error"], child["Y0_statistical_error"])
        self.assertEqual(child["YA_total_error"], child["Y0_total_error"])
        self.assertEqual([item.contents for item in fixture["baseline_objects"]], [[20.0], [6.0], [14.0]])
        self.assertEqual(source["lambda_integration_window"], [1.08, 1.18])
        self.assertNotIn("transient_factor_map", str(source))

    def test_nontrivial_parallel_branch_changes_only_the_private_pion_subtraction(self):
        source, fixture = self._build(1.5)
        child = source["children"][0]
        self.assertEqual(child["pion_input"].contents, [20.0])
        self.assertEqual(child["B_pi_0"].contents, [6.0])
        self.assertEqual(child["B_pi_A"].contents, [9.0])
        self.assertEqual(child["MM_0"].contents, [14.0])
        self.assertEqual(child["MM_A"].contents, [11.0])
        self.assertEqual(child["YA"], 11.0)
        self.assertEqual([item.contents for item in fixture["baseline_objects"]], [[20.0], [6.0], [14.0]])

    def test_effective_filler_coefficient_not_raw_cache_field_controls_live_parity(self):
        source, fixture = self._build(
            fixture_kwargs={
                "source_coefficient": 4.0,
                "base_coefficient": 2.0,
                "accepted_coefficient": 2.0,
            },
        )
        self.assertFalse(source["available"])
        self.assertIn("signed_source_coefficient", source["reason"])
        self.assertEqual(
            pion._component_cache_event_coefficient(
                {"cache_section": _cache(), "coefficient": 4.0, "base_coefficient": 2.0}, 0,
            ),
            4.0,
        )
        self.assertEqual([item.contents for item in fixture["baseline_objects"]], [[20.0], [6.0], [14.0]])

    def test_branch_runtime_failures_are_unavailable_without_partial_children(self):
        failure_patches = (
            lambda fixture: mock.patch.object(
                fixture["calc"], "clone_root_histogram", side_effect=RuntimeError("clone failure"),
            ),
            lambda fixture: mock.patch.object(
                fixture["calc"], "fill_simc_shape_pion_subtraction_templates", side_effect=RuntimeError("fill failure"),
            ),
            lambda _fixture: mock.patch.object(
                _YieldHistogram, "Add", side_effect=RuntimeError("add failure"),
            ),
        )
        for make_patch in failure_patches:
            with self.subTest(failure=make_patch.__code__.co_firstlineno):
                fixture = self._branch_fixture()
                original = [list(item.contents) for item in fixture["baseline_objects"]]
                with ExitStack() as stack:
                    for patcher in (*fixture["patches"], make_patch(fixture)):
                        stack.enter_context(patcher)
                    source = fixture["calc"]._build_f6_3_parallel_method_a_source(
                        fixture["hist"], fixture["processed"], fixture["cache"], [0.0, 1.0], [-180.0, 180.0],
                        fixture["inp"], fixture["measurements"], 1.0, 1.0, 4,
                    )
                self.assertFalse(source["available"])
                self.assertIn("f6_3_branch_exception:RuntimeError", source["reason"])
                self.assertNotIn("children", source)
                self.assertEqual([item.contents for item in fixture["baseline_objects"]], original)

    def test_unsupported_profile_or_mode_leaves_baseline_untouched_and_is_unavailable(self):
        fixture = self._branch_fixture()
        original = [list(item.contents) for item in fixture["baseline_objects"]]
        with mock.patch.object(fixture["calc"], "get_active_bg_profile_name", return_value="legacy"):
            source = fixture["calc"]._build_f6_3_parallel_method_a_source(
                fixture["hist"], fixture["processed"], fixture["cache"], [0.0, 1.0], [-180.0, 180.0],
                fixture["inp"], fixture["measurements"], 1.0, 1.0, 4,
            )
        self.assertFalse(source["available"])
        self.assertEqual([item.contents for item in fixture["baseline_objects"]], original)
        with mock.patch.object(fixture["calc"], "get_active_bg_profile_name", return_value="no_empirical_residual"), \
             mock.patch.object(fixture["calc"], "resolve_particle_subtraction_mode", return_value="legacy"):
            source = fixture["calc"]._build_f6_3_parallel_method_a_source(
                fixture["hist"], fixture["processed"], fixture["cache"], [0.0, 1.0], [-180.0, 180.0],
                fixture["inp"], fixture["measurements"], 1.0, 1.0, 4,
            )
        self.assertFalse(source["available"])

    def test_parallel_builder_uses_existing_cache_only_and_has_no_method_b_dependency(self):
        text = inspect.getsource(self.calculate_yield._build_f6_3_parallel_method_a_source).lower()
        self.assertNotIn("bin_data(", text)
        self.assertNotIn("pion_hgcer_method_b", text)
        self.assertNotIn("bg_fit(", text)

    @staticmethod
    def _snapshot_public_groups(groups):
        """Return a JSON-serializable complete snapshot of public yield groups."""
        child_keys = sorted((int(t_index), int(phi_index)) for t_index, phi_index in groups)
        return {
            "container_type": type(groups).__name__,
            "top_level_child_keys": [list(key) for key in child_keys],
            "children": {
                "{},{}".format(t_index, phi_index): {
                    "field_keys": sorted(str(field) for field in groups[(t_index, phi_index)]),
                    "values": {
                        str(field): float(value)
                        for field, value in groups[(t_index, phi_index)].items()
                    },
                }
                for t_index, phi_index in child_keys
            },
        }

    def _public_baseline_run(self, *, branch_runtime_failure=False):
        calc = self.calculate_yield
        fixture = _public_yield_fixture(calc)
        fixture["binned_dict"]["yield"] = fixture["binned_dict"].pop("kaon")
        component_payload = {
            "accepted": True,
            "fixture_schema": "f6_3_component_payload_sentinel/v1",
            "component_roles": ["pion", {"source": "prompt", "coefficient": -0.25}],
            "provenance": {"revision": 7, "labels": ["baseline", "preserve"]},
        }
        fixture["processed_entry"]["particle_subtraction_component_payload"] = component_payload
        component_payload_snapshot = copy.deepcopy(component_payload)
        component_payload_json_snapshot = json.dumps(
            component_payload, sort_keys=True, separators=(",", ":"),
        )
        integration = mock.Mock(side_effect=((10.0, 0.3), (5.0, 0.4)))
        common = (
            mock.patch.object(calc, "_resolve_bg_opt_prepass_data_cache", return_value=None),
            mock.patch.object(calc, "_get_cached_kaon_signal_shape_payload", return_value=None),
            mock.patch.object(calc, "_get_cached_kaon_sigma0_shape_payload", return_value=None),
            mock.patch.object(calc, "bin_data", return_value=(fixture["binned_dict"], None, None)),
            mock.patch.object(calc, "get_active_bg_profile_name", return_value="no_empirical_residual"),
            mock.patch.object(calc, "integral_with_stat_error", integration),
        )
        if branch_runtime_failure:
            branch = (
                mock.patch.object(calc, "resolve_particle_subtraction_mode", return_value="simc_shape_components"),
                mock.patch.object(calc, "resolve_pion_subtraction_scope", return_value="t_bin"),
                mock.patch.object(calc, "accepted_f6_3_artifact_paths", side_effect=RuntimeError("branch authority failure")),
            )
        else:
            branch = (mock.patch.object(calc, "resolve_particle_subtraction_mode", return_value="legacy"),)
        with ExitStack() as stack:
            for patcher in (*common, *branch):
                stack.enter_context(patcher)
            groups = calc.calculate_yield_data(
                "yield", fixture["hist"], [0.0, 1.0], [-180.0, 180.0], fixture["inp_dict"],
            )
        child = fixture["hist"]["_e8_2_baseline_stage_source"]["children"][0]
        return {
            "groups": self._snapshot_public_groups(groups),
            "stage_window_yields": dict(fixture["processed_entry"]["stage_window_yields"]),
            "scale_factor": fixture["binned_dict"]["yield"]["scale_factor"],
            "final_hist": list(fixture["final_hist"].contents),
            "final_errors": list(fixture["final_hist"].errors),
            "e8_2_yield": (
                child["final_yield"], child["statistical_error"], child["total_error"],
            ),
            "component_payload": fixture["processed_entry"]["particle_subtraction_component_payload"],
            "component_payload_original": component_payload,
            "component_payload_snapshot": component_payload_snapshot,
            "component_payload_json_snapshot": component_payload_json_snapshot,
            "f6_3": fixture["hist"]["_f6_3_parallel_method_a_source"],
        }

    def test_public_calculate_yield_baseline_outputs_survive_unavailable_and_runtime_failed_f6_3(self):
        unavailable = self._public_baseline_run()
        failed = self._public_baseline_run(branch_runtime_failure=True)
        self.assertEqual(unavailable["groups"], {
            "container_type": "defaultdict",
            "top_level_child_keys": [[0, 0]],
            "children": {
                "0,0": {
                    "field_keys": ["yield", "yield_err"],
                    "values": {"yield": 10.0, "yield_err": 0.3905124837953327},
                },
            },
        })
        self.assertEqual(unavailable["groups"], failed["groups"])
        self.assertEqual(unavailable["stage_window_yields"], {
            "raw_prompt": 21.0,
            "after_random_subtraction": 18.0,
            "after_dummy_subtraction": 16.0,
            "after_pion_subtraction": 10.0,
        })
        self.assertEqual(unavailable["stage_window_yields"], failed["stage_window_yields"])
        self.assertEqual(unavailable["scale_factor"], [[0.0]])
        self.assertEqual(unavailable["scale_factor"], failed["scale_factor"])
        self.assertEqual(unavailable["final_hist"], [10.0])
        self.assertEqual(unavailable["final_hist"], failed["final_hist"])
        self.assertEqual(unavailable["final_errors"], [0.3])
        self.assertEqual(unavailable["final_errors"], failed["final_errors"])
        self.assertEqual(unavailable["e8_2_yield"], (10.0, 0.3, 0.3905124837953327))
        self.assertEqual(unavailable["e8_2_yield"], failed["e8_2_yield"])
        for result in (unavailable, failed):
            self.assertIs(result["component_payload"], result["component_payload_original"])
            self.assertEqual(result["component_payload"], result["component_payload_snapshot"])
            self.assertEqual(
                json.dumps(result["component_payload"], sort_keys=True, separators=(",", ":")),
                result["component_payload_json_snapshot"],
            )
        self.assertFalse(unavailable["f6_3"]["available"])
        self.assertFalse(failed["f6_3"]["available"])
        self.assertIn("f6_3_branch_exception:RuntimeError", failed["f6_3"]["reason"])


class BranchBoundaryTests(unittest.TestCase):
    def test_unavailable_source_preserves_baseline_and_all_forbidden_boundaries(self):
        source = f63.unavailable_parallel_source("authority_missing", setting_id="Left-lowe")
        self.assertFalse(source["available"])
        self.assertTrue(source["baseline_public_output_unchanged"])
        self.assertFalse(source["production_promotion_performed"])
        self.assertFalse(source["method_b_numerical_dependency"])
        self.assertFalse(source["empirical_residual_used"])
        self.assertFalse(source["event_correction_persisted"])

    def test_helper_is_pure_python_and_does_not_traverse_trees(self):
        text = (ROOT / "src" / "cuts" / "pion_hgcer_method_a_parallel_full_procedure.py").read_text(encoding="utf-8")
        self.assertNotIn("import ROOT", text)
        self.assertNotIn("GetEntries", text)


if __name__ == "__main__":
    unittest.main()
