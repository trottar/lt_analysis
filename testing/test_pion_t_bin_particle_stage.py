"""Static integration contracts for particle-stage t-bin pion parents."""

from __future__ import annotations

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
RAND_SUB = REPO_ROOT / "src" / "cuts" / "rand_sub.py"
PARENT_STAGE = REPO_ROOT / "src" / "cuts" / "pion_t_bin_parents.py"
MAIN = REPO_ROOT / "src" / "main.py"
AVE = REPO_ROOT / "src" / "binning" / "ave_per_bin.py"


class ParticleStageTBinPionParentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rand_source = RAND_SUB.read_text(encoding="utf-8")
        cls.parent_source = PARENT_STAGE.read_text(encoding="utf-8")
        cls.main_source = MAIN.read_text(encoding="utf-8")
        cls.ave_source = AVE.read_text(encoding="utf-8")

    def test_rand_sub_builds_only_after_committed_proton_cleaning(self):
        builder_call = self.rand_source.index("build_setting_t_bin_pion_parents(\n                    histDict")
        committed_targets = self.rand_source.index(
            "active_component_targets = proton_cleaning_application.get(\"final_targets\")"
        )
        self.assertGreater(builder_call, committed_targets)
        self.assertIn("parent_inputs=parent_inputs", self.rand_source[builder_call:])
        self.assertIn("_build_authoritative_pion_control_source_cache", self.rand_source)
        self.assertIn("canonical_t_products", self.rand_source)
        self.assertIn("H_proton_cleaned_final_rf", self.rand_source)
        ave_source = (REPO_ROOT / "src" / "binning" / "ave_per_bin.py").read_text(encoding="utf-8")
        self.assertNotIn("t_integrated_fit_only", self.parent_source)
        self.assertNotIn("process_hist_data", self.parent_source)
        self.assertNotIn("prune_hist", self.parent_source)
        self.assertIn("fingerprint_histogram_content_error", self.parent_source)
        self.assertIn("retired_t_integrated_parent_comparison_route", ave_source)
        self.assertIn("if hgcer_cutg is None and ParticleType", ave_source)
        self.assertIn("final_cleaned_factor", ave_source)
        self.assertIn("strict=True", ave_source)
        for source_label in ("prompt", "rand", "dummy_prompt", "dummy_rand"):
            self.assertIn('source_label="{}"'.format(source_label), ave_source)

    def test_component_children_consume_the_particle_stage_control_cache(self):
        yield_source = (REPO_ROOT / "src" / "binning" / "calculate_yield.py").read_text(encoding="utf-8")
        self.assertIn("_authoritative_pion_control_cache_by_phi", yield_source)
        self.assertIn("child_event_cache", yield_source)
        self.assertIn("missing_authoritative_pion_control_cache_for_t_bin_consumer", yield_source)
        self.assertIn("one_signed_pion_control_source_cache_no_proton_weight", self.rand_source)
        for field in ("source_label", "entry_index", "coefficient", "canonical_t_counts"):
            self.assertIn('"{}"'.format(field), self.rand_source)

    def test_parent_and_setting_wide_diagnostics_apply_the_cache_not_trees(self):
        callback_start = self.rand_source.index("def _build_parent_diagnostic_application(")
        callback_end = self.rand_source.index("build_setting_t_bin_pion_parents(", callback_start)
        callback_source = self.rand_source[callback_start:callback_end]
        self.assertIn("_build_authoritative_parent_mm_diagnostic_proposal", callback_source)
        self.assertIn("_build_authoritative_parent_single_scale_final", callback_source)
        self.assertNotIn("build_component_pion_application_proposal", callback_source)
        self.assertNotIn("sub_tree_bundle", callback_source)
        self.assertIn("H_pion_control_global", self.rand_source)
        self.assertIn("if not authoritative_component_t_bin", self.rand_source)

    def test_parent_state_and_page_contract_are_explicit(self):
        for field in (
            '"fit_result": fit_result',
            '"proposed_diagnostic_application_result"',
            '"proposed_diagnostic_application_status"',
            '"final_diagnostic_application_result"',
            '"final_diagnostic_application_status"',
            '"production_application_policy"',
            '"diagnostic_application_result"',
            '"diagnostic_application_status"',
            '"application_result": None',
            '"input_selection": "no_rf_proton_cleaning_then_rf_restored"',
            '"source_target_state": "post_proton_post_rf"',
        ):
            self.assertIn(field, self.parent_source)
        for page_key in (
            "pion_control_fit",
            "protected_fit_or_status",
            "weight",
            "model_closure",
            "event_template_closure",
            "before_after",
            "lambda_comparison",
            "overview",
            "application_gate",
            "proposal_final_transition",
        ):
            self.assertIn(page_key, self.parent_source)

    def test_parent_diagnostic_preserves_proposal_but_applies_final_policy(self):
        self.assertIn("evaluate_component_pion_application_proposal", self.rand_source)
        self.assertIn("proposal_only=True", self.rand_source)
        self.assertIn("resolve_parent_diagnostic_final_application", self.rand_source)
        self.assertIn("_build_zero_parent_diagnostic_final", self.rand_source)
        self.assertIn("_build_single_scale_parent_diagnostic_final", self.rand_source)
        self.assertIn("pion_parent_diagnostic_strict", self.parent_source)
        self.assertIn("_parent_plot_contract", self.parent_source)

    def test_step_six_and_averages_are_consumer_only(self):
        self.assertNotIn("build_t_bin_pion_parents", self.main_source)
        self.assertIn("validate_frozen_t_bin_pion_parent_collection", self.main_source)
        self.assertNotIn("def build_t_bin_pion_parents", self.ave_source)
        self.assertNotIn("parent_only", self.ave_source)
        self.assertNotIn('hist["_pion_t_bin_parent_results"] =', self.ave_source)
        # The old averages template builder remains for non-component modes,
        # but canonical-t production must bypass both its construction and
        # its scalar-fallback application.
        self.assertGreaterEqual(
            self.ave_source.count('if ParticleType == "kaon" and not use_t_bin_parents'),
            3,
        )
        tbin_branch = self.ave_source[
            self.ave_source.index("if use_t_bin_parents and not t_integrated_fit_only:"):
            self.ave_source.index("elif component_shape_payload is not None:")
        ]
        self.assertIn("use_legacy_scalar_subtraction = False", tbin_branch)
        self.assertNotIn("particle_subtraction_ave(", tbin_branch)

    def test_particle_stage_is_the_single_normal_production_assignment(self):
        assignments = []
        pattern = re.compile(r'\["_pion_t_bin_parent_results"\]\s*=')
        for source_path in (REPO_ROOT / "src").rglob("*.py"):
            matches = list(pattern.finditer(source_path.read_text(encoding="utf-8")))
            assignments.extend((source_path.relative_to(REPO_ROOT), match.start()) for match in matches)
        self.assertEqual(
            [path.as_posix() for path, _ in assignments],
            ["src/cuts/pion_t_bin_parents.py"],
        )


if __name__ == "__main__":
    unittest.main()
