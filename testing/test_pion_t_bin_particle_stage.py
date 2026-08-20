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
        self.assertIn("proton_cleaning_result=proton_cleaning_result", self.rand_source[builder_call:])
        self.assertIn("t_integrated_fit_only=True", self.parent_source)
        ave_source = (REPO_ROOT / "src" / "binning" / "ave_per_bin.py").read_text(encoding="utf-8")
        self.assertIn("final_cleaned_factor", ave_source)
        self.assertIn("strict=True", ave_source)
        for source_label in ("prompt", "rand", "dummy_prompt", "dummy_rand"):
            self.assertIn('source_label="{}"'.format(source_label), ave_source)

    def test_parent_state_and_page_contract_are_explicit(self):
        for field in (
            '"fit_result": fit_result',
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
        ):
            self.assertIn('"{}"'.format(page_key), self.parent_source)

    def test_step_six_and_averages_are_consumer_only(self):
        self.assertNotIn("build_t_bin_pion_parents", self.main_source)
        self.assertIn("validate_frozen_t_bin_pion_parent_collection", self.main_source)
        self.assertNotIn("def build_t_bin_pion_parents", self.ave_source)
        self.assertNotIn("parent_only", self.ave_source)
        self.assertNotIn('hist["_pion_t_bin_parent_results"] =', self.ave_source)

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
