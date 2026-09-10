"""Static contracts for detached Phase-E runtime wiring."""

from __future__ import annotations

from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
RAND_SUB = REPO_ROOT / "src" / "cuts" / "rand_sub.py"
YIELD_SOURCE = REPO_ROOT / "src" / "binning" / "calculate_yield.py"
PION_COMPONENT_SOURCE = REPO_ROOT / "src" / "cuts" / "pion_component_subtraction.py"


class PionHGCerPhaseERuntimeContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = RAND_SUB.read_text(encoding="utf-8")

    def test_e2_is_built_after_phase_d_and_d11_then_rendered_keyword_only(self):
        source = self.source
        self.assertIn("build_full_background_subtraction_e2_payload", source)
        start = source.index(
            "# Phases D.6 through D.11 and E.2 through E.6 are detached terminal"
        )
        phase_d_write = source.index("write_pion_hgcer_phase_d_checkpoint_json(\n")
        d11_build = source.index(
            "build_full_background_subtraction_d11_payload(\n", start
        )
        e2_build = source.index(
            "build_full_background_subtraction_e2_payload(\n", start
        )
        render = source.index(
            "render_full_background_subtraction_procedure_pages(\n", start
        )
        self.assertLess(phase_d_write, d11_build)
        self.assertLess(d11_build, e2_build)
        self.assertLess(e2_build, render)

        e2_call = source[e2_build:source.index(
            "            full_background_subtraction_e3_payload", e2_build
        )]
        self.assertIn("pion_hgcer_tdelta_diagnostic", e2_call)
        self.assertIn("pion_hgcer_event_contract", e2_call)
        for forbidden in (
            "pion_hgcer_method_a",
            "pion_hgcer_method_b",
            "pion_hgcer_ab_comparison",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, e2_call)

        render_call = source[render:source.index("                    )\n", render) + 22]
        self.assertIn(
            "e2_payload=full_background_subtraction_e2_payload", render_call
        )

    def test_e7_1_follows_e7_and_uses_only_the_frozen_parent_inputs(self):
        source = self.source
        start = source.index(
            "# Phases D.6 through D.11 and E.2 through E.6 are detached terminal"
        )
        e6_build = source.index("build_full_background_subtraction_e6_payload(\n", start)
        prototype_build = source.index(
            "build_pion_hgcer_ab_combination_prototype(\n", start
        )
        correction_build = source.index(
            "build_pion_hgcer_parent_preserving_correction(\n", start
        )
        artifact_build = source.index(
            "build_pion_hgcer_parent_preserving_correction_artifact(\n", start
        )
        artifact_write = source.index(
            "write_pion_hgcer_parent_preserving_correction_json(\n", start
        )
        e7_presentation = source.index(
            "build_full_background_subtraction_e7_payload(\n", start
        )
        render = source.index(
            "render_full_background_subtraction_procedure_pages(\n", start
        )
        self.assertLess(e6_build, prototype_build)
        self.assertLess(prototype_build, correction_build)
        self.assertLess(correction_build, artifact_build)
        self.assertLess(artifact_build, artifact_write)
        self.assertLess(artifact_write, e7_presentation)
        self.assertLess(e7_presentation, render)

        correction_call = source[
            correction_build:source.index(
                '            histDict["_pion_hgcer_parent_preserving_correction"]',
                correction_build,
            )
        ]
        self.assertIn("pion_hgcer_ab_combination_prototype", correction_call)
        self.assertIn("pion_hgcer_event_contract", correction_call)
        for forbidden in (
            "histDict[", "pion_hgcer_method_a", "pion_hgcer_method_b",
            "pion_hgcer_ab_comparison", "phase_d_checkpoint_for_plots",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, correction_call)

        e7_call = source[e7_presentation:source.index("            if (", e7_presentation)]
        self.assertIn("pion_hgcer_ab_combination_prototype", e7_call)
        self.assertNotIn("pion_hgcer_parent_preserving_correction", e7_call)
        self.assertIn('"_pion_hgcer_parent_preserving_correction"', source)
        self.assertIn('"pion_hgcer_parent_preserving_correction_summary"', source)
        self.assertIn('"pion_hgcer_parent_preserving_correction_artifacts"', source)
        self.assertIn("parent_preserving_correction_artifact_write_exception", source)

    def test_e7_1_has_no_production_consumer_or_pion_component_input(self):
        correction_identifier = "pion_hgcer_parent_preserving_correction"
        self.assertNotIn(correction_identifier, YIELD_SOURCE.read_text(encoding="utf-8"))
        self.assertNotIn(
            correction_identifier,
            PION_COMPONENT_SOURCE.read_text(encoding="utf-8"),
        )


if __name__ == "__main__":
    unittest.main()
