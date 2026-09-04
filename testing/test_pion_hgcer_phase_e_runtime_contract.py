"""Static contract checks for detached Phase-E E.2 runtime wiring."""

from __future__ import annotations

from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
RAND_SUB = REPO_ROOT / "src" / "cuts" / "rand_sub.py"


class PionHGCerPhaseERuntimeContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = RAND_SUB.read_text(encoding="utf-8")

    def test_e2_is_built_after_phase_d_and_d11_then_rendered_keyword_only(self):
        source = self.source
        self.assertIn("build_full_background_subtraction_e2_payload", source)
        start = source.index(
            "# Phases D.6 through D.11 and E.2 are terminal presentation only."
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

        e2_call = source[e2_build:source.index("            if (", e2_build)]
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


if __name__ == "__main__":
    unittest.main()
