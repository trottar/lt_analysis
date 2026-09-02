"""Static contract checks for the detached Phase-D rand_sub integration."""

from __future__ import annotations

from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
RAND_SUB = REPO_ROOT / "src" / "cuts" / "rand_sub.py"


class PionHGCerPhaseDRuntimeContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = RAND_SUB.read_text(encoding="utf-8")

    def test_phase_d_imports_and_runtime_order_are_detached(self):
        source = self.source
        for name in (
            "build_pion_hgcer_comparison_input_contract",
            "build_pion_hgcer_method_a_comparison",
            "build_pion_hgcer_method_b_comparison",
            "build_pion_hgcer_ab_comparison",
            "build_pion_hgcer_phase_d_checkpoint",
            "write_pion_hgcer_phase_d_checkpoint_json",
            "render_pion_hgcer_ab_comparison_pages",
        ):
            self.assertIn(name, source)
        phase_c_write = source.index("write_pion_hgcer_refinement_checkpoint_json(\n")
        d1 = source.index("build_pion_hgcer_comparison_input_contract(\n", phase_c_write)
        d2 = source.index("build_pion_hgcer_method_a_comparison(\n", d1)
        d3 = source.index("build_pion_hgcer_method_b_comparison(\n", d1)
        d4 = source.index("build_pion_hgcer_ab_comparison(\n", d1)
        phase_d_build = source.index("build_pion_hgcer_phase_d_checkpoint(\n", d4)
        phase_d_write = source.index("write_pion_hgcer_phase_d_checkpoint_json(\n", phase_d_build)
        phase_d_render = source.index("render_pion_hgcer_ab_comparison_pages(\n", phase_d_write)
        self.assertLess(phase_c_write, d1)
        self.assertLess(d1, d2)
        self.assertLess(d1, d3)
        self.assertLess(max(d2, d3), d4)
        self.assertLess(d4, phase_d_build)
        self.assertLess(phase_d_build, phase_d_write)
        self.assertLess(phase_d_write, phase_d_render)

    def test_phase_d_block_has_no_application_or_correction_identifiers(self):
        source = self.source
        start = source.index("# Phase D.4 is a terminal detached diagnostic chain.")
        end = source.index(
            "                except Exception as exc:\n                    # A checkpoint is a persistent diagnostic artifact only;",
            start,
        )
        block = source[start:end]
        for forbidden in (
            "C_A",
            "C_B",
            "C_final",
            "refined_pion_weight",
            "applied_refinement_weight",
            "w_pi_refined",
        ):
            self.assertNotIn(forbidden, block)
        self.assertIn("pion_hgcer_method_a_comparison", block)
        self.assertIn("pion_hgcer_method_b_comparison", block)
        self.assertIn("pion_hgcer_ab_comparison", block)


if __name__ == "__main__":
    unittest.main()
