"""Static runtime contracts for detached Phase-F.1 wiring."""

from __future__ import annotations

from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
RAND_SUB = REPO_ROOT / "src" / "cuts" / "rand_sub.py"
YIELD_SOURCE = REPO_ROOT / "src" / "binning" / "calculate_yield.py"
PION_COMPONENT_SOURCE = REPO_ROOT / "src" / "cuts" / "pion_component_subtraction.py"


class PionHGCerPhaseFRuntimeContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = RAND_SUB.read_text(encoding="utf-8")

    def test_f1_is_built_after_method_a_from_frozen_part1_phase_a_cache_and_phi_edges(self):
        source = self.source
        phase_a = source.index("pion_hgcer_event_contract = build_pion_hgcer_event_contract(\n")
        f1 = source.index("pion_hgcer_method_a_acceptance_contract = (\n")
        method_a = source.index("pion_hgcer_method_a = build_pion_hgcer_method_a(\n")
        self.assertLess(phase_a, method_a)
        self.assertLess(method_a, f1)
        call = source[f1:source.index('                histDict["pion_hgcer_method_a_acceptance_contract_summary"]', f1)]
        self.assertIn("pion_hgcer_tdelta_diagnostic", call)
        self.assertIn("pion_hgcer_method_a", call)
        self.assertIn("pion_hgcer_event_contract", call)
        self.assertIn("pion_control_cache", call)
        self.assertIn("phi_edges=frozen_phi_bins", call)
        self.assertNotIn('histDict["_pion_hgcer_method_a_acceptance_contract"]', source)
        for forbidden in (
            "pion_hgcer_method_b", "pion_hgcer_ab_comparison",
            "phase_d_checkpoint_for_plots", "pion_hgcer_ab_combination_prototype",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, call)

    def test_f1_artifact_and_keyword_only_renderer_wiring_are_detached(self):
        source = self.source
        f1 = source.index("pion_hgcer_method_a_acceptance_contract = (\n")
        artifact = source.index("build_pion_hgcer_method_a_acceptance_contract_artifact(\n", f1)
        artifact_write = source.index("write_pion_hgcer_method_a_acceptance_contract_json(\n", artifact)
        render = source.index("render_full_background_subtraction_procedure_pages(\n", f1)
        self.assertLess(f1, artifact)
        self.assertLess(artifact, artifact_write)
        self.assertLess(artifact_write, render)
        render_call = source[render:source.index("                    )\n", render) + 22]
        self.assertIn("f1_payload=full_background_subtraction_f1_payload", render_call)
        self.assertNotIn('histDict["full_background_subtraction_f1_payload"]', source)
        self.assertNotIn('histDict["full_background_subtraction_f1"]', source)

    def test_cache_preserves_only_scalar_f1_provenance_and_production_has_no_consumer(self):
        source = self.source
        cache_start = source.index("def _build_authoritative_pion_control_source_cache(")
        cache_end = source.index("def _fill_mm_only_authoritative_pion_control_templates(", cache_start)
        cache = source[cache_start:cache_end]
        for field in (
            '"ssxptar"', '"ssyptar"', '"hsxptar"', '"hsyptar"',
            '"phi_degrees"', '"phi_index"', '"P_hgcer_xAtCer"',
            '"P_hgcer_yAtCer"', '"ssdelta"',
        ):
            with self.subTest(field=field):
                self.assertIn(field, cache)
        identifier = "pion_hgcer_method_a_acceptance_contract"
        self.assertNotIn(identifier, YIELD_SOURCE.read_text(encoding="utf-8"))
        self.assertNotIn(identifier, PION_COMPONENT_SOURCE.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
