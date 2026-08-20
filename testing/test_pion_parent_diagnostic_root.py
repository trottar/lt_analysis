"""PyROOT-gated proposal/final parent diagnostic behavior."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

try:
    import ROOT
except ImportError:  # pragma: no cover - reported as a skipped gated suite.
    ROOT = None


@unittest.skipUnless(ROOT is not None, "PyROOT is not available")
class PionParentDiagnosticRootTests(unittest.TestCase):
    @staticmethod
    def _hist(name, values):
        hist = ROOT.TH1D(name, name, len(values), 1.0, 1.3)
        hist.SetDirectory(0)
        for index, value in enumerate(values, start=1):
            hist.SetBinContent(index, float(value))
            hist.SetBinError(index, float(value) ** 0.5)
        return hist

    def setUp(self):
        try:
            import rand_sub
        except ImportError as exc:  # Production dependencies are required for this integration path.
            self.skipTest("rand_sub dependencies are unavailable: {}".format(exc))
        self.rand_sub = rand_sub

    def _proposal_payload(self):
        before = self._hist("parent_before", (8.0, 10.0, 12.0))
        full_before = self._hist("parent_full_before", (10.0, 12.0, 14.0))
        proposed_after = self._hist("parent_proposed_after", (9.0, 10.0, 11.0))
        model = self._hist("parent_model", (1.0, 2.0, 3.0))
        return {
            "analysis_scope": "t_bin1",
            "input_selection": "no_rf_proton_cleaning_then_rf_restored",
            "source_target_state": "post_proton_post_rf",
            "proposal_available": True,
            "H_MM_before_pion_subtraction": before,
            "H_MM_nosub_before_pion_subtraction": full_before,
            "H_MM_nosub_after_pion_subtraction": proposed_after,
            "H_pion_weight_vs_MM": model,
            "H_pion_control_model": model,
            "H_kaon_pion_model": model,
        }

    def test_rejected_evaluable_proposal_keeps_proposal_and_draws_zero_final(self):
        proposal = self._proposal_payload()
        production = {"accepted": False, "fallback_mode": "zero", "reason": "quality gate"}
        final, status = self.rand_sub.resolve_parent_diagnostic_final_application(
            proposal, production
        )
        self.assertEqual(status["final_status"], "zero")
        self.assertIsNotNone(final)
        self.assertNotEqual(
            proposal["H_pion_control_model"].Integral(),
            final["H_pion_subtraction_template_MM"].Integral(),
        )
        self.assertAlmostEqual(
            final["H_MM_nosub_before_pion_subtraction"].Integral(),
            final["H_MM_nosub_after_pion_subtraction"].Integral(),
        )

    def test_accepted_and_nonapplying_final_policy_states_are_distinct(self):
        proposal = self._proposal_payload()
        final, status = self.rand_sub.resolve_parent_diagnostic_final_application(
            proposal, {"accepted": True, "fallback_mode": "error", "reason": ""}
        )
        self.assertEqual(status["final_status"], "applied_component")
        self.assertIs(
            final["H_MM_nosub_after_pion_subtraction"],
            proposal["H_MM_nosub_after_pion_subtraction"],
        )

        final, status = self.rand_sub.resolve_parent_diagnostic_final_application(
            proposal, {"accepted": False, "fallback_mode": "skip_bin", "reason": "rejected"}
        )
        self.assertIsNone(final)
        self.assertEqual(status["final_status"], "skip_bin")

    def test_single_scale_state_uses_only_the_explicit_fallback_factory(self):
        proposal = self._proposal_payload()
        expected = {"accepted": True, "final_application_status": "applied_fallback"}
        factory = mock.Mock(return_value=expected)
        final, status = self.rand_sub.resolve_parent_diagnostic_final_application(
            proposal,
            {"accepted": False, "fallback_mode": "single_scale", "reason": "rejected"},
            fallback_context=factory,
        )
        factory.assert_called_once_with()
        self.assertIs(final, expected)
        self.assertEqual(status["final_status"], "applied_fallback")


if __name__ == "__main__":
    unittest.main()
