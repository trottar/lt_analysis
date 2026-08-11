"""Focused non-ROOT coverage for the timing-|t| Lambda-preservation gate."""

from __future__ import annotations

import sys
import types
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "utility"))
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))


try:
    import ROOT  # noqa: F401
    import proton_contamination_weights as proton_cleaning
except ImportError:
    class _RootImportStub:
        def __init__(self, *_args, **_kwargs):
            pass

    root_stub = types.ModuleType("ROOT")
    for name in (
        "TCanvas", "TLatex", "TLegend", "TLine", "TPad", "TPaveText", "TH1D", "TH2D",
        "TF1", "TGraphErrors",
    ):
        setattr(root_stub, name, _RootImportStub)
    root_stub.gPad = types.SimpleNamespace()
    root_stub.gStyle = types.SimpleNamespace()
    for name in ("kBlack", "kBlue", "kGray", "kGreen", "kMagenta", "kOrange", "kRed", "kViolet"):
        setattr(root_stub, name, 1)
    with mock.patch.dict(sys.modules, {"ROOT": root_stub}):
        sys.modules.pop("proton_contamination_weights", None)
        import proton_contamination_weights as proton_cleaning


def _config(**gate_overrides):
    gate = {
        "enabled": True,
        "validation_window_key": "lambda_peak",
        "maximum_lambda_removed_fraction": 0.10,
        "minimum_raw_prompt_events": 0,
        "minimum_absolute_support": None,
        "minimum_positive_signed_yield": None,
        "insufficient_support_policy": "bypass",
        "affects_event_weights": True,
        "affects_fit_acceptance": False,
    }
    gate.update(gate_overrides)
    return {
        "validation_windows": {"lambda_peak": (1.105, 1.125)},
        "lambda_preservation_gate": gate,
    }


def _result():
    return {
        "accepted": True,
        "diagnostics": {"setting_support": {"support_label": "supported"}},
    }


def _row(weight, probability, *, prompt=True, entry=0, rf_accept=True, mm=1.115):
    return {
        "source_label": "prompt" if prompt else "random",
        "source_entry_index": entry,
        "is_prompt_source": prompt,
        "physical_weight": weight,
        "proposed_proton_probability": probability,
        "proposed_cleaned_factor": 1.0 - probability,
        "proposed_final_cleaned_factor": (1.0 - probability) if rf_accept else 0.0,
        "rf_accept": rf_accept,
        "adj_mm": mm,
        "adj_t": 0.5,
        "selected_timing": 0.1,
    }


class LambdaPreservationGateTests(unittest.TestCase):
    def test_resolves_only_the_shared_validation_window(self):
        config = _config()
        resolved = proton_cleaning._resolve_lambda_preservation_gate_config(config)
        self.assertEqual(resolved["validation_window_key"], "lambda_peak")
        self.assertEqual((resolved["window_low"], resolved["window_high"]), (1.105, 1.125))
        self.assertNotIn("lambda_window", config["lambda_preservation_gate"])
        with self.assertRaisesRegex(ValueError, "validation window key"):
            proton_cleaning._resolve_lambda_preservation_gate_config(
                _config(validation_window_key="not_a_window")
            )

    def test_inclusive_threshold_and_negative_fraction_policy(self):
        at_limit = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(10.0, 0.10)], _config()
        )
        self.assertEqual(at_limit["status"], "pass")
        self.assertAlmostEqual(at_limit["proposed_removed_fraction"], 0.10)
        above_limit = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(10.0, 0.1000001)], _config()
        )
        self.assertEqual(above_limit["status"], "fail")
        negative_fraction = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(10.0, 0.0, prompt=True), _row(-9.0, 1.0, prompt=False, entry=1)],
            _config(),
        )
        self.assertEqual(negative_fraction["status"], "pass")
        self.assertLess(negative_fraction["proposed_removed_fraction"], 0.0)
        self.assertIn("negative_proposed_lambda_yield", negative_fraction["observational_warnings"])
        self.assertIn("negative_lambda_removed_fraction", negative_fraction["observational_warnings"])

    def test_support_requirements_and_cancellation_ratio(self):
        insufficient_prompt = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(4.0, 0.02)], _config(minimum_raw_prompt_events=2)
        )
        self.assertEqual(insufficient_prompt["status"], "insufficient_support")
        self.assertIn("minimum_raw_prompt_events_not_met", insufficient_prompt["support_reasons"])
        zero_signed = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(1.0, 0.1), _row(-1.0, 0.1, prompt=False, entry=1)], _config()
        )
        self.assertEqual(zero_signed["status"], "insufficient_support")
        self.assertIn("raw_signed_yield_not_positive_or_nonfinite", zero_signed["support_reasons"])
        valid = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(3.0, 0.05), _row(-1.0, 0.05, prompt=False, entry=1)], _config()
        )
        self.assertEqual(valid["status"], "pass")
        self.assertAlmostEqual(valid["raw_signed_to_absolute_support_ratio"], 0.5)

    def test_optional_support_thresholds_and_invalid_policy(self):
        below_absolute = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(2.0, 0.01)], _config(minimum_absolute_support=3.0)
        )
        self.assertEqual(below_absolute["status"], "insufficient_support")
        self.assertIn("minimum_absolute_support_not_met", below_absolute["support_reasons"])
        below_signed = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(2.0, 0.01)], _config(minimum_positive_signed_yield=3.0)
        )
        self.assertEqual(below_signed["status"], "insufficient_support")
        self.assertIn("minimum_positive_signed_yield_not_met", below_signed["support_reasons"])
        with self.assertRaisesRegex(ValueError, "insufficient-support policy"):
            proton_cleaning._resolve_lambda_preservation_gate_config(
                _config(insufficient_support_policy="apply")
            )

    def test_gate_uses_pre_rf_proposed_probability(self):
        result = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(10.0, 0.20, rf_accept=False)], _config()
        )
        self.assertEqual(result["status"], "fail")
        self.assertAlmostEqual(result["proposed_removed_fraction"], 0.20)

    def test_disabled_gate_preserves_normal_commit_path(self):
        result = proton_cleaning._evaluate_lambda_preservation_gate(
            _result(), [_row(10.0, 0.20)], _config(enabled=False)
        )
        self.assertEqual(result["status"], "disabled")
        self.assertEqual(result["production_action"], "apply")
        self.assertTrue(result["proton_cleaning_committed"])

    def test_commit_bypass_preserves_proposal_and_rf_membership(self):
        first = _row(1.0, 0.35, entry=0, rf_accept=False)
        second = _row(1.0, 0.25, entry=1, rf_accept=True)
        lookup = {
            "prompt:0": dict(first),
            "prompt:1": dict(second),
        }
        gate = {"status": "fail", "production_action": "bypass"}
        proton_cleaning._commit_lambda_preservation_gate(lookup, [first, second], gate)
        self.assertEqual(lookup["prompt:0"]["proposed_proton_probability"], 0.35)
        self.assertEqual(lookup["prompt:0"]["applied_proton_probability"], 0.0)
        self.assertEqual(lookup["prompt:0"]["applied_cleaned_factor"], 1.0)
        self.assertFalse(lookup["prompt:0"]["rf_accept"])
        self.assertEqual(lookup["prompt:0"]["applied_final_cleaned_factor"], 0.0)
        self.assertEqual(lookup["prompt:1"]["final_cleaned_factor"], 1.0)
        self.assertEqual(first["proton_weight"], 0.0)
        self.assertEqual(first["cleaned_factor"], 1.0)
        self.assertEqual(gate["final_applied_closure"]["pre_rf_closure_difference"], 0.0)
        self.assertTrue(gate["proposed_model_closure"]["cell_delta_closure_preserved"])

    def test_commit_pass_copies_the_proposed_lookup(self):
        row = _row(2.0, 0.35, entry=3)
        lookup = {"prompt:3": dict(row)}
        gate = {"status": "pass", "production_action": "apply"}
        proton_cleaning._commit_lambda_preservation_gate(lookup, [row], gate)
        payload = lookup["prompt:3"]
        self.assertEqual(payload["applied_proton_probability"], payload["proposed_proton_probability"])
        self.assertEqual(payload["cleaned_factor"], payload["proposed_cleaned_factor"])
        self.assertEqual(gate["final_applied_closure"]["pre_rf_closure_difference"], 0.0)

    def test_audit_is_bounded_and_deterministic(self):
        rows = [_row(1.0, 0.1, entry=index) for index in range(8)]
        first = proton_cleaning._lambda_gate_audit_rows(rows, limit=3)
        second = proton_cleaning._lambda_gate_audit_rows(list(reversed(rows)), limit=3)
        self.assertEqual(first["sampled_entry_count"], 3)
        self.assertEqual(first["rows"], second["rows"])
        self.assertIn("applied_zero_reason", first["rows"][0])


if __name__ == "__main__":
    unittest.main()
