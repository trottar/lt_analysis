"""Deterministic audit tests; no ROOT or production event inputs."""
import copy
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "cuts"))
import pion_control_anchor_audit as audit


class Histogram:
    def __init__(self, contents=(-2.0, 0.0), errors=(math.sqrt(30), 0.0)):
        self.contents = list(contents)
        self.errors = list(errors)
        self.edges = [0.88, 0.91, 0.94]

    def GetNbinsX(self):
        return len(self.contents)

    def GetXaxis(self):
        return self

    def GetBinLowEdge(self, i):
        return self.edges[i - 1]

    def GetBinContent(self, i):
        if i == 0 or i == self.GetNbinsX() + 1:
            return 0.
        return self.contents[i - 1]

    def GetBinError(self, i):
        if i == 0 or i == self.GetNbinsX() + 1:
            return 0.
        return self.errors[i - 1]


def inputs(y=-2.0):
    records = [dict(source_label=s, coefficient=c, adj_MM=0.88, nommcuts=True)
               for s, c in zip(audit.SOURCES, (1., -2., -4., 3.))]
    d = dict(fit_bin_indices=[1], fit_y=[y], fit_sigma=[math.sqrt(30)],
             fit_template=[1.], weighted_numerator=y / 30,
             weighted_denominator=1 / 30, estimator="least_squares")
    step = dict(pass_index=1, component_name="pi_n", amplitude=max(0., y),
                sigma=math.sqrt(30), fit_diagnostics=d, anchor_windows=[[.88, .94]])
    solve = dict(success=True, message="", amplitude=max(0., y))
    parent = dict(t_bin_index=1, t_edges=[.567, .733], pion_parent_id="parent-t2",
                  coordinate_fingerprint="coord", H_pion_control=Histogram(),
                  fit_result=dict(H_pion_fit_step_overlays=[step], B_n=0.0002021,
                      diagnostics=dict(pion=dict(staged_pass_history=[dict(pi_n=solve)]))))
    hist = dict(phi_setting="Center", _authoritative_pion_control_source_cache=dict(
        by_t=[dict(t_index=1, records=records)]))
    inp = dict(ParticleType="kaon", Q2="4p4", W="2p74", EPSSET="low")
    return hist, inp, [parent], {}


class AnchorAuditTests(unittest.TestCase):
    def result(self, data):
        return audit.build_anchor_audit(*data)["parents"][0]

    def test_negative_signed_numerator_and_variance(self):
        result = self.result(inputs())
        self.assertEqual(result["decision"], "zero_due_to_nonpositive_numerator")
        decomp = result["first_pass"]["source_decomposition"]
        self.assertEqual(decomp["closure"], "PASS")
        self.assertEqual(decomp["bins"][0]["combined_sum"], -2)
        self.assertAlmostEqual(decomp["bins"][0]["combined_error"], math.sqrt(30))
        self.assertEqual(set(decomp["bins"][0]["sources"]), set(audit.SOURCES))

    def test_positive_numerator(self):
        self.assertEqual(self.result(inputs(2.0))["decision"], "positive_staged_amplitude")

    def test_zero_numerator(self):
        self.assertEqual(self.result(inputs(0.0))["decision"], "zero_due_to_nonpositive_numerator")

    def test_failed_support_and_invalid_template(self):
        for message in ("template has zero support inside anchor window", "no valid fit bins",
                        "missing SIMC template shape for pi_n", "bin-count mismatch for pi_n",
                        "axis-range mismatch for pi_n", "non-positive integral for pi_n"):
            data = inputs()
            fit = data[2][0]["fit_result"]
            fit["diagnostics"]["pion"]["staged_pass_history"][0]["pi_n"].update(success=False, message=message)
            fit["H_pion_fit_step_overlays"][0]["fit_diagnostics"] = {}
            self.assertEqual(self.result(data)["decision"], "zero_due_to_invalid_or_unsupported_template")

    def test_missing_diagnostics_does_not_infer_numerator(self):
        data = inputs()
        data[2][0]["fit_result"]["H_pion_fit_step_overlays"][0]["fit_diagnostics"] = {}
        self.assertEqual(self.result(data)["decision"], "not_determined")
        data[2][0]["fit_result"]["H_pion_fit_step_overlays"] = []
        self.assertIn("unavailable_reason", self.result(data)["first_pass"])

    def test_moment_mismatch_or_nonfinite_is_not_classified(self):
        for value in (123., float("nan")):
            data = inputs()
            data[2][0]["fit_result"]["H_pion_fit_step_overlays"][0]["fit_diagnostics"]["weighted_numerator"] = value
            self.assertEqual(self.result(data)["decision"], "not_determined")

    def test_original_solve_survives_postfit_overlay_scaling(self):
        data = inputs(2.0)
        data[2][0]["fit_result"]["H_pion_fit_step_overlays"][0].update(
            amplitude=0., raw_amplitude=2., postfit_scale_factor=0.)
        result = self.result(data)
        self.assertEqual(result["decision"], "positive_staged_amplitude")
        self.assertEqual(result["first_pass"]["amplitude"], 2.)
        self.assertEqual(result["first_pass"]["overlay_amplitude"], 0.)

    def test_bin_edges_and_nommcuts(self):
        records = [dict(source_label="prompt", coefficient=1., adj_MM=x, nommcuts=True)
                   for x in (.879, .88, .909, .91, .939, .94)]
        records.append(dict(source_label="prompt", coefficient=100., adj_MM=.9, nommcuts=False))
        result = audit.decompose_anchor(Histogram((2., 2.), (math.sqrt(2), math.sqrt(2))), records, [1, 2])
        self.assertEqual(result["closure"], "PASS")
        self.assertEqual([r["sources"]["prompt"]["count"] for r in result["bins"]], [2, 2])

    def test_mutated_records_histogram_and_variance_detected(self):
        for kind in ("record", "content", "error"):
            data = inputs()
            h = data[2][0]["H_pion_control"]
            if kind == "record":
                data[0]["_authoritative_pion_control_source_cache"]["by_t"][0]["records"][0]["coefficient"] += 1
            elif kind == "content":
                h.contents[0] += 1
            else:
                h.errors[0] += 1
            self.assertEqual(self.result(data)["first_pass"]["source_decomposition"]["closure"], "FAIL")

    def test_invalid_record_not_silently_dropped(self):
        data = inputs()
        data[0]["_authoritative_pion_control_source_cache"]["by_t"][0]["records"][0]["adj_MM"] = float("nan")
        d = self.result(data)["first_pass"]["source_decomposition"]
        self.assertEqual(d["closure"], "FAIL")
        self.assertEqual(len(d["invalid_records"]), 1)

    def test_multipass_keeps_first_and_adjacent_parent(self):
        data = inputs()
        fit = data[2][0]["fit_result"]
        second = copy.deepcopy(fit["H_pion_fit_step_overlays"][0])
        second.update(pass_index=2, amplitude=2.)
        fit["H_pion_fit_step_overlays"].append(second)
        fit["diagnostics"]["pion"]["staged_pass_history"].append(dict(pi_n=dict(success=True)))
        adjacent = copy.deepcopy(data[2][0]); adjacent.update(t_bin_index=0)
        data[2].append(adjacent)
        result = audit.build_anchor_audit(*data)
        self.assertEqual(len(result["parents"]), 2)
        parent = result["parents"][0]
        self.assertEqual(parent["first_pass"]["amplitude"], 0.)
        self.assertEqual(len(parent["passes"]), 2)

    def test_target_gating_no_io(self):
        for key, value in (("ParticleType", "pion"), ("Q2", "3p0"), ("W", "2p32"), ("EPSSET", "high"), ("EPSSET", "lowe")):
            data = inputs(); data[1][key] = value
            self.assertIsNone(audit.build_anchor_audit(*data))
        data = inputs(); data[0]["phi_setting"] = "Left"
        with patch.object(audit, "write_anchor_audit", side_effect=AssertionError("unexpected IO")):
            audit.observe_anchor_audit(*data, "unused")

    def test_json_nonfinite_explicit(self):
        result = audit.json_safe(dict(x=[float("nan"), float("inf"), None]))
        json.dumps(result, allow_nan=False)
        self.assertTrue(all("unavailable_reason" in x for x in result["x"]))

    def test_atomic_idempotent_collision_safe(self):
        result = audit.build_anchor_audit(*inputs())
        with tempfile.TemporaryDirectory() as directory:
            path = audit.write_anchor_audit(directory, result)
            original = path.read_bytes()
            self.assertEqual(audit.write_anchor_audit(directory, result), path)
            result["attribution"] = "different run"
            with self.assertRaises(FileExistsError):
                audit.write_anchor_audit(directory, result)
            self.assertEqual(path.read_bytes(), original)
            self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_failed_publication_leaves_no_partial_file(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(audit.os, "link", side_effect=OSError("IO failure")):
                with self.assertRaises(OSError):
                    audit.write_anchor_audit(directory, audit.build_anchor_audit(*inputs()))
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_no_mutation_and_diagnostic_failure_isolated(self):
        data = inputs(); before = copy.deepcopy(data)
        with tempfile.TemporaryDirectory() as directory, redirect_stdout(io.StringIO()):
            audit.observe_anchor_audit(*data, directory)
            with patch.object(audit, "write_anchor_audit", side_effect=OSError("denied")):
                audit.observe_anchor_audit(*data, directory)
        self.assertEqual(data[0], before[0])
        self.assertEqual(data[1], before[1])
        self.assertEqual(data[2][0]["fit_result"], before[2][0]["fit_result"])
        self.assertEqual(data[2][0]["H_pion_control"].__dict__, before[2][0]["H_pion_control"].__dict__)

    def test_missing_target_and_cache_explicit(self):
        data = inputs(); data[2].clear()
        self.assertIn("unavailable_reason", audit.build_anchor_audit(*data))
        data = inputs(); data[0].pop("_authoritative_pion_control_source_cache")
        self.assertIn("unavailable_reason", self.result(data)["first_pass"]["source_decomposition"])


if __name__ == "__main__":
    unittest.main()
