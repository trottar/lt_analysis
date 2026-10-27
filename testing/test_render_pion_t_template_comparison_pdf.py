"""Stored presentation replay tests; no ROOT/farm/PDF acceptance."""
from copy import deepcopy
import json
import math
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch

from testing import render_pion_t_template_comparison_pdf as replay
from testing.test_pion_t_template_comparison import Hist

SHA = "a"*40


def fixture():
    value = replay.comparison.histogram(Hist(values=(3., 2., 6., 4.), errors=(.5, 1., 2., .75)))
    fit = {"H_pion_control_input": deepcopy(value), "H_kaon_nosub_input": deepcopy(value)}
    for side in ("pion", "kaon"):
        for suffix in ("", "_refined"):
            fit["H_{}_fit_total{}".format(side, suffix)] = deepcopy(value)
            fit["H_{}_fit_pi_n_scaled{}".format(side, suffix)] = deepcopy(value)
            fit["H_{}_fit_k_sigma0_scaled{}".format(side, suffix)] = {"status": "unavailable"}
    return {"schema_version": replay.comparison.SCHEMA, "source_HEAD": SHA, "status": "available",
        "diagnostic_only": True, "production_applied": False, "scientifically_accepted": False,
        "source_identity_rejection_reasons": [], "phi_setting": "Center",
        "setting": {"Q2": "4p4", "W": "2p74", "EPSSET": "low", "ParticleType": "kaon"},
        "parents": [{"t_bin_index": i, "status": "available", "authoritative_parent_unchanged": True,
                     "baseline_A": deepcopy(fit), "candidate_B": deepcopy(fit)} for i in range(3)]}


class StoredHist(Hist):
    def __init__(self, name, title, n, edges):
        super().__init__([0.]*(n+2), [0.]*(n+2), edges)
        class Sumw2(list):
            def SetAt(self, value, index): self[index] = value
        self.variances = Sumw2([0.]*(n+2))
    def GetSumw2(self): return self.variances
    def GetBinError(self, i): return math.sqrt(self.variances[i])
    def SetBinError(self, i, value): self.variances[i] = value*value


class ReplayTests(unittest.TestCase):
    def test_restore_flow_contents_errors_and_sumw2_without_mutating_payload(self):
        value = fixture()["parents"][0]["baseline_A"]["H_pion_control_input"]
        before = deepcopy(value)
        hist = replay.restore_histogram(value, "fixture", types.SimpleNamespace(TH1D=StoredHist))
        retained = replay.comparison.histogram(hist)
        self.assertEqual(retained, value)
        self.assertEqual(value, before)
        self.assertEqual(hist.GetBinContent(0), 3.)
        self.assertEqual(hist.GetBinContent(3), 4.)

    def test_all_required_stored_inputs_and_identity_are_validated(self):
        document = fixture()
        before = deepcopy(document)
        rows = replay.validate_document(document, SHA)
        self.assertEqual(len(rows), 3)
        self.assertEqual(document, before)
        for mutation in (
            lambda d: d.update(source_HEAD="b"*40),
            lambda d: d.update(production_applied=True),
            lambda d: d.update(phi_setting="Right"),
            lambda d: d["parents"].pop(),
            lambda d: d["parents"][0]["baseline_A"].pop("H_pion_fit_total_refined"),
            lambda d: d["parents"][0]["candidate_B"]["H_pion_control_input"]["values"].pop(),
        ):
            bad = deepcopy(document); mutation(bad)
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                replay.validate_document(bad, SHA)

    def test_invalid_histogram_geometry_numbers_errors_and_hashes_fail(self):
        original = fixture()["parents"][0]["baseline_A"]["H_pion_control_input"]
        for field, value in (("edges", [0., 0., 2.]), ("values", [0., 1.]),
                             ("errors", [0., -1., 2., 0.]), ("values", [0., float("nan"), 1., 0.]),
                             ("axis_signature", "0"*64), ("content_error_fingerprint", "0"*64),
                             ("stored_sumw2", [0., 1., 4., 0.])):
            bad = deepcopy(original); bad[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                replay.validate_histogram(bad)

    def test_new_output_replay_records_separate_input_and_renderer_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory); source = base/"input.json"
            source.write_text(json.dumps(fixture()))
            before = source.read_bytes()
            root = types.SimpleNamespace(TH1D=StoredHist, gROOT=Mock())
            def render(path, rows, title):
                self.assertEqual(len(rows), 3)
                path.write_bytes(b"%PDF-1.4\nfixture only, not ROOT\n")
            with patch.dict(sys.modules, {"ROOT": root}), patch.object(replay.comparison, "render_comparison_pdf", side_effect=render):
                report = replay.replay(source, base/"new.pdf", base/"report.json", SHA)
            self.assertEqual(source.read_bytes(), before)
            self.assertEqual(report["scientific_input"]["source_HEAD"], SHA)
            self.assertNotEqual(report["renderer"]["observed_git_HEAD"], SHA)
            self.assertEqual(report["expected_pages"], 6)
            self.assertFalse(report["refitted"])
            self.assertFalse(report["pdf_legibility_verified"])
            with self.assertRaisesRegex(ValueError, "new files"):
                replay.replay(source, base/"new.pdf", base/"report.json", SHA)

    def test_paths_and_missing_pdf_fail_without_report_acceptance(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory); source=base/"input.json"
            source.write_text(json.dumps(fixture()))
            with self.assertRaisesRegex(ValueError, "collision"):
                replay.replay(source, base/"new.pdf", source, SHA)
            with self.assertRaisesRegex(ValueError, "outside source"):
                replay.replay(source, replay.REPO/"new.pdf", replay.REPO/"report.json", SHA)
            root=types.SimpleNamespace(TH1D=StoredHist, gROOT=Mock())
            with patch.dict(sys.modules, {"ROOT": root}), patch.object(replay.comparison, "render_comparison_pdf"):
                with self.assertRaisesRegex(ValueError, "PDF missing"):
                    replay.replay(source, base/"new.pdf", base/"report.json", SHA)
            self.assertFalse((base/"report.json").exists())

    def test_invalid_document_stops_before_root_or_renderer(self):
        with tempfile.TemporaryDirectory() as directory:
            base=Path(directory); source=base/"input.json"
            source.write_text(json.dumps(fixture()))
            with patch.object(replay.comparison, "render_comparison_pdf") as render:
                with self.assertRaisesRegex(ValueError, "identity differs"):
                    replay.replay(source, base/"new.pdf", base/"report.json", "b"*40)
            render.assert_not_called()
            self.assertFalse((base/"new.pdf").exists())


if __name__ == "__main__":
    unittest.main()
