"""PyROOT-gated ownership and lifetime regressions for particle subtraction."""

from __future__ import annotations

import gc
import os
import sys
import unittest


ROOT = None
try:
    import ROOT  # type: ignore
except ImportError:
    ROOT = None


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for relative_path in ("src/cuts", "src/utility"):
    path = os.path.join(REPO_ROOT, relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)


@unittest.skipUnless(ROOT is not None, "PyROOT is required for ROOT ownership tests")
class ParticleSubtractionRootOwnershipTests(unittest.TestCase):
    def setUp(self):
        ROOT.gROOT.SetBatch(True)
        from pion_component_subtraction import (
            FIT_OWNED_APPLICATION_PAYLOAD_KEYS,
            assert_component_subtraction_payload_ownership,
        )
        from root_histogram_ownership import (
            clone_root_histogram,
            configure_particle_subtraction_root_ownership_debug,
            get_particle_subtraction_root_ownership_debug_records,
        )

        self.forbidden_keys = FIT_OWNED_APPLICATION_PAYLOAD_KEYS
        self.assert_payload = assert_component_subtraction_payload_ownership
        self.clone = clone_root_histogram
        self.configure_debug = configure_particle_subtraction_root_ownership_debug
        self.debug_records = get_particle_subtraction_root_ownership_debug_records

    @staticmethod
    def _source(name):
        hist = ROOT.TH1D(name, name, 12, 0.8, 1.4)
        hist.Sumw2()
        for index in range(1, hist.GetNbinsX() + 1):
            hist.SetBinContent(index, float(index))
            hist.SetBinError(index, 0.1 * float(index))
        return hist

    def test_detached_clone_preserves_contents_errors_and_reset_statistics(self):
        source = self._source("ownership_source")
        copied = self.clone(source, scope="unit", role="copy", sumw2=True)
        reset = self.clone(source, scope="unit", role="reset", reset=True, sumw2=True)

        self.assertIsNone(copied.GetDirectory())
        self.assertIsNone(reset.GetDirectory())
        self.assertNotEqual(copied.GetName(), reset.GetName())
        for index in range(1, source.GetNbinsX() + 1):
            self.assertAlmostEqual(copied.GetBinContent(index), source.GetBinContent(index))
            self.assertAlmostEqual(copied.GetBinError(index), source.GetBinError(index))
            self.assertEqual(reset.GetBinContent(index), 0.0)
            self.assertEqual(reset.GetBinError(index), 0.0)

        reset.Add(source)
        reset.Scale(0.5)
        for index in range(1, source.GetNbinsX() + 1):
            self.assertAlmostEqual(reset.GetBinContent(index), 0.5 * source.GetBinContent(index))
            self.assertAlmostEqual(reset.GetBinError(index), 0.5 * source.GetBinError(index))

    def test_retained_matrix_survives_gc_with_unique_names_and_debug_records(self):
        self.configure_debug(True)

        def build_retained_matrix():
            source = self._source("matrix_source")
            retained = []
            for t_index in range(4):
                for phi_index in range(8):
                    scope = "t{}_phi{}".format(t_index + 1, phi_index + 1)
                    fit_result = {
                        "H_pi_delta_protected_fit_input": self.clone(
                            source, scope=scope, role="fit_input", sumw2=False
                        ),
                        "H_pi_delta_lambda_gauge": self.clone(
                            source, scope=scope, role="lambda_gauge", sumw2=False
                        ),
                    }
                    application_payload = {
                        "H_pion_weight_vs_MM": self.clone(
                            source, scope=scope, role="weight", reset=True, sumw2=True
                        ),
                        "H_MM_before_pion_subtraction": self.clone(
                            source, scope=scope, role="before", sumw2=False
                        ),
                    }
                    application_payload["H_pion_weight_vs_MM"].Add(source)
                    self.assert_payload(application_payload)
                    retained.append((fit_result, application_payload))
            return retained

        retained = build_retained_matrix()
        gc.collect()
        names = []
        for fit_result, application_payload in retained:
            for histogram in (*fit_result.values(), *application_payload.values()):
                self.assertIsNone(histogram.GetDirectory())
                self.assertGreater(histogram.GetNbinsX(), 0)
                names.append(histogram.GetName())
        self.assertEqual(len(names), len(set(names)))

        records = self.debug_records()
        self.assertEqual(len(records), len(names))
        self.assertTrue(all(record["detached"] is True for record in records))
        self.assertTrue(all(not hasattr(record, "GetName") for record in records))

        canvas = ROOT.TCanvas("ownership_matrix_render", "ownership matrix render", 600, 400)
        retained[0][0]["H_pi_delta_protected_fit_input"].Draw("hist")
        retained[-1][1]["H_pion_weight_vs_MM"].Draw("hist same")
        canvas.Update()
        canvas.Close()

        gc.collect()
        for fit_result, application_payload in retained:
            self.assertGreater(fit_result["H_pi_delta_lambda_gauge"].Integral(), 0.0)
            self.assertGreater(application_payload["H_MM_before_pion_subtraction"].Integral(), 0.0)
        self.configure_debug(False)

    def test_application_payload_rejects_fit_owned_histograms(self):
        self.assert_payload({"H_pion_weight_vs_MM": self._source("application_ok")})
        for key in self.forbidden_keys:
            with self.subTest(key=key):
                with self.assertRaises(AssertionError):
                    self.assert_payload({key: self._source("forbidden_{}".format(key))})
        with self.assertRaises(AssertionError):
            self.assert_payload({"H_pi_delta_protected_fit_total": self._source("forbidden_prefix")})


if __name__ == "__main__":
    unittest.main()
