"""Deterministic arithmetic/provenance tests; no farm or ROOT acceptance."""
import contextlib
import io
import json
import math
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from testing import audit_pion_simc_parent_support as audit


def config():
    return dict(Q2="4p4", W="2p74", EPSSET="low", ParticleType="kaon", POL="positive",
                tmin=0., tmax=2., mm_min=.8, mm_max=1.4, NumtBins=2, NumPhiBins=2,
                cut_mode="poly", poly_points=[[0., 0.], [1., 0.], [0., 1.]],
                analysis_runtime_config={"config_hash": "fixture"},
                analysis_runtime_config_hash="fixture", source_head="fixture-head")


class MomentsAndAxes(unittest.TestCase):
    def test_mm_flows_and_exact_boundaries(self):
        h = audit.Histogram([0., 1., 2.])
        for value in [-1., 0., math.nextafter(1., 0.), 1., math.nextafter(2., 0.), 2., 3.]:
            h.fill(value, 1.)
        r = h.report()
        self.assertEqual([r["underflow"]["count"], *[b["count"] for b in r["normal_bins"]],
                          r["overflow"]["count"]], [1, 2, 2, 2])
        self.assertEqual(r["normalization_integral"], 4.)
        self.assertEqual(r["normalized_content"], [.5, .5])

    def test_axis_validation(self):
        for axis in ([], [1], [0, 0], [2, 1], [0, math.nan], [0, math.inf]):
            with self.subTest(axis=axis), self.assertRaises(ValueError):
                audit.Histogram(axis)

    def test_t_phi_upper_edge_excluded(self):
        p = audit.Population([0, 2], [0, 1, 2], [-180, 0, 180])
        for t, phi in [(0, -180), (1, 0), (2, 0), (1, 180), (-1, 0)]:
            p.fill(1., t, 1., phi)
        r = p.report()
        self.assertEqual(r["unassigned_t"]["total"]["count"], 2)
        self.assertEqual(r["parents"][0]["children"][0]["histogram"]["total"]["count"], 1)
        self.assertEqual(r["parents"][1]["phi_unassigned"]["total"]["count"], 1)
        self.assertEqual(r["parents"][1]["children"][1]["histogram"]["total"]["count"], 1)

    def test_signed_random_like_moments(self):
        h = audit.Histogram([0, 2])
        for weight in (2., -1., -.5, .5):
            h.fill(1., weight)
        r = h.report()["total"]
        self.assertEqual(r["count"], 4)
        self.assertEqual((r["sum_w"], r["sum_abs_w"], r["sum_w2"]), (1., 4., 5.5))
        self.assertAlmostEqual(r["neff_signed"], 1 / 5.5)
        self.assertAlmostEqual(r["neff_absolute"], 16 / 5.5)
        self.assertEqual(r["cancellation_fraction"], .75)

    def test_zero_and_negative_integral_never_normalized(self):
        for weights in ([], [1, -1], [-2]):
            h = audit.Histogram([0, 2])
            for w in weights:
                h.fill(1, w)
            r = h.report()
            self.assertFalse(r["normalization_meaningful"])
            self.assertIsNone(r["normalized_content"])
        empty = audit.Histogram([0, 2]).report()["total"]
        self.assertIsNone(empty["neff_absolute"])
        self.assertIsNone(empty["cancellation_fraction"])

    def test_finite_rejection_and_overflow(self):
        for mm, w in ((math.nan, 1), (math.inf, 1), (1, math.nan), (1, math.inf), (1, 1e308)):
            with self.subTest(mm=mm, w=w), self.assertRaises(ValueError):
                audit.Histogram([0, 2]).fill(mm, w)

    def test_invalid_coordinates_accounted_without_losing_mm(self):
        p = audit.Population([0, 2], [0, 2], [-180, 180])
        p.fill(1, math.nan, -1, 0)
        p.fill(1, 1, 2, math.inf)
        r = p.report()
        self.assertEqual((r["invalid_t"], r["invalid_phi"]), (1, 1))
        self.assertEqual(r["setting"]["total"]["count"], 2)
        self.assertEqual(r["unassigned_t"]["total"]["sum_w"], -1)
        self.assertEqual(r["parents"][0]["phi_unassigned"]["total"]["sum_w"], 2)

    def test_regions_event_bounds_and_overlap_counted_once(self):
        h = audit.Histogram([0, 1, 2], {"anchor": [[.5, 1], [1, 1.5]]})
        for mm in (.5, 1, 1.5, 1.6):
            h.fill(mm, 2)
        self.assertEqual(h.report()["regions"]["anchor"]["count"], 3)
        self.assertEqual(h.report()["regions"]["anchor"]["sum_w2"], 12)
        with self.assertRaises(ValueError):
            audit.Histogram([0, 1], {"bad": [[2, 1]]})


class ClosureAndShapes(unittest.TestCase):
    def population(self):
        p = audit.Population([0, 1, 2], [0, 1, 2], [-180, 0, 180], {"anchor": [[.5, 1.5]]})
        for mm, t, phi, w in [(.5, .5, -90, 3), (1.5, 1.5, 90, 2),
                              (2., .5, 180, -1), (-1., 3., 0, .5), (1., .5, 90, -.5)]:
            p.fill(mm, t, w, phi)
        return p

    def test_full_bin_region_and_child_closure(self):
        r = self.population().report()
        self.assertEqual(r["setting"]["total"]["count"], 5)
        self.assertEqual(r["unassigned_t"]["underflow"]["sum_w"], .5)
        for parent in r["parents"]:
            for key in ("child_closure", "complement_closure"):
                self.assertTrue(all(all(v == 0 for v in row.values())
                                    for row in parent[key]["bin_then_region_residuals"]))
        self.assertEqual(r["parents"][0]["disjoint_complement"]["total"]["sum_w"], 2.5)

    def test_closure_failure_detects_sumw2_corruption(self):
        p = self.population()
        p.parents[0].bins[1]["sum_w2"] += 1
        with self.assertRaisesRegex(ValueError, "sum_w2"):
            p.report()

    def test_zero_children_reported(self):
        p = audit.Population([0, 2], [0, 1, 2], [-180, 0, 180])
        p.fill(1, .5, 1, -90)
        r = p.report()
        self.assertEqual(r["parents"][1]["children"][1]["histogram"]["total"]["count"], 0)
        self.assertEqual(r["parents"][1]["parent_vs_complement"]["status"], "NOT VERIFIED")
        self.assertEqual(r["parents"][0]["histogram"]["scientific_adequacy"], "NOT VERIFIED")

    def test_no_phi_arrays_fabricated(self):
        p = audit.Population([0, 2], [0, 1])
        p.fill(1, .5, 1)
        r = p.report()
        self.assertEqual(r["phi_status"], "NOT VERIFIED")
        self.assertNotIn("children", r["parents"][0])

    def test_disjoint_shape_variances(self):
        a, b = audit.Histogram([0, 1, 2]), audit.Histogram([0, 1, 2])
        for mm in (.5, .5, 1.5, 1.5):
            a.fill(mm, 1)
            b.fill(mm, 1)
        r = audit.shape_comparison(a, b)
        self.assertEqual(r["bins"][0]["left_variance"], .0625)
        self.assertEqual(r["bins"][0]["difference_variance"], .125)
        self.assertEqual(r["bins"][0]["shared_covariance"], 0)
        self.assertIsNone(r["p_value"])

    def test_inclusive_identical_population_covariance(self):
        h = audit.Histogram([0, 1, 2])
        h.fill(.5, 1)
        h.fill(1.5, 2)
        r = audit.shape_comparison(h, h, shared=True)
        self.assertTrue(all(row["difference_variance"] == 0 for row in r["bins"]))
        self.assertGreater(r["bins"][0]["shared_covariance"], 0)

    def test_axis_mismatch_and_unsupported_shapes(self):
        h = audit.Histogram([0, 1])
        with self.assertRaises(ValueError):
            audit.shape_comparison(h, audit.Histogram([0, 2]))
        self.assertEqual(audit.shape_comparison(h, h)["status"], "NOT VERIFIED")

    def test_finder_is_used_for_children(self):
        calls = []
        def finder(t, phi, te, pe):
            calls.append((t, phi, te, pe))
            return 0, None
        p = audit.Population([0, 2], [0, 1], [-180, 180])
        p.fill(1, .5, 1, 0, finder)
        self.assertEqual(len(calls), 1)
        self.assertEqual(p.report()["parents"][0]["phi_unassigned"]["total"]["count"], 1)


class ProvenanceAndCLI(unittest.TestCase):
    def test_population_selection_default_weight_hole_and_phi_units(self):
        def event(mm=1., **kwargs):
            return SimpleNamespace(missmass=mm, t=-.5, phipq=math.pi / 2,
                                   phgcer_x_det=0, phgcer_y_det=0, **kwargs)
        events = [event(iter_weight=2), event(1.8, iter_weight=-.5), event(),
                  event(iter_weight=math.nan), event(iter_weight="invalid"),
                  event(1.2, iter_weight=99, in_hole=True), event(.9, reject=True),
                  event(math.nan)]
        full = audit.Population([0, 1, 2], [0, 1], [-180, 0, 180])
        mmcut = audit.Population([0, 1, 2], [0, 1], [-180, 0, 180])
        calls, mirror = [], []
        def subcuts(e, lo, hi):
            calls.append((lo, hi))
            return not getattr(e, "reject", False)
        def cuts(e, lo, hi):
            return not getattr(e, "reject", False) and lo <= e.missmass <= hi
        def finder(t, phi, te, pe):
            return audit.bin_index(t, te), audit.bin_index(phi, pe)
        hole = SimpleNamespace(IsInside=lambda x, y: False)
        # Give the hole distinct detector coordinates without changing its API.
        events[5].phgcer_x_det = 99
        hole.IsInside = lambda x, y: x == 99
        counts = audit.accumulate_events(events, config(), full, mmcut, hole, cuts, subcuts,
                                         finder, SimpleNamespace(Fill=lambda mm, w: mirror.append((mm, w))))
        self.assertEqual(counts, dict(seen=8, invalid_weight=1, nonfinite_weight=1,
                                     selected_invalid_mm=1, rejected_cuts_or_hole=2))
        self.assertEqual(len(calls), 6)
        self.assertEqual(mirror, [(1., 2.), (1.8, -.5), (1., 1.)])
        self.assertEqual(full.report()["setting"]["total"]["sum_w"], 2.5)
        self.assertEqual(mmcut.report()["setting"]["total"]["sum_w"], 3.)
        self.assertAlmostEqual(full.report()["coverage"]["phi_min"], 90.)
        self.assertEqual(full.report()["parents"][0]["children"][1]["histogram"]["total"]["count"], 3)

    def test_config_identity_and_missing_prerequisites(self):
        audit.validate_config(config(), "fixture-head")
        for key in config():
            bad = config()
            del bad[key]
            with self.subTest(key=key), self.assertRaises(ValueError):
                audit.validate_config(bad, "fixture-head")
        with self.assertRaisesRegex(ValueError, "HEAD"):
            audit.validate_config(config(), "other-head")
        bad = config()
        bad["analysis_runtime_config_hash"] = "other"
        with self.assertRaisesRegex(ValueError, "identity"):
            audit.validate_config(bad, "fixture-head")

    def test_out_of_scope_and_invalid_config(self):
        for key, value in (("Q2", "2p1"), ("EPSSET", "high"), ("poly_points", []),
                           ("tmax", math.inf), ("poly_points", [[0, 0], [0, 1], [math.nan, 2]])):
            bad = config()
            bad[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                audit.validate_config(bad, "fixture-head")

    def test_branches_required_phi_optional_weight_default(self):
        names = {"missmass", "t", "ssdelta", "ssxptar", "ssyptar", "hsdelta", "hsxptar",
                 "hsyptar", "Q2", "W", "phgcer_x_det", "phgcer_y_det"}
        audit.validate_branches(names, phi=False)
        with self.assertRaisesRegex(ValueError, "phipq"):
            audit.validate_branches(names, phi=True)
        audit.validate_branches(names | {"phipq"}, phi=True)
        with self.assertRaisesRegex(ValueError, "missmass"):
            audit.validate_branches(names - {"missmass"}, phi=False)

    def test_optional_center_identity_check(self):
        identities = [{"pion_parent_id": "test", "t_edges": [.1, .2]}]
        audit.cross_check_evidence({"parents": identities}, identities)
        for document in ({}, {"parents": [{"pion_parent_id": "test", "t_edges": [.2, .3]}]}):
            with self.assertRaises(ValueError):
                audit.cross_check_evidence(document, identities)

    def test_production_parity_and_hard_mismatch(self):
        h = audit.Histogram([0, 1, 2])
        h.fill(.5, 2)
        h.fill(2, 1)
        class RootShape:
            def GetNbinsX(self): return 2
            def GetXaxis(self): return self
            def GetBinLowEdge(self, i): return [0, 1, 2][i - 1]
            def GetBinContent(self, i): return [0, 1, 0, .5][i]
            def GetBinError(self, i): return [0, 1, 0, .5][i]
        d = dict(fallback_used=False, normalized=True, weighted_integral_before_norm=2.)
        audit.production_parity(h, RootShape(), d)
        with self.assertRaisesRegex(ValueError, "fallback"):
            audit.production_parity(h, RootShape(), dict(d, fallback_used=True))
        h.fill(1.5, 1)
        with self.assertRaisesRegex(ValueError, "integral mismatch"):
            audit.production_parity(h, RootShape(), d)

    def test_root_free_cli_failure_does_not_write_or_overwrite(self):
        with tempfile.TemporaryDirectory() as td:
            source, report = Path(td) / "config.json", Path(td) / "report.json"
            source.write_text(json.dumps(config()), encoding="utf-8")
            sha = audit.file_identity(source)["sha256"]
            argv = ["--config", str(source), "--config-sha256", sha, "--settings", "Left", "--report", str(report)]
            with patch.object(audit, "scan", side_effect=ValueError("missing provenance")), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(audit.main(argv), 1)
            self.assertFalse(report.exists())
            report.write_text("preserved", encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(audit.main(argv), 1)
            self.assertEqual(report.read_text(), "preserved")

    def test_cli_sha_failure_and_success_serialization(self):
        with tempfile.TemporaryDirectory() as td:
            source, report = Path(td) / "config.json", Path(td) / "report.json"
            source.write_text(json.dumps(config()), encoding="utf-8")
            argv = ["--config", str(source), "--config-sha256", "wrong", "--settings", "Center", "--report", str(report)]
            with patch.object(audit, "scan") as scan, contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(audit.main(argv), 1)
                scan.assert_not_called()
            argv[3] = audit.file_identity(source)["sha256"]
            with patch.object(audit, "scan", return_value={"schema": audit.SCHEMA, "scientific_adequacy": "NOT VERIFIED"}):
                self.assertEqual(audit.main(argv), 0)
            self.assertEqual(json.loads(report.read_text())["scientific_adequacy"], "NOT VERIFIED")

    def test_synthetic_json_is_finite_and_not_acceptance(self):
        p = audit.Population([0, 1, 2], [0, 1, 2])
        p.fill(.5, .5, -1)
        p.fill(.5, .5, 1)
        result = json.loads(audit.encoded(p.report()))
        self.assertEqual(result["setting"]["scientific_adequacy"], "NOT VERIFIED")
        self.assertEqual(result["setting"]["total"]["cancellation_fraction"], 1)


if __name__ == "__main__":
    unittest.main()
