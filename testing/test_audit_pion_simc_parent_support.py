"""Deterministic arithmetic/provenance tests; no farm or ROOT acceptance."""
import contextlib
import copy
import io
import json
import math
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from testing import audit_pion_simc_parent_support as audit
from src.utility import background_config as source_bg


def config():
    return dict(Q2="4p4", W="2p74", EPSSET="low", ParticleType="kaon", POL="positive",
                tmin=0., tmax=2., mm_min=.8, mm_max=1.4, NumtBins=2, NumPhiBins=2,
                cut_mode="poly", poly_points=[[0., 0.], [1., 0.], [0., 1.]],
                analysis_runtime_config={"config_hash": "fixture"},
                analysis_runtime_config_hash="fixture", source_head="fixture-head")


class JSONTokenBoundaries(unittest.TestCase):
    def assert_number_across_boundary(self, number, marker, offset):
        text = " " * (offset - number.index(marker)) + number + "]"
        self.assertEqual(text[offset], marker)
        tokens = audit._JSONTokens(io.StringIO(text))
        self.assertEqual(tokens.next(), number)
        self.assertEqual(tokens.next(), "]")
        self.assertIsNone(tokens.next())

    def test_decimal_point_at_offset_65535_and_adjacent_boundaries(self):
        for number in ("1.25", "-0.125", "12.345"):
            for offset in (65534, 65535, 65536, 131071):
                with self.subTest(number=number, offset=offset):
                    self.assert_number_across_boundary(number, ".", offset)

    def test_exponent_marker_and_sign_across_read_boundaries(self):
        for number, marker in (("1e12", "e"), ("1E-12", "E"), ("1E-12", "-"),
                               ("1.25e+12", "e"), ("1.25e+12", "+")):
            for offset in (65534, 65535, 65536, 131071):
                with self.subTest(number=number, marker=marker, offset=offset):
                    self.assert_number_across_boundary(number, marker, offset)

    def test_scalar_requires_separator_or_complete_eof(self):
        for text in ("1.", "1e", "1e+", "1e-", "1.25E+", "1x", "truex"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                audit._JSONTokens(io.StringIO(text)).next()
        for text in ("1.25", "1e-12", "true", "null"):
            tokens = audit._JSONTokens(io.StringIO(text))
            self.assertEqual(tokens.next(), text)
            self.assertIsNone(tokens.next())

    def test_numeric_token_limit_remains_one_mib(self):
        with self.assertRaisesRegex(ValueError, "1 MiB"):
            audit._JSONTokens(io.StringIO("0." + "1" * (1024 * 1024))).next()


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
    def test_scan_source_integrity_accepts_resolved_repository_alias(self):
        script = Path(audit.__file__).resolve()
        root = script.parents[1]
        alias = root.with_name(root.name + "-alias") / "testing" / script.name
        resolve = Path.resolve

        def resolve_alias(path, *args, **kwargs):
            return script if path == alias else resolve(path, *args, **kwargs)

        # Model /group -> /u/group without requiring symlink privileges.
        # Stop at the unchanged integrity guard, before any ROOT import.
        with patch.object(audit, "__file__", str(alias)), \
                patch.object(Path, "resolve", autospec=True, side_effect=resolve_alias), \
                patch.object(audit.subprocess, "check_output", return_value="fixture-head\n"), \
                patch.object(audit.subprocess, "run", side_effect=RuntimeError("integrity guard reached")) as guard:
            with self.assertRaisesRegex(RuntimeError, "integrity guard reached"):
                audit.scan(config(), ["Center"])
        guard.assert_called_once_with(
            ["git", "diff", "--quiet", "HEAD", "--", str(Path("testing") / script.name),
             str(Path("src/utility/background_config.py")), str(Path("src/binning/find_bins.py")),
             str(Path("src/cuts/apply_cuts.py")), str(Path("src/cuts/hgcer_hole.py")),
             str(Path("src/utility/pion_component_shapes.py")),
             str(Path("src/cuts/pion_component_subtraction.py")),
             str(Path("src/binning/binning_helpers.py")), str(Path("src/main.py"))], cwd=root)

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


class OriginalAnalysisInput(unittest.TestCase):
    """Small authentic export envelope; selectors stay source-owned, ROOT-free."""
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.source = self.root / "FullAnalysis.json"
        self.bundle = self.root / "transfer"
        self.bundle.mkdir()
        self.report = self.bundle / "fresh.json"
        self.config = config()
        del self.config["source_head"]
        self.config["analysis_runtime_config"] = source_bg.resolve_analysis_runtime_config("4p4", "2p74")
        self.config["analysis_runtime_config_hash"] = self.config["analysis_runtime_config"]["config_hash"]
        self.config.update(particle_subtraction_mode="simc_shape_components",
                           bg_active_profile="no_empirical_residual", simc_tree_name="h10",
                           t_bins=[0., 1., 2.], phi_bins=[-180., 0., 180.],
                           canonical_t_binning=dict(t_edges=[0., 1., 2.], phi_edges=[-180., 0., 180.],
                               canonical_interval_pair_id="fixture-pair", canonical_interval_pair_hash="pair-hash"),
                           particle_subtraction_component_fit_windows={"fixture": [[.88, .94]]})
        self.config["simc_pion_component_files"] = {}
        for setting in ("Center", "Left"):
            paths = self.config["simc_pion_component_files"][setting] = {}
            for component in audit.COMPONENTS:
                path = self.root / (setting + component + ".root")
                path.write_bytes(b"path-existence fixture only; not ROOT evidence")
                paths[component] = str(path)
        self.document = {"inpDict": dict(self.config, _authoritative_pion_control_cache_by_phi={
            "Center": [{"source_head": "decoy", "missmass": 999., "iter_weight": 999.}]}),
            "histlist": [{"inpDict": {"Q2": "decoy"}, "source_head": "decoy"}]}
        self.save()
        self.bg = SimpleNamespace(
            resolve_particle_subtraction_mode=source_bg.resolve_particle_subtraction_mode,
            resolve_analysis_runtime_config=source_bg.resolve_analysis_runtime_config,
            resolve_simc_pion_component_files=source_bg.resolve_simc_pion_component_files)
        self.tcheck = dict(valid=True, t_edges=[0., 1., 2.], validation_rejection_reasons=[])
        self.pcheck = dict(valid=True, phi_edges=[-180., 0., 180.], validation_rejection_reasons=[])
        self.pair = dict(valid=True, canonical_interval_pair_id="fixture-pair",
                         canonical_interval_pair_hash="pair-hash", validation_rejection_reasons=[])
        self.finder = SimpleNamespace(
            _validate_authoritative_t_interval=lambda c, e: self.tcheck,
            _validate_authoritative_phi_interval=lambda c, e: self.pcheck,
            _validate_canonical_interval_pair=lambda t, p: self.pair)

    def save(self):
        self.source.write_text(json.dumps(self.document), encoding="utf-8")

    def argv(self):
        return ["--analysis-json", str(self.source), "--analysis-json-sha256",
                audit.file_identity(self.source)["sha256"], "--settings", "Center", "Left",
                "--children", "--report", str(self.report)]

    def compatible(self, config=None):
        return audit.validate_input_compatibility(config or self.config, self.bg, self.finder,
                                                   historical_input=True)

    def test_extracts_original_selection_without_histlist_or_cache(self):
        selected, producer = audit.read_analysis_config(self.source)
        self.assertEqual(selected, self.config)
        self.assertEqual(producer["recorded_producer_HEAD"], "NOT VERIFIED")
        self.assertEqual(producer["producer_HEAD_status"], "NOT VERIFIED")
        audit.validate_config(selected, "current-diagnostic-head", historical_input=True)
        self.assertNotIn("source_head", selected)
        self.compatible(selected)

    def test_streaming_crosses_chunk_boundaries_and_discards_nonfinite_decoys(self):
        self.document["histlist"] = [{"payload": "x" * 65535 + '\\"\nλ', "unused": math.nan}] * 3
        self.document["inpDict"]["public_override"] = "y" * 65537 + '\\"\nλ'
        self.save()
        selected, _ = audit.read_analysis_config(self.source)
        self.assertEqual(selected["public_override"], self.document["inpDict"]["public_override"])
        self.assertNotIn("histlist", selected)

    def test_combined_json_numeric_boundaries_in_config_and_discarded_cache(self):
        for number, marker in (("1.25", "."), ("1e-12", "e"), ("1.25E+12", "+")):
            for cached in (False, True):
                prefix = '{"inpDict":' + json.dumps(self.config)[:-1]
                prefix += ',"_authoritative_pion_control_cache_by_phi":[' if cached else ',"boundary_value":'
                text = prefix + " " * (65535 - len(prefix) - number.index(marker)) + number
                text += ' ]},"histlist":[]}' if cached else '},"histlist":[]}'
                self.assertEqual(text[65535], marker)
                self.source.write_text(text, encoding="utf-8")
                with self.subTest(number=number, cached=cached):
                    selected, _ = audit.read_analysis_config(self.source)
                    expected = json.loads(text)["inpDict"]
                    expected.pop("_authoritative_pion_control_cache_by_phi", None)
                    self.assertEqual(selected, expected)

    def test_missing_or_invalid_envelope_and_malformed_ignored_data(self):
        for text in ('[]', '{}', '{"inpDict": [], "histlist": []}',
                     '{"inpDict": {}, "histlist": {}}', '{"inpDict": {}}',
                     '{"inpDict": {}, "histlist": [1,]}',
                     '{"inpDict": {}, "histlist": ["bad\\q"]}',
                     '{"inpDict": {}, "histlist": [01]}',
                     '{"inpDict": {}, "inpDict": {}, "histlist": []}',
                     '{"inpDict": {}, "histlist": []} trailing',
                     '{"inpDict": {}, "histlist": ['):
            self.source.write_text(text, encoding="utf-8")
            with self.subTest(text=text), self.assertRaises(ValueError):
                audit.read_analysis_config(self.source)

    def test_fabricated_head_cannot_be_producer_evidence(self):
        argv = self.argv()  # Pin original bytes before someone backfills a HEAD.
        self.document["inpDict"]["source_head"] = "a" * 40
        self.save()
        with patch.object(audit, "scan") as scan, contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(audit.main(argv), 1)
            scan.assert_not_called()
        self.assertFalse(self.report.exists())
        # Even a claimed HEAD in separately checksummed input is only a claim.
        selected, producer = audit.read_analysis_config(self.source)
        self.assertEqual(producer["recorded_producer_HEAD"], "a" * 40)
        self.assertEqual(producer["producer_HEAD_status"], "NOT VERIFIED")
        audit.validate_config(selected, "b" * 40, historical_input=True)
        self.document["inpDict"]["source_head"] = "today's HEAD"
        self.save()
        with self.assertRaisesRegex(ValueError, "producer evidence"):
            audit.read_analysis_config(self.source)

    def test_conflicting_recorded_heads_rejected(self):
        self.document["source_HEAD"] = "a" * 40
        self.document["inpDict"]["source_head"] = "b" * 40
        self.save()
        with self.assertRaisesRegex(ValueError, "conflicting"):
            audit.read_analysis_config(self.source)

    def test_dual_modes_missing_or_crossed_checksums_rejected(self):
        cases = [self.argv() + ["--config", str(self.source)],
                 self.argv() + ["--config-sha256", "extra"],
                 [x for i, x in enumerate(self.argv()) if i not in (2, 3)],
                 ["--config", str(self.source), "--analysis-json-sha256", "wrong",
                  "--settings", "Left", "--report", str(self.report)]]
        for argv in cases:
            with self.subTest(argv=argv), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as exc:
                audit.main(argv)
            self.assertEqual(exc.exception.code, 2)
        self.assertFalse(self.report.exists())

    def test_scope_polygon_runtime_mode_and_profile_fail_closed(self):
        for key, bad in (("Q2", "bad"), ("W", "bad"), ("EPSSET", "high"),
                         ("ParticleType", "pion"), ("poly_points", []),
                         ("analysis_runtime_config_hash", "wrong"), ("analysis_runtime_config", [])):
            selected = dict(self.config, **{key: bad})
            with self.subTest(key=key), self.assertRaises(ValueError):
                audit.validate_config(selected, "current-head", historical_input=True)
        for key, bad in (("analysis_runtime_config", {"config_hash": "fixture", "unexpected": True}),
                         ("particle_subtraction_mode", "single_scale"),
                         ("bg_active_profile", "other"), ("simc_tree_name", "other")):
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.compatible(dict(self.config, **{key: bad}))
        selected = dict(self.config)
        del selected["particle_subtraction_mode"]
        with self.assertRaisesRegex(ValueError, "original production mode"):
            self.compatible(selected)

    def test_canonical_edges_pair_and_sidecar_failures(self):
        for key, bad in (("t_bins", [0, .9, 2]), ("phi_bins", [-180, 1, 180])):
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "sidecar"):
                self.compatible(dict(self.config, **{key: bad}))
        for key, bad in (("t_edges", [0, .9, 2]), ("phi_edges", [-180, 1, 180]),
                         ("canonical_interval_pair_id", "stale"),
                         ("canonical_interval_pair_hash", "stale")):
            selected = copy.deepcopy(self.config)
            selected["canonical_t_binning"][key] = bad
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.compatible(selected)
        for sidecar in (self.tcheck, self.pcheck, self.pair):
            sidecar.update(valid=False, validation_rejection_reasons=["fixture_stale_sidecar"])
            with self.assertRaisesRegex(ValueError, "fixture_stale_sidecar"):
                self.compatible()
            sidecar["valid"] = True
        selected = copy.deepcopy(self.config)
        selected["canonical_t_binning"]["authoritative_metadata"] = {"canonical_interval_pair_hash": "stale"}
        with self.assertRaisesRegex(ValueError, "sidecar"):
            self.compatible(selected)

    def test_all_six_original_paths_required_without_fallback(self):
        for setting in ("Center", "Left"):
            for component in audit.COMPONENTS:
                selected = copy.deepcopy(self.config)
                del selected["simc_pion_component_files"][setting][component]
                with self.subTest(setting=setting, component=component), self.assertRaisesRegex(ValueError, "original SIMC path"):
                    self.compatible(selected)
        selected = copy.deepcopy(self.config)
        selected["simc_pion_component_files"]["Left"]["pi_delta"] = str(self.root / "absent.root")
        with self.assertRaisesRegex(ValueError, "original SIMC input missing"):
            self.compatible(selected)
        self.bg.resolve_simc_pion_component_files = lambda c, s: dict.fromkeys(audit.COMPONENTS, "fallback")
        with self.assertRaisesRegex(ValueError, "differs from original"):
            self.compatible()

    def test_cli_direct_handoff_preserves_input_and_separates_identities(self):
        before = self.source.read_bytes()
        def scan(selected, settings, children, evidence, *, historical_input):
            self.assertTrue(historical_input)
            self.assertEqual(selected, self.config)
            self.assertEqual(settings, ["Center", "Left"])
            self.assertTrue(children)
            self.assertIsNone(evidence)
            self.compatible(selected)
            return {"diagnostic_checkout": {"HEAD": "current-head"}, "scientific_adequacy": "NOT VERIFIED"}
        with patch.object(audit, "scan", side_effect=scan), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(audit.main(self.argv()), 0)
        result = json.loads(self.report.read_text())
        self.assertEqual(result["input_mode"], "analysis-json")
        self.assertEqual(result["original_analysis_input"]["sha256"], audit.file_identity(self.source)["sha256"])
        self.assertEqual(result["original_analysis_input"]["bytes"], len(before))
        self.assertEqual(result["original_analysis_input"]["recorded_producer_HEAD"], "NOT VERIFIED")
        self.assertEqual(result["diagnostic_checkout"]["HEAD"], "current-head")
        self.assertEqual(self.source.read_bytes(), before)

    def test_optional_evidence_passed_and_identity_rechecked(self):
        path = self.root / "center.json"
        evidence = {"parents": [{"pion_parent_id": "fixture", "t_edges": [0, 1]}]}
        path.write_text(json.dumps(evidence), encoding="utf-8")
        argv = self.argv() + ["--center-evidence", str(path)]
        def scan(config, settings, children, document, **kwargs):
            self.assertEqual(document, evidence)
            path.write_text("changed", encoding="utf-8")
            return {}
        with patch.object(audit, "scan", side_effect=scan), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(audit.main(argv), 1)
        self.assertFalse(self.report.exists())

    def test_token_and_nesting_bounds_fail_without_output(self):
        self.source.write_text('{"inpDict": {}, "histlist": ["' + 'x' * (1024 * 1024 + 1) + '"]}',
                               encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "1 MiB"):
            audit.read_analysis_config(self.source)
        self.source.write_text('{"inpDict": {}, "histlist": ' + '[' * 130 + '0' + ']' * 130 + '}',
                               encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "depth"):
            audit.read_analysis_config(self.source)
        self.assertFalse(self.report.exists())

    def test_failed_scan_changed_input_and_existing_report_leave_no_new_evidence(self):
        for failure in (ValueError("parity failure"), OSError("input missing")):
            with patch.object(audit, "scan", side_effect=failure), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(audit.main(self.argv()), 1)
            self.assertFalse(self.report.exists())
        def mutate(*args, **kwargs):
            self.source.write_text("changed", encoding="utf-8")
            return {}
        argv = self.argv()
        with patch.object(audit, "scan", side_effect=mutate), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(audit.main(argv), 1)
        self.assertFalse(self.report.exists())
        self.save()
        def race(*args, **kwargs):
            self.report.write_text("prior evidence", encoding="utf-8")
            return {}
        with patch.object(audit, "scan", side_effect=race), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(audit.main(self.argv()), 1)
        self.assertEqual(self.report.read_text(), "prior evidence")

    def test_serialization_and_write_failures_do_not_leave_partial_report(self):
        with self.assertRaises(ValueError):
            audit.write_report(self.report, {"bad": math.nan})
        self.assertFalse(self.report.exists())
        stream = self.report.open("x", encoding="utf-8")
        class FailingWrite:
            def __enter__(self): return self
            def __exit__(self, *args): stream.close()
            def write(self, payload):
                stream.write(payload[:5])
                raise OSError("fixture disk failure")
        with patch.object(Path, "open", return_value=FailingWrite()), self.assertRaises(OSError):
            audit.write_report(self.report, {"complete": True})
        self.assertFalse(self.report.exists())
        with self.assertRaises(FileNotFoundError):
            audit.write_report(self.root / "missing-directory" / "report.json", {})
        self.assertFalse((self.root / "missing-directory").exists())


if __name__ == "__main__":
    unittest.main()
