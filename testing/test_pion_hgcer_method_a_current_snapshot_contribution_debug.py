"""Coherent raw-byte fixtures and adversarial snapshot/owner checks; no farm.

The new live gate is never patched to accept synthetic authority. Fixtures use
the unchanged F.3/F.4 builders and actual serialized-byte digests. Historical
pins, production authority and scientific validators remain untouched.
"""

from contextlib import redirect_stdout, redirect_stderr
import copy
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src" / "cuts"), str(ROOT / "testing")]
import pion_hgcer_method_a_current_snapshot_contribution_debug as debug
import analyze_pion_hgcer_method_a_current_snapshot_contribution_debug as cli
import test_pion_hgcer_method_a_current_lineage_diagnostic as fixtures
import pion_hgcer_method_a_acceptance_contract as f1


class SnapshotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixtures.DiagnosticTests.setUpClass()
        cls.f = fixtures.DiagnosticTests
        cls.hashes = {**cls.f.hashes, "f3": cls.f.map_sha, "f4": cls.f.correction_sha}
        cls.snapshot, cls.parents = cls.reproduce()
        cls.names = {s: s+".json" for s in debug.ROLES}
        cls.provenance = {"git_head": "0"*40, "git_branch": "test", "worktree_dirty": False,
                          "input_basenames": cls.names}
        cls.diagnostic = debug.build_diagnostic(cls.snapshot, cls.parents, provenance=cls.provenance)

    @classmethod
    def reproduce(cls, **overrides):
        args = {"f1_artifacts": cls.f.inputs, "f3_artifact": cls.f.map,
                "f4_artifact": cls.f.correction, "hashes": cls.hashes}
        args.update(overrides)
        return debug.reproduce_snapshot(**args)

    def write_locator(self, root):
        payloads = {s: a for s, a in zip(debug.SETTINGS, self.f.inputs)}
        payloads.update(f3=self.f.map, f4=self.f.correction)
        for role, value in payloads.items():
            (root/self.names[role]).write_bytes(fixtures.encoded(value))
        manifest = {"f1": {s: self.names[s] for s in debug.SETTINGS},
                    "f3": self.names["f3"], "f4": self.names["f4"]}
        path = root/"inputs.json"
        path.write_text(json.dumps(manifest))
        return path, manifest

    def args(self, root, locator):
        return ["--inputs-manifest", str(locator), "--expected-snapshot-fingerprint", self.snapshot["fingerprint"],
                "--outdir", str(root/"artifacts"), "--bundle-dir", str(root/"bundles")]

    def prepare(self, root):
        (root/"artifacts").mkdir(); (root/"bundles").mkdir()
        return self.write_locator(root)[0]

    def test_coherent_live_gate_without_pin_patches_and_preservation(self):
        before = copy.deepcopy((self.f.inputs, self.f.map, self.f.correction,
            debug.aggregate.F1_HASHES, debug.aggregate.RECORDED_PARENTS,
            debug.lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC))
        snapshot, parents = self.reproduce()
        self.assertEqual((snapshot, parents), (self.snapshot, self.parents))
        self.assertEqual(snapshot["input_file_hashes"], self.hashes)
        self.assertEqual((self.f.inputs, self.f.map, self.f.correction,
            debug.aggregate.F1_HASHES, debug.aggregate.RECORDED_PARENTS,
            debug.lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC), before)
        self.assertTrue(debug.validate_diagnostic(self.diagnostic))
        self.assertEqual(self.diagnostic["canonical_owner_equivalence"], "BLOCKED")
        self.assertFalse(self.diagnostic["canonical_owner_equivalence_evaluated"])
        self.assertNotIn("accepted_authority_match", cli.json_bytes(self.diagnostic).decode())

    def test_raw_f1_binding_including_observed_center_highe_mismatch(self):
        for role in debug.SETTINGS:
            hashes = {**self.hashes, role: "40d9175bc1afaaa5cacc965c5ece88b3de08aacd4d20a07ea9bf248e3c21ef96"}
            with self.subTest(role=role), mock.patch.object(debug.f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data") as builder:
                with self.assertRaisesRegex(debug.ContributionDebugError, "f3_.*mismatch"):
                    self.reproduce(hashes=hashes)
                builder.assert_not_called()

    def test_f3_raw_binding_checked_before_stored_record_use(self):
        with mock.patch.object(debug.f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data") as builder:
            with self.assertRaisesRegex(debug.ContributionDebugError, "f3_runtime_authority_source_file_sha256_mismatch"):
                self.reproduce(hashes={**self.hashes, "f3": "a"*64})
            builder.assert_not_called()

    def test_tampered_schemas_wrappers_algorithm_and_body_fingerprints(self):
        for role, path in (("f3", ("schema_version",)), ("f3", ("artifact_fingerprint",)),
                           ("f3", ("acceptance_map", "fingerprint")), ("f3", ("acceptance_map", "algorithm_fingerprint")),
                           ("f4", ("schema_version",)), ("f4", ("artifact_fingerprint",)),
                           ("f4", ("correction", "fingerprint")), ("f4", ("production_objects_mutated",))):
            altered = copy.deepcopy(self.f.map if role == "f3" else self.f.correction)
            node = altered
            for key in path[:-1]: node = node[key]
            node[path[-1]] = "tampered"
            with self.subTest(role=role, path=path), self.assertRaises(debug.ContributionDebugError):
                self.reproduce(**{role+"_artifact": altered})

    def test_f4_f1_f3_links_cannot_be_self_asserted(self):
        for field in ("f3_source_file_sha256", "f3_map_fingerprint", "f3_algorithm_fingerprint", "f3_artifact_fingerprint"):
            altered = copy.deepcopy(self.f.correction)
            altered["correction"][field] = "a"*64
            with self.subTest(field=field), self.assertRaisesRegex(debug.ContributionDebugError, "snapshot_f4_link_mismatch"):
                self.reproduce(f4_artifact=altered)
        altered = copy.deepcopy(self.f.correction)
        altered["correction"]["input_fingerprints"][0]["source_file_sha256"] = "a"*64
        with self.assertRaisesRegex(debug.ContributionDebugError, "snapshot_f4_f1_binding_mismatch"):
            self.reproduce(f4_artifact=altered)
        altered = copy.deepcopy(self.f.correction)
        altered["correction"]["f3_runtime_authority"]["accepted"]["map_fingerprint"] = "a"*64
        with self.assertRaisesRegex(debug.ContributionDebugError, "f3_runtime_authority_map_fingerprint_mismatch"):
            self.reproduce(f4_artifact=altered)

    def test_exact_f4_equality_includes_float_type_and_unknown_fields(self):
        shared = debug.f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data
        for mode in ("float", "type", "unknown"):
            def change(*args, **kwargs):
                body, reviews = shared(*args, **kwargs)
                if mode == "float":
                    p = body["parents"][0]
                    p["baseline_parent_sum"] = math.nextafter(p["baseline_parent_sum"], math.inf)
                elif mode == "type": body["production_objects_mutated"] = 0
                else: body["unknown_scientific_field"] = 1
                return body, reviews
            with self.subTest(mode=mode), mock.patch.object(debug.f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", side_effect=change):
                with self.assertRaisesRegex(debug.ContributionDebugError, "exact_reproduction_mismatch"):
                    self.reproduce()

    def test_missing_parent_review_factor_phi_and_row_mismatch(self):
        shared = debug.f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data
        for mode in ("missing", "duplicate", "factor", "phi"):
            def change(*args, **kwargs):
                body, reviews = shared(*args, **kwargs)
                if mode == "missing": reviews.pop()
                elif mode == "duplicate": reviews.append(reviews[0])
                elif mode == "factor": reviews[0]["correction_factors"] = reviews[0]["correction_factors"][:-1]
                else: reviews[0]["correction_factors"][0] *= 2
                return body, reviews
            with self.subTest(mode=mode), mock.patch.object(debug.f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data", side_effect=change):
                with self.assertRaises(debug.ContributionDebugError): self.reproduce()
        altered = copy.deepcopy(self.f.inputs)
        altered[0]["contract"]["application_records"].append(copy.deepcopy(altered[0]["contract"]["application_records"][0]))
        with self.assertRaises(debug.ContributionDebugError): self.reproduce(f1_artifacts=altered)

    def test_resealed_joint_moment_source_tail_support_and_acceptance_corruption(self):
        def joint(v): v["parents"][0]["joint"][0][0]["sum_b2"] += 1
        def source(v): v["parents"][0]["sources"]["prompt"]["totals"]["count"] += 1
        def tail(v): v["parents"][0]["tails"]["C"][0]["selected_count"] += 1
        def support(v): v["parents"][0]["support"]["ood_count"] += 1
        def owner(v): v["canonical_owner_equivalence"] = "PASS"
        def fp(v): v["snapshot"]["scientific_fingerprints"]["f4"]["body"] = "a"*64
        def row(v): v["rows"] = []
        def path(v): v["provenance"]["input_basenames"]["f3"] = "/private/map.json"
        for change in (joint, source, tail, support, owner, fp, row, path):
            altered = copy.deepcopy(self.diagnostic); altered.pop("fingerprint")
            change(altered)
            with self.subTest(change=change.__name__), self.assertRaises(debug.ContributionDebugError):
                debug.validate_diagnostic(debug.seal(altered))

    def test_zero_negative_large_weights_signed_cancellation_and_squared_moments(self):
        a = debug.aggregate
        s = a.summary([0, 2, 2, 1e8], [3, 2, 2, 1], [0, 2, -2, 1e8])
        a.validate_summary(s)
        self.assertEqual(s["zero_weight_count"], 1)
        self.assertEqual(s["sum_b2"], math.fsum([0, 4, 4, 1e16]))
        self.assertEqual(s["sum_b_rw2"], math.fsum([0, 16, 16, 1e16]))
        self.assertEqual(s["B"], 1e8)
        self.assertEqual(s["A"], 1e8+4)
        cancel = a.summary([1, 1], [1, 1], [1, -1])
        self.assertFalse(cancel["signed_relative_delta"]["available"])
        self.assertEqual(cancel["neff"]["abs_b"]["value"], 2)
        for w, c in (([-1], [1]), ([1], [0]), ([1e308], [2])):
            with self.assertRaises((ValueError, OverflowError)): a.summary(w, c, w)

    def test_fixed_mm_flows_geometry_and_all_source_joint_closure(self):
        a = debug.aggregate
        for value, expected in ((0.7, 0), (0.8, 1), (1.4, 40), (1.5, 41)):
            self.assertEqual(a.bin_index(value, a.MM_EDGES), expected)
        self.assertEqual(len(self.parents), 3)
        for p in self.parents:
            self.assertTrue(debug.validate_parent(p))
            self.assertEqual(sum(cell["count"] for row in p["joint"] for cell in row), p["totals"]["count"])

    def test_coherent_zero_weight_mm_flows_and_ood_use_frozen_factors(self):
        inputs = copy.deepcopy(self.f.inputs)
        application = inputs[0]["contract"]["application_records"]
        application[0]["baseline_pion_weight_w0"] = 0.0
        application[0]["signed_baseline_event_contribution"] = 0.0
        application[0]["analysis_MM"] = 0.79
        application[1]["analysis_MM"] = 1.41
        outside = next(r for r in application if r["source_label"] != "prompt")
        outside["P_hgcer_xAtCer"] = 1e5
        for artifact in inputs:
            fixtures.fixtures.f1_fixtures._seal_f1_artifact(artifact)
        hashes = {s: fixtures.digest(a) for s, a in zip(debug.SETTINGS, inputs)}
        map_value = fixtures.fixtures.f4_fixtures._f3(inputs, hashes)
        map_sha = fixtures.digest(map_value)
        stored = fixtures.fixtures.f4_fixtures._authority(map_value, map_sha)
        correction = debug.f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(
            inputs, map_value, f1_input_file_hashes=hashes, f3_input_file_sha256=map_sha,
            accepted_f3_runtime_authority_by_kinematic=stored, input_paths={"f1": {}, "f3": "synthetic"})
        self.assertTrue(correction["correction"]["available"])
        snapshot, parents = debug.reproduce_snapshot(inputs, map_value, correction,
            hashes={**hashes, "f3": map_sha, "f4": fixtures.digest(correction)})
        self.assertNotEqual(snapshot["fingerprint"], self.snapshot["fingerprint"])
        self.assertEqual(sum(p["totals"]["zero_weight_count"] for p in parents), 1)
        self.assertEqual(sum(p["mm"][0]["count"] for p in parents), 1)
        self.assertEqual(sum(p["mm"][-1]["count"] for p in parents), 1)
        self.assertGreater(sum(p["support"]["ood_count"] for p in parents), 0)
        for parent in parents:
            self.assertEqual(parent["support"]["ood_final_correction"], 1/parent["f4_parent"]["parent_normalization"])

    def test_preflight_pass_is_canonical_no_output_or_renderer_or_git_claim(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); locator, _ = self.write_locator(root)
            before = {p.name: p.read_bytes() for p in root.iterdir()}
            with mock.patch.object(cli, "write_review_pdf") as render, redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(cli.main(["--inputs-manifest", str(locator), "--preflight-only"]), 0)
            render.assert_not_called()
            receipt = json.loads(stdout.getvalue())
            self.assertEqual(receipt["snapshot_fingerprint"], self.snapshot["fingerprint"])
            self.assertEqual(receipt["preflight"], "PASS")
            self.assertFalse(receipt["outputs_written"])
            self.assertEqual(set(receipt["observations"]), set(debug.ROLES))
            self.assertTrue(all(r["fingerprint_validated"] for r in receipt["observations"].values()))
            self.assertNotIn(str(root), stdout.getvalue())
            self.assertEqual(before, {p.name: p.read_bytes() for p in root.iterdir()})

    def test_preflight_failure_retains_all_observed_raw_hashes_no_figures(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); locator, _ = self.write_locator(root)
            with (root/self.names["Center-highe"]).open("ab") as f: f.write(b"\n")
            with mock.patch.object(cli, "write_review_pdf") as render, redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(cli.main(["--inputs-manifest", str(locator), "--preflight-only"]), 1)
            receipt = json.loads(stdout.getvalue())
            self.assertIn("f3_fingerprint_input_content_mismatch", receipt["first_failure"])
            self.assertEqual(set(receipt["observations"]), set(debug.ROLES))
            self.assertIsNone(receipt["snapshot_fingerprint"])
            self.assertFalse(receipt["outputs_written"])
            render.assert_not_called()

    def test_locator_missing_roles_missing_file_duplicate_and_swapped_setting(self):
        for mode in ("roles", "missing", "duplicate", "swap", "extra"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); locator, manifest = self.write_locator(root)
                if mode == "roles": manifest["f1"].pop("Right-highe")
                elif mode == "missing": manifest["f3"] = "absent.json"
                elif mode == "duplicate": manifest["f4"] = manifest["f3"]
                elif mode == "swap": manifest["f1"]["Left-lowe"], manifest["f1"]["Left-highe"] = manifest["f1"]["Left-highe"], manifest["f1"]["Left-lowe"]
                else: manifest["accepted_authority_match"] = True
                locator.write_text(json.dumps(manifest))
                receipt, data = cli.preflight(locator)
                self.assertIsNone(data); self.assertEqual(receipt["preflight"], "FAIL")
                self.assertFalse(receipt["outputs_written"])

    def test_duplicate_json_and_nonfinite_rejected(self):
        for raw in (b'{"f1":{},"f1":{}}', b'{"x":NaN}', b'[]'):
            with self.assertRaises(ValueError): cli.parse_json(raw)

    def test_expected_fingerprint_required_and_replaced_file_cannot_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); locator = self.prepare(root)
            with mock.patch.object(cli, "write_review_pdf") as render, redirect_stderr(io.StringIO()):
                self.assertEqual(cli.main(["--inputs-manifest", str(locator)]), 1)
                args = self.args(root, locator); args[3] = "a"*64
                self.assertEqual(cli.main(args), 1)
                with (root/self.names["Left-lowe"]).open("ab") as f: f.write(b" ")
                self.assertEqual(cli.main(self.args(root, locator)), 1)
            render.assert_not_called()
            self.assertEqual(list((root/"artifacts").iterdir()), [])
            self.assertEqual(list((root/"bundles").iterdir()), [])

    def test_change_during_preflight_or_immediately_before_production(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); locator = self.prepare(root)
            shared = debug.reproduce_snapshot
            def change(*args, **kwargs):
                result = shared(*args, **kwargs)
                with (root/self.names["f4"]).open("ab") as f: f.write(b"\n")
                return result
            with mock.patch.object(debug, "reproduce_snapshot", side_effect=change):
                receipt, data = cli.preflight(locator)
            self.assertIsNone(data); self.assertIn("snapshot_changed_during_preflight:f4", receipt["first_failure"])
            locator, _ = self.write_locator(root)
            provenance = cli.git_provenance
            def change_before(names):
                with (root/self.names["f3"]).open("ab") as f: f.write(b"\n")
                return provenance(names)
            with mock.patch.object(cli, "git_provenance", side_effect=change_before), redirect_stderr(io.StringIO()) as stderr:
                self.assertEqual(cli.main(self.args(root, locator)), 1)
            self.assertIn("snapshot_changed_before_production:f3", stderr.getvalue())
            self.assertEqual(list((root/"artifacts").iterdir()), [])

    def mock_pdf(self, path, persisted):
        # Validate serialized aggregates independently; mock no scientific gate.
        debug.validate_diagnostic(cli.parse_json(Path(persisted).read_bytes()))
        raw = b"%PDF-1.4\n" + b"/Type /Page\n"*8 + b"/Type /Pages /Count 8\n%%EOF\n"
        path.write_bytes(raw)
        return {"pages": 8, "page_titles": list(cli.PAGE_TITLES), "renderer_input": "persisted aggregates only"}

    def test_complete_owner_mock_pdf_verified_zip_determinism_and_tamper(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); locator = self.prepare(root)
            before = {p.name: p.read_bytes() for p in root.glob("*.json")}
            with mock.patch.object(cli, "write_review_pdf", side_effect=self.mock_pdf), redirect_stdout(io.StringIO()) as stdout:
                self.assertEqual(cli.main(self.args(root, locator)), 0)
            path = Path(stdout.getvalue().strip()); manifest = cli.verify_zip(path)
            self.assertEqual(manifest["snapshot_fingerprint"], self.snapshot["fingerprint"])
            self.assertEqual(manifest["pdf"]["pages"], 8)
            self.assertEqual(before, {p.name: p.read_bytes() for p in root.glob("*.json")})
            self.assertEqual(list((root/"artifacts").iterdir()), [])
            with zipfile.ZipFile(path) as archive:
                members = {name: archive.read(name) for name in cli.MEMBERS}
                self.assertNotIn(str(root).encode(), members["diagnostic.json"])
            repeat = root/"repeat.zip"; cli.create_zip(repeat, members)
            self.assertEqual(path.read_bytes(), repeat.read_bytes())
            bad = root/"bad.zip"; cli.create_zip(bad, {**members, "diagnostic.json": b"{}"})
            with self.assertRaisesRegex(ValueError, "digest"): cli.verify_zip(bad)

    def test_persisted_only_renderer_uses_candidate_labels(self):
        import matplotlib
        matplotlib.use("Agg", force=True)
        from matplotlib.backends.backend_pdf import PdfPages
        import numpy as np
        figures = []
        original_save = PdfPages.savefig
        def capture(pdf, figure, *args, **kwargs):
            figures.append(figure)
            return original_save(pdf, figure, *args, **kwargs)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); js = root/"diagnostic.json"; pdf = root/"diagnostic.pdf"
            js.write_bytes(cli.json_bytes(self.diagnostic))
            before = js.read_bytes()
            with mock.patch.object(debug, "reproduce_snapshot", side_effect=AssertionError("no transient inputs")), \
                 mock.patch.object(cli, "load_inputs", side_effect=AssertionError("no input reload")), \
                 mock.patch.object(debug.f4, "build_pion_hgcer_method_a_parent_preserving_correction_with_review_data",
                                   side_effect=AssertionError("no factor rebuild")), \
                 mock.patch.object(PdfPages, "savefig", new=capture):
                receipt = cli.write_review_pdf(pdf, js)
            self.assertEqual(js.read_bytes(), before)
            self.assertEqual(cli.check_pdf(pdf.read_bytes()), 8)
            self.assertEqual(receipt, {"pages": 8, "page_titles": list(cli.PAGE_TITLES),
                                       "renderer_input": "persisted aggregates only"})
            self.assertEqual([fig._suptitle.get_text() for fig in figures],
                             [f"{i+1}. {title}" for i, title in enumerate(cli.PAGE_TITLES)])
            edges = self.diagnostic["display_geometry"]["mm_edges"]
            centers, limits, flows = cli.mm_display_geometry(edges)
            self.assertEqual(len(centers), 40)
            self.assertAlmostEqual(centers[0], .8075)
            self.assertAlmostEqual(centers[-1], 1.3925)
            for t, parent in enumerate(self.diagnostic["parents"]):
                axis = figures[4].axes[2*t+1]
                for j, field in enumerate(("B", "B_rw", "signed_delta")):
                    normal, flow = axis.lines[2*j:2*j+2]
                    np.testing.assert_array_equal(normal.get_xdata(), centers)
                    np.testing.assert_array_equal(normal.get_ydata(), [s[field] for s in parent["mm"][1:-1]])
                    np.testing.assert_array_equal(flow.get_xdata(), flows)
                    np.testing.assert_array_equal(flow.get_ydata(), [parent["mm"][0][field], parent["mm"][-1][field]])
                    self.assertEqual(flow.get_linestyle(), "None")
                    self.assertEqual(len(normal.get_ydata())+len(flow.get_ydata()), 42)
                self.assert_mm_axis(axis, edges, limits)
                for axis, field in zip(figures[5+t].axes[:4], ("signed_delta", "A_rw", "neff", "count")):
                    image, = axis.images
                    expected = [[s["neff"]["abs_b_rw"]["value"] if s["neff"]["abs_b_rw"]["available"] else np.nan
                                 for s in row] for row in parent["joint"]] if field == "neff" else [[s[field] for s in row] for row in parent["joint"]]
                    self.assertEqual(image.get_array().shape, (9, 42))
                    np.testing.assert_array_equal(image.get_array().filled(np.nan), expected)
                    np.testing.assert_allclose(image.get_extent(), (*limits, -.5, 8.5))
                    # Every interior pixel boundary maps to its persisted physical edge.
                    np.testing.assert_allclose(np.linspace(*limits, 43)[1:-1], edges, atol=1e-15)
                    self.assertEqual(image.get_interpolation(), "nearest")
                    self.assertEqual(image.origin, "lower")
                    if field == "signed_delta":
                        limit = max(float(np.max(np.abs(expected))), 1e-30)
                        self.assertEqual(image.get_clim(), (-limit, limit))
                        self.assertEqual(image.get_cmap().name, "RdBu_r")
                    self.assert_mm_axis(axis, edges, limits)
            import inspect
            renderer = inspect.getsource(cli.write_review_pdf)
            self.assertIn("current-snapshot candidate — detached diagnostic", renderer)
            self.assertNotIn("October 4", renderer)

    def assert_mm_axis(self, axis, edges, limits):
        import numpy as np
        np.testing.assert_allclose(axis.get_xlim(), limits)
        np.testing.assert_allclose(axis.get_xticks(), [.8, .9, 1., 1.1, 1.2, 1.3, 1.4])
        self.assertEqual([tick.get_text() for tick in axis.get_xticklabels()],
                         ["0.8", "0.9", "1.0", "1.1", "1.2", "1.3", "1.4"])
        self.assertEqual([text.get_text() for text in axis.texts], ["UF", "OF"])
        self.assertIn("Analysis missing mass [GeV]", axis.get_xlabel())
        self.assertIn("outside-range categories", axis.get_xlabel())
        np.testing.assert_array_equal(axis.lines[-2].get_xdata(), [edges[0]]*2)
        np.testing.assert_array_equal(axis.lines[-1].get_xdata(), [edges[-1]]*2)

    def test_failed_pdf_or_final_zip_gate_cleans_only_owned_outputs(self):
        for mode in ("pdf", "final"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); locator = self.prepare(root)
                keep = root/"artifacts"/"keep.txt"; keep.write_text("keep")
                original_verify = cli.verify_zip
                def verify(path, expected=None):
                    if mode == "final" and path.parent == root/"bundles": raise ValueError("final_tamper")
                    return original_verify(path, expected)
                def render(path, persisted):
                    receipt = self.mock_pdf(path, persisted)
                    if mode == "pdf": path.write_bytes(b"%PDF-1.4\n/Type /Page\n/Type /Pages /Count 1\n%%EOF")
                    return receipt
                with mock.patch.object(cli, "write_review_pdf", side_effect=render), mock.patch.object(cli, "verify_zip", side_effect=verify), redirect_stderr(io.StringIO()):
                    self.assertEqual(cli.main(self.args(root, locator)), 1)
                self.assertEqual(keep.read_text(), "keep")
                self.assertEqual(list((root/"artifacts").iterdir()), [keep])
                self.assertEqual(list((root/"bundles").iterdir()), [])

    def test_collision_and_checkout_output_are_rejected_without_touching_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); locator = self.prepare(root)
            collision = root/"bundles"/"kaonlt-current-snapshot-contribution-fixed.zip"; collision.write_text("keep")
            with mock.patch.object(cli.uuid, "uuid4", return_value=mock.Mock(hex="fixed")), mock.patch.object(cli, "preflight") as preflight, redirect_stderr(io.StringIO()):
                self.assertEqual(cli.main(self.args(root, locator)), 1)
            preflight.assert_not_called(); self.assertEqual(collision.read_text(), "keep")
            with self.assertRaises(ValueError): cli.roots(ROOT, root/"bundles")

    def test_source_known_basename_functions_and_no_input_search(self):
        self.assertEqual(f1.pion_hgcer_method_a_acceptance_contract_filename("Left", "kaon", debug.KINEMATIC, "lowe"),
                         "Left_kaon_pion-background_hgcer_method-a-acceptance-contract_Q4p4W2p74_lowe.json")
        self.assertTrue(debug.lineage.F6_3_CANDIDATE_BASENAMES_BY_KINEMATIC[debug.KINEMATIC]["f3"].endswith("current-baseline-candidate.json"))


if __name__ == "__main__":
    unittest.main()
