"""Local fixtures only: never invoke the analysis launcher or farm."""
from contextlib import contextmanager, ExitStack
from copy import deepcopy
import ast
import io
import json
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import Mock, patch
import zipfile

from testing import run_e8_2_left_lowe_scientific_audit_gate as owner

ROOT = Path(__file__).resolve().parents[1]
SHA = "a" * 40
MARKERS = ("Left / lowe debug analysis completed.\n"
           "Full high-epsilon processing is intentionally skipped in -d debug mode.\n")


def pages():
    return {"schema_version": "full_background_subtraction_page_manifest/v1",
        "pdf_basename": owner.PDF,
        "setting": {"kinematic_token": owner.KINEMATIC, "epsilon_setting": "low",
            "epsilon_filename_token": "lowe", "phi_setting": "Left", "particle_type": "kaon",
            "Q2": "4p4", "W": "2p74"}, "renderer_failures": [],
        "pages": [{"page_id": pid, "scope": f"t{n}", "t_index": n - 1,
            "semantic_stage": semantic,
            "represented_phi_inventory": [{"phi_index": i} for i in range(9)],
            "invalid_unavailable_children": [], "authoritative": False, "presentation_only": True}
            for pid, (n, semantic) in owner.PAGE_IDS.items()] +
            [{"page_id": "full_background.e8.handoff"}]}


def write_json(path, payload):
    path.write_text(json.dumps(payload), encoding="utf-8")


def artifacts(outdir):
    (outdir / owner.PDF).write_bytes(b"%PDF-1.4\nfixture\n")
    write_json(outdir / owner.PDF.replace(".pdf", "-manifest.json"), pages())
    write_json(outdir / f"kaon_FullAnalysis_{owner.KINEMATIC}_lowe.json",
        {"inpDict": {"ParticleType": "kaon", "EPSSET": "low", "Q2": "4p4", "W": "2p74",
                     "OutFilename": f"FullAnalysis_{owner.KINEMATIC}_lowe"},
         "histlist": [{"phi_setting": "Left"}]})
    base = outdir / f"kaon_FullAnalysis_{owner.KINEMATIC}_lowe_correction_ledger_no_empirical_residual"
    write_json(base.with_suffix(".json"), {"active_profile": "no_empirical_residual",
        "particle_type": "kaon", "epsset": "low", "q2": "4p4", "w": "2p74",
        "settings": [{"phi_setting": "Left"}]})
    base.with_suffix(".csv").write_text("row_kind,phi_setting\nsetting_total,Left\n", encoding="utf-8")


class OwnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.outdir = self.directory / "volatile/OUTPUT/Analysis/KaonLT"
        self.outdir.mkdir(parents=True)
        self.output = self.directory / "audit.zip"
        self.profile = owner.resolved_profile(ROOT, SHA)

    def test_profile_and_filtered_resolution(self):
        template = owner.collector.load_validation_profile(ROOT / owner.PROFILE)
        self.assertEqual(template["schema_version"], "pion_hgcer_validation_bundle_profile/v4")
        self.assertEqual(template["validation_profile"], owner.PROFILE_ID)
        self.assertEqual(template["settings"], owner.DECLARED_SETTINGS)
        self.assertEqual(template["artifacts"], owner.ARTIFACTS)
        self.assertEqual(owner.collector.resolve_settings("Left", "lowe", template), (("Left", "lowe"),))
        self.assertEqual(self.profile["source_identity"]["required_analysis_commit"], SHA)
        self.assertEqual(template["source_identity"]["required_analysis_commit"], owner.BASE_HEAD)

    def test_one_setting_profile_rejected_by_frozen_collector(self):
        payload = json.loads((ROOT / owner.PROFILE).read_text())
        payload["settings"] = owner.REQUESTED_SETTINGS
        path = self.directory / "single.json"
        write_json(path, payload)
        with self.assertRaisesRegex(ValueError, "validation_bundle_profile_invalid"):
            owner.collector.load_validation_profile(path)

    def test_profile_resolution_fail_closed(self):
        for sha in ("short", "A" * 40, "g" * 40):
            with self.subTest(sha=sha), self.assertRaisesRegex((ValueError, RuntimeError), "full_sha"):
                owner.resolved_profile(ROOT, sha)
        for field in ("validation_profile", "source_identity", "artifacts"):
            payload = json.loads((ROOT / owner.PROFILE).read_text())
            if field == "artifacts":
                payload[field]["settings"][0]["basename_template"] = "wrong.pdf"
            elif field == "source_identity":
                payload[field]["allowed_committed_files"] = ["src/example.py"]
            else:
                payload[field] = "wrong/v1"
            with self.subTest(field=field), patch.object(owner.collector, "load_validation_profile", return_value=payload):
                with self.assertRaises((ValueError, RuntimeError)):
                    owner.resolved_profile(ROOT, SHA)

    def test_exact_eighteen_pages(self):
        self.assertEqual(len(owner.PAGE_IDS), 18)
        self.assertEqual(owner.verify_pages(pages()), 19)

    def test_missing_and_duplicate_pages(self):
        for index in range(18):
            for duplicate in (False, True):
                payload = pages()
                if duplicate:
                    payload["pages"].insert(0, deepcopy(payload["pages"][index]))
                else:
                    del payload["pages"][index]
                with self.subTest(index=index, duplicate=duplicate), self.assertRaisesRegex((ValueError, RuntimeError), "missing_or_duplicate"):
                    owner.verify_pages(payload)

    def test_setting_and_manifest_identity(self):
        for key in pages()["setting"]:
            if key in ("Q2", "W"):
                continue
            payload = pages()
            payload["setting"][key] = "wrong"
            with self.subTest(key=key), self.assertRaisesRegex((ValueError, RuntimeError), "setting_invalid"):
                owner.verify_pages(payload)
        for key in ("schema_version", "pdf_basename"):
            payload = pages()
            payload[key] = "wrong"
            with self.subTest(key=key), self.assertRaises((ValueError, RuntimeError)):
                owner.verify_pages(payload)

    def test_parent_semantics_children_and_flags(self):
        changes = {"scope": "t0", "t_index": True, "semantic_stage": "wrong",
            "represented_phi_inventory": [{"phi_index": i} for i in reversed(range(9))],
            "invalid_unavailable_children": [0], "authoritative": True, "presentation_only": False}
        for key, value in changes.items():
            for index in range(18):
                payload = pages()
                payload["pages"][index][key] = value
                with self.subTest(key=key, index=index), self.assertRaises((ValueError, RuntimeError)):
                    owner.verify_pages(payload)
        for children in ([], [{"phi_index": i} for i in range(8)],
                         [{"phi_index": str(i)} for i in range(9)]):
            payload = pages()
            payload["pages"][0]["represented_phi_inventory"] = children
            with self.assertRaises((ValueError, RuntimeError)):
                owner.verify_pages(payload)

    def test_handoff_and_renderer_failures(self):
        for mode in ("missing", "early", "duplicate", "renderer"):
            payload = pages()
            if mode == "missing":
                payload["pages"].pop()
            elif mode == "early":
                payload["pages"].insert(0, payload["pages"].pop())
            elif mode == "duplicate":
                payload["pages"].append(payload["pages"][-1])
            else:
                payload["renderer_failures"] = [{"error": "render failed"}]
            with self.subTest(mode=mode), self.assertRaises((ValueError, RuntimeError)):
                owner.verify_pages(payload)

    def test_later_unavailable_pages_and_no_candidates_are_allowed(self):
        payload = pages()
        payload["pages"].insert(18, {"page_id": "full_background.e8_4.unavailable"})
        self.assertEqual(owner.verify_pages(payload), 20)
        artifacts(self.outdir)
        records, _ = owner.verify_artifacts(self.outdir, self.profile, 0, {})
        self.assertEqual(len(records), 5)
        self.assertFalse(any("f1" in name or "f4" in name for name in records))

    def test_missing_empty_stale_and_unchanged_artifacts(self):
        for name in owner.artifact_names(self.profile) - {owner.SUMMARY}:
            for mode in ("missing", "empty", "old", "unchanged"):
                artifacts(self.outdir)
                path = self.outdir / name
                before = {}
                started = 0
                if mode == "missing":
                    path.unlink()
                elif mode == "empty":
                    path.write_bytes(b"")
                elif mode == "old":
                    started = time.time_ns() + 10**9
                else:
                    before[name] = (path.stat().st_mtime_ns, path.stat().st_size, owner.sha256(path))
                with self.subTest(name=name, mode=mode), self.assertRaises((ValueError, RuntimeError)):
                    owner.verify_artifacts(self.outdir, self.profile, started, before)

    def test_invalid_pdf_and_strict_json(self):
        artifacts(self.outdir)
        (self.outdir / owner.PDF).write_bytes(b"wrong")
        with self.assertRaisesRegex((ValueError, RuntimeError), "pdf_invalid"):
            owner.verify_artifacts(self.outdir, self.profile, 0, {})
        artifacts(self.outdir)
        path = self.outdir / owner.PDF.replace(".pdf", "-manifest.json")
        path.write_text('{"x": NaN}')
        with self.assertRaises(ValueError):
            owner.verify_artifacts(self.outdir, self.profile, 0, {})

    def test_full_analysis_identity(self):
        name = f"kaon_FullAnalysis_{owner.KINEMATIC}_lowe.json"
        for key in ("ParticleType", "EPSSET", "Q2", "W", "OutFilename", "histlist"):
            artifacts(self.outdir)
            path = self.outdir / name
            payload = json.loads(path.read_text())
            if key == "histlist":
                payload[key].append({"phi_setting": "Center"})
            else:
                payload["inpDict"][key] = "wrong"
            write_json(path, payload)
            with self.subTest(key=key), self.assertRaisesRegex((ValueError, RuntimeError), "full_analysis"):
                owner.verify_artifacts(self.outdir, self.profile, 0, {})

    def test_ledger_identity_and_csv(self):
        name = f"kaon_FullAnalysis_{owner.KINEMATIC}_lowe_correction_ledger_no_empirical_residual"
        for key in ("active_profile", "particle_type", "epsset", "q2", "w", "settings", "csv"):
            artifacts(self.outdir)
            path = self.outdir / (name + ".json")
            payload = json.loads(path.read_text())
            if key == "csv":
                (self.outdir / (name + ".csv")).write_text("row_kind,phi_setting\nsetting_total,Right\n")
            else:
                payload[key] = [{"phi_setting": "Right"}] if key == "settings" else "wrong"
                write_json(path, payload)
            with self.subTest(key=key), self.assertRaisesRegex((ValueError, RuntimeError), "ledger"):
                owner.verify_artifacts(self.outdir, self.profile, 0, {})

    def collect(self):
        artifacts(self.outdir)
        write_json(self.outdir / owner.SUMMARY, {"fixture": True})
        profile_path = self.directory / "effective.json"
        write_json(profile_path, self.profile)
        with patch.object(owner.collector, "collect_source_state", return_value=("fixture state", SHA)), \
             patch.object(owner.collector, "collect_source_checks", return_value=("fixture checks", [])), \
             patch.object(owner.collector, "_committed_identity", return_value=(True, [], [])):
            result = owner.collector.collect_validation_bundle(outdir=self.outdir,
                kinematic=owner.KINEMATIC, output=self.output, profile_path=profile_path,
                repo_root=ROOT, phi="Left", epsilon="lowe")
        self.assertEqual(result["returncode"], 0)
        records = {name: {"sha256": owner.sha256(self.outdir / name),
                   "bytes": (self.outdir / name).stat().st_size} for name in owner.artifact_names(self.profile)}
        return records

    def test_actual_generic_collector_only_packages_filtered_setting(self):
        records = self.collect()
        owner.verify_zip(self.output, SHA, records, self.profile)
        with zipfile.ZipFile(self.output) as archive:
            manifest = json.loads(archive.read("manifest.json"))
            self.assertEqual(manifest["requested_settings"], owner.REQUESTED_SETTINGS)
            self.assertEqual(len(manifest["settings"]), 1)
            self.assertEqual(len(records), 6)
            self.assertFalse(any(s in n for n in archive.namelist() for s in ("highe", "Center", "Right")))

    def test_zip_rejects_identity_inventory_hash_and_extra_settings(self):
        records = self.collect()
        with zipfile.ZipFile(self.output) as archive:
            original = {name: archive.read(name) for name in archive.namelist()}
        for mode in ("head", "required", "profile", "requested", "settings", "missing", "extra", "hash", "size"):
            content = dict(original)
            manifest = json.loads(content["manifest.json"])
            if mode in ("head", "required", "profile"):
                manifest[{"head": "git_head", "required": "required_analysis_commit", "profile": "validation_profile"}[mode]] = "wrong"
            elif mode == "requested":
                manifest["requested_settings"].append({"phi": "Center", "epsilon": "highe"})
            elif mode == "settings":
                manifest["settings"].append(deepcopy(manifest["settings"][0]))
            elif mode == "missing":
                del content["Left_lowe/" + owner.PDF]
            elif mode == "extra":
                content["Center_highe/unwanted.pdf"] = b"%PDF-fixture"
            else:
                entry = manifest["settings"][0]["artifacts"]["procedure_pdf"]
                entry["sha256" if mode == "hash" else "byte_size"] = "wrong" if mode == "hash" else 999
            content["manifest.json"] = json.dumps(manifest).encode()
            path = self.directory / (mode + ".zip")
            with zipfile.ZipFile(path, "w") as archive:
                for name, raw in content.items():
                    archive.writestr(name, raw)
            with self.subTest(mode=mode), self.assertRaises((ValueError, RuntimeError, KeyError)):
                owner.verify_zip(path, SHA, records, self.profile)

    def test_completion_markers_and_high_epsilon_rejection(self):
        path = self.directory / "child.log"
        for text, rc in ((MARKERS, 1), ("", 0), (MARKERS + "High Epsilon Completed!", 0)):
            path.write_text(text)
            with self.assertRaises((ValueError, RuntimeError)):
                owner.verify_completion(path, rc)
        path.write_text(MARKERS)
        owner.verify_completion(path, 0)

    def test_launcher_exactly_once_with_owned_cwd_and_environment(self):
        child = Mock(stdout=io.BytesIO(MARKERS.encode()))
        child.wait.return_value = 0
        worktree = self.directory / "owned-analysis"
        env = {"PYTHONPATH": "copied-ltsep"}
        with patch.object(owner.subprocess, "Popen", return_value=child) as launch:
            self.assertEqual(owner.run_analysis(worktree, self.directory / "launch.log", env), 0)
        launch.assert_called_once_with(["./run_Prod_Analysis.sh", "-d", "4p4", "2p74"],
            cwd=worktree, env=env, stdout=owner.subprocess.PIPE, stderr=owner.subprocess.STDOUT)

    @contextmanager
    def flow(self, failure=None):
        worktree = self.directory / "owned-analysis"
        worktree.mkdir(exist_ok=True)
        overlay = {"baseline_paths": {}, "fixture": "copied-ltsep"}
        probe = [{"paths": {"VOLATILEPATH": str(self.directory / "volatile"),
            "LTANAPATH": str(worktree), "ANATYPE": "Kaon",
            "OUTPATH": str(worktree / "OUTPUT/Analysis/KaonLT")}}]
        events = []
        def runtime_overlay(*args):
            events.append("overlay")
            return {"fixture": "overlay"}, overlay
        def path_probes(*args):
            events.append("probe")
            return probe
        def external_preflight(*args):
            events.append("external")
            self.assertTrue((worktree / "OUTPUT").is_symlink())
            return {"passed": True}
        @contextmanager
        def owned(repo, sha, kind):
            events.append("enter:" + kind)
            try:
                yield worktree
            finally:
                events.append("cleanup:" + kind)
                if kind == "analysis" and (worktree / "OUTPUT").is_symlink():
                    (worktree / "OUTPUT").unlink()
                if failure == "cleanup" and kind == "analysis":
                    raise RuntimeError("owned_cleanup_failed")
        def run(cwd, log, env):
            self.assertEqual(cwd, worktree)
            self.assertNotEqual(cwd, ROOT)
            self.assertEqual(env, {"fixture": "overlay"})
            self.assertTrue((cwd / "OUTPUT").is_symlink())
            receipt = json.loads(owner.isolation.accepted.gate_status_path(self.outdir, self.output).read_text())
            self.assertIs(receipt["debug_output_link"]["passed"], True)
            events.append("analysis")
            artifacts(self.outdir)
            # Explicit fresh fixture timestamps avoid filesystem clock rounding.
            fixture_ns = time.time_ns() + 10**9
            for name in owner.artifact_names(self.profile) - {owner.SUMMARY}:
                os.utime(self.outdir / name, ns=(fixture_ns, fixture_ns))
            log.write_text(MARKERS)
            return 1 if failure == "child" else 0
        def preserve(*args):
            events.append("preserve")
            if failure == "checkout":
                raise RuntimeError("ordinary_checkout_changed")
            return {"passed": True}
        def ltsep(*args):
            events.append("ltsep")
            if failure == "ltsep":
                raise RuntimeError("installed_ltsep_changed")
            return {"passed": True}
        with ExitStack() as stack:
            original_output_link = owner.prepare_debug_output_link
            def output_link(*args):
                events.append("output-link")
                if failure == "link":
                    raise ValueError("debug_output_path_already_exists")
                return original_output_link(*args)
            for name, value in (("preflight", Mock(return_value={"head": SHA})),
                ("owned_worktree", owned), ("prepare_runtime_overlay", Mock(side_effect=runtime_overlay)),
                ("probe_paths", Mock(side_effect=path_probes)),
                ("prepare_debug_output_link", Mock(side_effect=output_link)),
                ("external_symlink_preflight", Mock(side_effect=external_preflight)),
                ("collector_source_preflight", Mock(side_effect=RuntimeError("source_check_failed")) if failure == "source" else Mock(return_value=[])),
                ("run_analysis", Mock(side_effect=run)), ("verify_preservation", Mock(side_effect=preserve)),
                ("verify_ltsep_preservation", Mock(side_effect=ltsep))):
                stack.enter_context(patch.object(owner, name, value))
            stack.enter_context(patch.object(owner, "collection_module", lambda source: self.module()))
            stack.enter_context(patch.object(owner.collector, "collect_source_state", return_value=("state", SHA)))
            stack.enter_context(patch.object(owner.collector, "collect_source_checks", return_value=("checks", [])))
            stack.enter_context(patch.object(owner.collector, "_committed_identity", return_value=(True, [], [])))
            collection = stack.enter_context(patch.object(owner.collector, "collect_validation_bundle", wraps=owner.collector.collect_validation_bundle))
            yield events, collection

    @contextmanager
    def module(self):
        yield owner.collector

    def test_owner_isolated_flow_filtered_collection_and_companions(self):
        with self.flow() as (events, collection):
            self.assertEqual(owner.execute_gate(ROOT, self.outdir, self.output, SHA), self.output)
            owner.run_analysis.assert_called_once()
            owner.prepare_debug_output_link.assert_called_once()
        self.assertEqual([event for event in events if event in
            ("overlay", "probe", "output-link", "external", "analysis")],
            ["overlay", "probe", "output-link", "external", "analysis"])
        self.assertLess(events.index("cleanup:analysis"), events.index("preserve"))
        self.assertLess(events.index("ltsep"), events.index("enter:source"))
        collection.assert_called_once()
        self.assertEqual(collection.call_args.kwargs["phi"], "Left")
        self.assertEqual(collection.call_args.kwargs["epsilon"], "lowe")
        status = json.loads(owner.isolation.accepted.gate_status_path(self.outdir, self.output).read_text())
        self.assertEqual(status["status"], "success")
        for flag in ("analysis_started", "analysis_completed", "artifact_verification_completed",
                     "collector_source_preflight_completed", "collection_completed", "zip_verification_completed"):
            self.assertIs(status[flag], True)
        for suffix in (".log", "-run-summary.json", "-gate-status.json"):
            self.assertTrue(self.output.with_name(self.output.stem + suffix).is_file())
        summary = json.loads((self.outdir / owner.SUMMARY).read_text())
        self.assertEqual(summary["boundaries"], owner.BOUNDARIES)
        self.assertEqual(summary["debug_output_link"], status["debug_output_link"])
        self.assertIs(summary["debug_output_link"]["passed"], True)

    def output_paths(self):
        worktree = self.directory / "fake-worktree"
        worktree.mkdir(exist_ok=True)
        return worktree, {"LTANAPATH": str(worktree),
            "VOLATILEPATH": str(self.directory / "volatile"), "ANATYPE": "Kaon",
            "OUTPATH": str(worktree / "OUTPUT/Analysis/KaonLT")}

    def test_debug_output_link_and_launcher_mkdir_keep_external_input_visible(self):
        worktree, paths = self.output_paths()
        fixture = self.outdir / "Center_kaon_Q4p4W2p74_input.fixture"
        fixture.write_bytes(b"filesystem visibility fixture only")
        record = owner.prepare_debug_output_link(worktree, paths, self.outdir)
        link, child = worktree / "OUTPUT", Path(paths["OUTPATH"])
        self.assertTrue(link.is_symlink())
        self.assertEqual(os.readlink(link), str(self.directory / "volatile/OUTPUT"))
        self.assertEqual(child.resolve(), self.outdir)
        self.assertEqual((child / fixture.name).read_bytes(), fixture.read_bytes())
        child.mkdir(parents=True, exist_ok=True)  # launcher's mkdir -p
        self.assertTrue(link.is_symlink())
        self.assertEqual(link.resolve(), (self.directory / "volatile/OUTPUT").resolve())
        self.assertEqual(child.resolve(), self.outdir)
        self.assertEqual((child / fixture.name).read_bytes(), fixture.read_bytes())
        self.assertEqual(record, {"link_path": str(link),
            "literal_target": str(self.directory / "volatile/OUTPUT"),
            "resolved_target": str((self.directory / "volatile/OUTPUT").resolve()),
            "resolved_child_outpath": str(self.outdir), "passed": True})

    def test_debug_output_link_never_replaces_preexisting_paths(self):
        worktree, paths = self.output_paths()
        for kind in ("directory", "file", "wrong", "dangling", "correct"):
            with self.subTest(kind=kind):
                link = worktree / "OUTPUT"
                if kind == "directory":
                    link.mkdir()
                elif kind == "file":
                    link.write_bytes(b"untouched")
                else:
                    target = self.directory / ("missing" if kind == "dangling" else
                             "volatile/OUTPUT" if kind == "correct" else "volatile")
                    link.symlink_to(target, target_is_directory=True)
                before = os.readlink(link) if link.is_symlink() else None
                with self.assertRaisesRegex(ValueError, "path_already_exists"):
                    owner.prepare_debug_output_link(worktree, paths, self.outdir)
                self.assertTrue(os.path.lexists(link))
                if before is not None:
                    self.assertTrue(link.is_symlink())
                    self.assertEqual(os.readlink(link), before)
                    link.unlink()  # remove this test's fixture only
                elif kind == "directory":
                    self.assertTrue(link.is_dir())
                    link.rmdir()
                else:
                    self.assertEqual(link.read_bytes(), b"untouched")
                    link.unlink()

    def test_debug_output_link_rejects_unvalidated_path_authority(self):
        worktree, paths = self.output_paths()
        cases = [("LTANAPATH", str(self.directory), "ltanapath_mismatch"),
            ("OUTPATH", str(self.outdir), "lexical_outpath_mismatch"),
            ("VOLATILEPATH", "relative", "not_absolute"),
            ("VOLATILEPATH", "", "not_absolute"),
            ("VOLATILEPATH", str(self.directory / "missing"), "external_output_missing"),
            ("ANATYPE", "../Kaon", "anatype_invalid")]
        for key, value, reason in cases:
            bad = {**paths, key: value}
            if key == "VOLATILEPATH" and value == "":
                del bad[key]
            with self.subTest(key=key, value=value), self.assertRaisesRegex(ValueError, reason):
                owner.prepare_debug_output_link(worktree, bad, self.outdir)
            self.assertFalse(os.path.lexists(worktree / "OUTPUT"))
        with self.assertRaisesRegex(ValueError, "artifact_root_mismatch"):
            owner.prepare_debug_output_link(worktree, paths, self.directory)
        with self.assertRaisesRegex(ValueError, "worktree_invalid"):
            owner.prepare_debug_output_link(Path("relative"), paths, self.outdir)
        self.assertFalse(os.path.lexists(worktree / "OUTPUT"))

    def test_child_failure_checks_preservation_and_does_not_collect(self):
        with self.flow("child") as (events, collection), self.assertRaisesRegex((ValueError, RuntimeError), "analysis_failed"):
            owner.execute_gate(ROOT, self.outdir, self.output, SHA)
        self.assertIn("cleanup:analysis", events)
        self.assertIn("preserve", events)
        self.assertIn("ltsep", events)
        collection.assert_not_called()
        self.assertFalse(self.output.exists())
        self.assertEqual(json.loads(owner.isolation.accepted.gate_status_path(self.outdir, self.output).read_text())["status"], "failed")

    def test_preservation_and_cleanup_fail_closed(self):
        for failure in ("checkout", "ltsep", "cleanup"):
            # Each attempt owns a distinct status/log name.
            self.output = self.directory / (failure + ".zip")
            with self.subTest(failure=failure), self.flow(failure) as (events, collection), self.assertRaises((ValueError, RuntimeError)):
                owner.execute_gate(ROOT, self.outdir, self.output, SHA)
            self.assertIn("ltsep", events)
            collection.assert_not_called()
            self.assertFalse(self.output.exists())

    def test_source_failure_precedes_analysis(self):
        with self.flow("source") as (events, collection), self.assertRaisesRegex((ValueError, RuntimeError), "source_check_failed"):
            owner.execute_gate(ROOT, self.outdir, self.output, SHA)
        self.assertNotIn("analysis", events)
        collection.assert_not_called()

    def test_output_link_failure_precedes_child_and_preserves_checkout(self):
        with self.flow("link") as (events, collection), self.assertRaisesRegex(ValueError, "path_already_exists"):
            owner.execute_gate(ROOT, self.outdir, self.output, SHA)
        self.assertNotIn("external", events)
        self.assertNotIn("analysis", events)
        self.assertIn("preserve", events)
        self.assertIn("ltsep", events)
        collection.assert_not_called()
        receipt = json.loads(owner.isolation.accepted.gate_status_path(self.outdir, self.output).read_text())
        self.assertIs(receipt["analysis_started"], False)
        self.assertEqual(receipt["status"], "failed")

    def test_preflight_identity_and_dirty_source(self):
        good = {"branch": "test", "head": SHA, "origin_test": SHA, "porcelain": ""}
        with patch.object(owner.isolation, "snapshot", return_value=good), patch.object(owner, "git"):
            self.assertEqual(owner.preflight(ROOT, SHA), good)
        for key, value in (("branch", "main"), ("head", "b" * 40), ("origin_test", "b" * 40),
                           ("porcelain", " M src/utility/background_config.py")):
            state = {**good, key: value}
            with self.subTest(key=key), patch.object(owner.isolation, "snapshot", return_value=state), \
                 patch.object(owner, "git"), self.assertRaises((ValueError, RuntimeError)):
                owner.preflight(ROOT, SHA)
        state = {**good, "head": owner.BASE_HEAD, "origin_test": owner.BASE_HEAD}
        with patch.object(owner.isolation, "snapshot", return_value=state), self.assertRaisesRegex((ValueError, RuntimeError), "not_committed"):
            owner.preflight(ROOT, owner.BASE_HEAD)

    def test_no_scientific_materialization_or_lineage_calls(self):
        tree = ast.parse(Path(owner.__file__).read_text())
        calls = {ast.unparse(node.func) for node in ast.walk(tree) if isinstance(node, ast.Call)}
        self.assertFalse(any(token in name for name in calls for token in
            ("materialize", "candidate", "f6_3", "stage_candidates", "verify_f4_reproduction")))
        self.assertEqual(sum(isinstance(node, ast.Call) and ast.unparse(node.func) == "run_analysis"
                             for node in ast.walk(tree)), 1)


if __name__ == "__main__":
    unittest.main()
