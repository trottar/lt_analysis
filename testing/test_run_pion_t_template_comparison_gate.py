"""Synthetic standard-library fixtures; no production launch, ROOT or network."""
from contextlib import contextmanager, ExitStack, redirect_stderr
import ast
from copy import deepcopy
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import unittest
from unittest.mock import Mock, patch
import zipfile

from testing import run_pion_t_template_comparison_gate as owner

ROOT = Path(__file__).resolve().parents[1]
SHA = "a" * 40
FP = "b" * 64
MARKERS = ("All available low-epsilon settings completed in debug mode.\n"
           "Full high-epsilon processing is intentionally skipped in -d debug mode.\n")
IDENTITIES = {name: {"sha256": FP, "bytes": 12} for name in owner.SOURCE_FILES}


def histogram():
    return {"status": "available", "stored_sumw2": [0., 1., 0.],
            "content_error_fingerprint": FP}


def payload(setting, pdf):
    parents = []
    edges = [.4, .55, .7, .9]
    for ti in range(3):
        fit = {"H_pion_control_input": histogram(), "H_kaon_nosub_input": histogram(),
               "pion_component_alignment": {"components": {c: {} for c in owner.COMPONENTS}},
               "resolved_subtraction_config": {"profile": "no_empirical_residual"}}
        parents.append({"t_bin_index": ti, "t_edges": edges[ti:ti+2],
            "pion_parent_id": setting + str(ti), "canonical_interval_pair_id": "pair",
            "canonical_interval_pair_hash": FP, "coordinate_fingerprint": FP,
            "baseline_parent_fingerprint": FP, "authoritative_parent_unchanged": True,
            "production_applied": False, "scientifically_accepted": False,
            "alignment_mode": "fixed_parent_alignment_template_shape_test",
            "baseline_A": deepcopy(fit), "candidate_B": deepcopy(fit),
            "frozen_inputs": {"pion": histogram(), "kaon": histogram()},
            "status": "available", "rejection_reasons": [],
            "templates": {c: {"status": "genuinely_local", "normalization_count": 1,
                "normalized_integral": 1., "canonical_selection_closure": True,
                "weighted_child_to_parent_closure": True, "production_global_selection_closure": True,
                "t_edges": edges[ti:ti+2], "provenance": {"status": "available", "sha256": FP,
                    "path": "synthetic.root", "tree_name": "h666",
                    "weight_expression": "getattr(evt, 'iter_weight', 1.0)",
                    "coordinate": "missmass; -t; phipq*180/pi"}}
                for c in owner.COMPONENTS}})
    return {"schema_version": "pion_t_template_comparison/v1", "diagnostic_only": True,
        "production_applied": False, "scientifically_accepted": False, "source_HEAD": SHA,
        "source_identity_rejection_reasons": [], "source_files": {
            name: {"status": "available", **record} for name, record in IDENTITIES.items()},
        "phi_setting": setting, "setting": {"Q2": "4p4", "W": "2p74", "EPSSET": "low", "ParticleType": "kaon"},
        "enable_environment_variable": owner.ENABLE_ENV,
        "pdf": {"status": "available", "path": str(pdf), "sha256": owner.sha256(pdf), "bytes": pdf.stat().st_size},
        "canonical_t_edges": edges, "parents": parents, "status": "available"}


def artifacts(outdir):
    for setting in ("Center", "Left"):
        stem = f"{setting}_kaon_Q4p4W2p74_lowe_pion-t-template-comparison"
        pdf = outdir / (stem + ".pdf")
        pdf.write_bytes(b"%PDF-1.4\nsynthetic fixture\n")
        (outdir / (stem + ".json")).write_text(json.dumps(payload(setting, pdf)), encoding="utf-8")


class GateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name).resolve()
        self.repo = self.directory / "repo"
        self.repo.mkdir()
        self.outdir = self.directory / "volatile/OUTPUT/Analysis/KaonLT"
        self.outdir.mkdir(parents=True)
        self.transfer = self.directory / "transfer"
        self.transfer.mkdir()
        self.output = self.transfer / "comparison.zip"
        self.profile = owner.resolved_profile(ROOT, SHA)
        self.before = dict.fromkeys(owner.names())

    def verify(self):
        with patch.object(owner, "pdf_pages", return_value=6):
            return owner.verify_artifacts(self.outdir, SHA, IDENTITIES, 0, self.before)

    def test_exact_profile_and_unchanged_collector_acceptance(self):
        tracked = owner.collector.load_validation_profile(ROOT / owner.PROFILE)
        self.assertEqual(tracked["artifacts"], owner.ARTIFACTS)
        self.assertEqual(tracked["settings"], owner.SETTINGS)
        self.assertEqual(tracked["artifacts"]["settings"], [])
        self.assertEqual(len(owner.names()), 4)
        self.assertEqual(tracked["source_identity"]["required_analysis_commit"], owner.REVIEWED_HEAD)
        self.assertEqual(self.profile["source_identity"]["required_analysis_commit"], SHA)
        self.assertTrue(all("{" not in name for name in owner.names()))
        self.assertEqual(len(owner.collector.resolve_settings(None, None, tracked)), 5)

    def test_profile_edits_rejected(self):
        for mutate in (lambda p: p["artifacts"]["settings"].append(p["artifacts"]["global"].pop()),
                       lambda p: p["source_identity"].update(required_analysis_commit=SHA),
                       lambda p: p.update(validation_profile="other")):
            changed = owner.collector.load_validation_profile(ROOT / owner.PROFILE)
            mutate(changed)
            with self.subTest(changed=changed), patch.object(owner.collector, "load_validation_profile", return_value=changed):
                with self.assertRaisesRegex(ValueError, "profile_contract_mismatch"):
                    owner.resolved_profile(ROOT, SHA)

    def test_source_identity_dirty_and_unexpected_successor_rejected(self):
        state = {"branch": "test", "head": SHA, "origin_test": SHA, "porcelain": ""}
        cases = [("branch", "main", "wrong_branch"), ("head", "c"*40, "wrong_head"),
                 ("origin_test", "c"*40, "source_not_observed_pushed"),
                 ("porcelain", " M src/main.py", "dirty_gate_source"),
                 ("porcelain", "?? injected.py", "dirty_gate_source")]
        for key, value, reason in cases:
            with self.subTest(key=key), patch.object(owner.isolation, "snapshot", return_value={**state, key: value}), patch.object(owner, "git", return_value=""):
                with self.assertRaisesRegex(ValueError, reason):
                    owner.preflight(self.repo, SHA)
        with patch.object(owner.isolation, "snapshot", return_value=state), patch.object(owner, "git", side_effect=["", "", "", "src/main.py"]):
            with self.assertRaisesRegex(ValueError, "unexpected_committed_source"):
                owner.preflight(self.repo, SHA)
        for sha in ("short", "A"*40, owner.BASE_HEAD):
            with patch.object(owner.isolation, "snapshot", return_value={**state, "head": sha, "origin_test": sha}):
                with self.assertRaises(ValueError):
                    owner.preflight(self.repo, sha)

    def test_reviewed_checkpoint_and_all_thirteen_future_changes(self):
        changed = owner.git(ROOT, "diff", "--name-only", owner.BASE_HEAD, owner.REVIEWED_HEAD).splitlines()
        allowed = {owner.PROFILE, "testing/run_pion_t_template_comparison_gate.py",
                   "testing/test_run_pion_t_template_comparison_gate.py"}
        extra = [p for p in changed if p not in allowed and not p.startswith("docs/memory/")]
        self.assertEqual(len(extra), 13)  # reproduce the original committed-range rejection
        state = {"branch": "test", "head": SHA, "origin_test": SHA, "porcelain": ""}
        with patch.object(owner.isolation, "snapshot", return_value=state), patch.object(owner, "git", return_value="") as run:
            self.assertEqual(owner.preflight(ROOT, SHA), state)
            self.assertIn((ROOT, "merge-base", "--is-ancestor", owner.BASE_HEAD, owner.REVIEWED_HEAD),
                          [call.args for call in run.call_args_list])
            self.assertEqual(run.call_args.args, (ROOT, "diff", "--name-only", owner.REVIEWED_HEAD, SHA))
        for path in extra + ["unrelated.py"]:
            with self.subTest(path=path), patch.object(owner.isolation, "snapshot", return_value=state), patch.object(
                    owner, "git", side_effect=["", "", "", path]):
                with self.assertRaisesRegex(ValueError, "unexpected_committed_source"):
                    owner.preflight(ROOT, SHA)

    def test_original_owner_reproduces_history_rejection_and_four_target_collision(self):
        raw = subprocess.check_output(["git", "show", owner.REVIEWED_HEAD + ":testing/run_pion_t_template_comparison_gate.py"], cwd=ROOT)
        old = {"__name__": "testing._original_ab_owner", "__package__": "testing", "__file__": owner.__file__}
        exec(compile(raw, "<reviewed-original-owner>", "exec"), old)
        state = {"branch": "test", "head": owner.REVIEWED_HEAD,
                 "origin_test": owner.REVIEWED_HEAD, "porcelain": ""}
        with patch.object(owner.isolation, "snapshot", return_value=state):
            with self.assertRaisesRegex(ValueError, "unexpected_committed_source"):
                old["preflight"](ROOT, owner.REVIEWED_HEAD)
        artifacts(self.outdir)
        before = owner.accepted_identities(self.outdir)
        old["preflight"] = lambda *args: state
        old["resolved_profile"] = lambda *args: self.profile
        with patch.object(owner.isolation.accepted, "REPO", self.repo), patch.object(
                owner.isolation.accepted, "ARTIFACTS", self.outdir), patch.object(
                owner.isolation.accepted, "GLOBUS", self.transfer), patch.object(
                owner.isolation, "verify_preservation", return_value={"passed": True}), patch.object(
                owner.debug, "run_analysis") as child:
            with self.assertRaisesRegex(ValueError, "comparison_artifact_preexists"):
                old["execute_gate"](self.repo, self.outdir, self.output, SHA)
        child.assert_not_called()
        self.assertEqual(before, owner.accepted_identities(self.outdir))

    def test_committed_source_bytes_and_missing_source(self):
        result = subprocess.CompletedProcess([], 0, stdout=b"committed bytes", stderr=b"")
        with patch.object(owner.subprocess, "run", return_value=result) as run:
            records = owner.committed_identities(self.repo, SHA)
        self.assertEqual(records[owner.SOURCE_FILES[0]]["sha256"], hashlib.sha256(result.stdout).hexdigest())
        self.assertEqual(run.call_count, len(owner.SOURCE_FILES))
        with patch.object(owner.subprocess, "run", return_value=subprocess.CompletedProcess([], 1, stdout=b"", stderr=b"missing")):
            with self.assertRaisesRegex(ValueError, "committed_source_missing"):
                owner.committed_identities(self.repo, SHA)

    def test_environment_scoped_and_conflicts_rejected(self):
        env = {"PATH": "existing", "PYTHONPATH": "overlay", "LT_ANALYSIS_DEBUG_LEFT_LOW": "1",
               "LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE": "old", "LT_ANALYSIS_ALLOW_UNPAIRED_CANONICAL_BINNING": "1"}
        original = dict(env)
        child = owner.child_environment(env)
        self.assertEqual(env, original)
        self.assertEqual(child, {"PATH": "existing", "PYTHONPATH": "overlay", owner.ENABLE_ENV: "1"})
        for extra in ({"LT_ANALYSIS_NEW_OPT_IN": "1"}, {owner.ENABLE_ENV: "0"}):
            with self.assertRaises(ValueError):
                owner.child_environment({**env, **extra})

    def test_single_mock_subprocess_debug_child(self):
        process = Mock(stdout=io.BytesIO(MARKERS.encode()))
        process.wait.return_value = 0
        env = owner.child_environment({})
        path = self.outdir / "child.log"
        with patch.object(owner.debug.subprocess, "Popen", return_value=process) as popen, patch("sys.stderr", io.StringIO()):
            result = owner.debug.run_analysis(self.repo, path, env)
        self.assertEqual(result, 0)
        popen.assert_called_once_with(owner.COMMAND, cwd=self.repo, env=env,
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        owner.verify_completion(path, result)

    def test_completion_failures(self):
        path = self.outdir / "child.log"
        for text, code in ((MARKERS, 1), ("Left / lowe debug analysis completed.", 0),
                           (MARKERS + "High Epsilon Completed!", 0), ("", 0)):
            path.write_text(text)
            with self.subTest(text=text, code=code), self.assertRaises(ValueError):
                owner.verify_completion(path, code)

    def test_six_page_artifacts_and_explicit_negative_evidence(self):
        artifacts(self.outdir)
        records, outcomes = self.verify()
        self.assertEqual(len(records), 4)
        self.assertEqual(outcomes["Center"]["status"], "available")
        path = next(self.outdir.glob("Center*.json"))
        p = json.loads(path.read_text())
        p["parents"][1].update(status="unavailable", candidate_B=None,
                                rejection_reasons=["pi_n:invalid_local_template_or_provenance"])
        p["parents"][1]["templates"]["pi_n"].update(status="unavailable", normalization_count=0)
        p["status"] = "partial_unavailable"
        path.write_text(json.dumps(p))
        _, outcomes = self.verify()
        self.assertEqual(outcomes["Center"]["status"], "partial_unavailable")
        self.assertTrue(outcomes["Center"]["rejection_reasons_by_t"]["1"])

    def test_structural_identity_authority_provenance_negative_cases(self):
        artifacts(self.outdir)
        path = next(self.outdir.glob("Center*.json"))
        good = json.loads(path.read_text())
        mutations = [
            lambda p: p.update(schema_version="wrong"), lambda p: p.update(production_applied=True),
            lambda p: p.update(source_HEAD="c"*40), lambda p: p.update(phi_setting="Left"),
            lambda p: p["setting"].update(EPSSET="high"), lambda p: p.update(status="unavailable"),
            lambda p: p["parents"].pop(), lambda p: p["parents"][1].update(t_bin_index=0),
            lambda p: p["parents"][0].update(t_edges=[.1,.2]),
            lambda p: p["parents"][0].update(authoritative_parent_unchanged=False),
            lambda p: p["parents"][0].update(alignment_mode="rescanned"),
            lambda p: p["parents"][0].update(baseline_A={}),
            lambda p: p["parents"][0]["candidate_B"]["H_pion_control_input"].update(content_error_fingerprint="c"*64),
            lambda p: p["parents"][0]["candidate_B"].update(resolved_subtraction_config={}),
            lambda p: p["parents"][0]["candidate_B"].update(pion_component_alignment={}),
            lambda p: p["parents"][0]["templates"]["pi_n"].update(normalization_count=2),
            lambda p: p["parents"][0]["templates"]["pi_n"].update(production_global_selection_closure=False),
            lambda p: p["parents"][0]["templates"]["pi_n"].update(provenance={}),
            lambda p: p["parents"][0]["templates"]["pi_n"]["provenance"].update(sha256="bad"),
            lambda p: p["source_files"][owner.SOURCE_FILES[0]].update(sha256="c"*64),
            lambda p: p["source_files"].pop(owner.SOURCE_FILES[0]),
            lambda p: p["pdf"].update(sha256="c"*64),
            lambda p: p["parents"][0].update(status="unavailable", candidate_B=None),
        ]
        for index, mutate in enumerate(mutations):
            changed = deepcopy(good)
            mutate(changed)
            path.write_text(json.dumps(changed))
            with self.subTest(index=index), self.assertRaises((ValueError, TypeError)):
                self.verify()

    def test_malformed_nonfinite_empty_missing_stale_and_pdf_failure(self):
        artifacts(self.outdir)
        path = next(self.outdir.glob("Center*.json"))
        good = path.read_bytes()
        for bad in (b"", b"{", b'{"x":NaN}', b'{"x":1e999}', b'{"x":1,"x":2}'):
            path.write_bytes(bad)
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                self.verify()
        path.write_bytes(good)
        with patch.object(owner, "pdf_pages", return_value=5), self.assertRaisesRegex(ValueError, "page_count"):
            owner.verify_artifacts(self.outdir, SHA, IDENTITIES, 0, self.before)
        pdf = next(self.outdir.glob("Center*.pdf"))
        pdf.write_bytes(b"not PDF")
        with self.assertRaisesRegex(ValueError, "signature"):
            self.verify()
        artifacts(self.outdir)
        with self.assertRaisesRegex(ValueError, "stale"):
            owner.verify_artifacts(self.outdir, SHA, IDENTITIES, time.time_ns()+10**9, self.before)
        path.unlink()
        with self.assertRaisesRegex(ValueError, "missing_or_empty"):
            self.verify()

    def test_pdf_backend_uses_structure_and_fails_closed(self):
        for stdout, code in (("Pages: 6\n", 0), ("Pages: 5\n", 0)):
            with patch.object(owner.subprocess, "run", return_value=subprocess.CompletedProcess([], code, stdout, "")):
                self.assertEqual(owner.pdf_pages(self.outdir / "x.pdf"), int(stdout.split()[1]))
        for result in (subprocess.CompletedProcess([], 1, "", "broken"),
                       subprocess.CompletedProcess([], 0, "no page structure", "")):
            with patch.object(owner.subprocess, "run", return_value=result), self.assertRaises(ValueError):
                owner.pdf_pages(self.outdir / "x.pdf")
        with patch.object(owner.subprocess, "run", side_effect=FileNotFoundError), self.assertRaises(FileNotFoundError):
            owner.pdf_pages(self.outdir / "x.pdf")

    def collect(self, output=None, profile=None, outdir=None):
        profile_path = self.directory / "effective.json"
        profile_path.write_text(json.dumps(profile or self.profile))
        module = owner.collector
        with patch.object(module, "collect_source_state", return_value=("source fixture", SHA)), patch.object(module, "collect_source_checks", return_value=("checks fixture", [])), patch.object(module, "_committed_identity", return_value=(True, [], [])):
            return module.collect_validation_bundle(outdir=outdir or self.outdir, kinematic="Q4p4W2p74",
                output=output or self.output, profile_path=profile_path, repo_root=self.repo)

    def test_global_only_real_collector_zip_and_hashes(self):
        artifacts(self.outdir)
        records, _ = self.verify()
        result = self.collect()
        self.assertEqual(result["returncode"], 0)
        owner.verify_zip(self.output, SHA, records)
        with zipfile.ZipFile(self.output) as archive:
            self.assertEqual(len(archive.namelist()), 13)
            self.assertEqual(sum(not name.endswith("/") and name.startswith("global/") for name in archive.namelist()), 4)
        self.assertEqual(len(result["manifest"]["settings"]), 5)

    def test_zip_inventory_manifest_and_hash_tampering(self):
        artifacts(self.outdir)
        records, _ = self.verify()
        self.collect()
        with zipfile.ZipFile(self.output) as archive:
            good = {name: archive.read(name) for name in archive.namelist()}
        variants = []
        for mutate in (lambda p: p.update(complete=False), lambda p: p.update(errors=["bad"]),
                       lambda p: p.update(git_head="c"*40), lambda p: p.update(required_analysis_commit="c"*40),
                       lambda p: p.update(validation_profile="other")):
            contents = dict(good)
            manifest = json.loads(contents["manifest.json"])
            mutate(manifest)
            contents["manifest.json"] = json.dumps(manifest).encode()
            variants.append(contents)
        variants += [{**good, "Right_lowe/extra.json": b"{}"}, {**good, "global/" + next(iter(owner.names())): b"changed"}]
        for index, contents in enumerate(variants):
            path = self.transfer / f"bad{index}.zip"
            with zipfile.ZipFile(path, "w") as archive:
                for name, data in contents.items():
                    archive.writestr(name, data)
            with self.subTest(index=index), self.assertRaises(ValueError):
                owner.verify_zip(path, SHA, records)

    @contextmanager
    def flow(self, fail=None):
        """Keep actual status, child mock, artifacts, collector and delivery logic."""
        parent = self.directory / "owned-runtime"
        parent.mkdir(mode=0o700, exist_ok=True)
        worktree = parent / "analysis"
        source = self.directory / "clean-source"
        config = self.directory / "installed.path"
        config.write_bytes(b"installed unchanged")
        installed = {"original_source_sha256": {str(config): owner.sha256(config)}}
        state = {"branch": "test", "head": SHA, "origin_test": SHA, "porcelain": ""}
        baseline = {"LTANAPATH": str(self.repo), "OUTPATH": str(self.outdir),
                    "VOLATILEPATH": str(self.outdir.parents[2]), "ANATYPE": "Kaon"}
        runtime_paths = {**baseline, "LTANAPATH": str(worktree), "OUTPATH": str(worktree / "OUTPUT/Analysis/KaonLT")}
        populate_inputs(self.outdir)
        ordinary_before = {p.name: owner.file_identity(p) for p in self.outdir.iterdir()
                           if p.suffix == ".root" or p.name in owner.names()}
        @contextmanager
        def owned(repo, sha, leaf):
            target = worktree if leaf == "analysis" else source
            target.mkdir()
            try:
                yield target
            finally:
                if leaf == "analysis":
                    shutil.rmtree(parent)
                else:
                    shutil.rmtree(source)
                if fail == "cleanup":
                    raise ValueError("cleanup_failed")
        @contextmanager
        def collection(clean):
            self.assertEqual(clean, source)
            self.assertFalse((clean / "OUTPUT").exists())
            module = Mock()
            def collect(**kw):
                self.assertEqual(kw["repo_root"], source)
                private = worktree.parent / "ab-output/Analysis/KaonLT"
                self.assertEqual(kw["outdir"], private)
                self.assertTrue(private.is_dir())
                if fail == "collector":
                    raise ValueError("collector_failed")
                return self.collect(outdir=private)
            module.collect_validation_bundle.side_effect = collect
            yield module
        def child(tree, log, env):
            self.assertEqual(tree, worktree)
            self.assertEqual(env[owner.ENABLE_ENV], "1")
            log.write_text(MARKERS)
            # Let the real filesystem clock advance before genuinely new writes;
            # never retime artifacts or change the owner's freshness boundary.
            time.sleep(.02)
            private = tree.parent / "ab-output/Analysis/KaonLT"
            artifacts(private)
            (private / "Center_kaon_Analysed_Data_Q4p4W2p74_lowe.root").write_bytes(b"private derived branches")
            return 1 if fail == "child" else 0
        with ExitStack() as stack:
            stack.enter_context(redirect_stderr(io.StringIO()))
            def mock(obj, name, **kwargs):
                return stack.enter_context(patch.object(obj, name, **kwargs))
            mock(owner.isolation.accepted, "REPO", new=self.repo)
            mock(owner.isolation.accepted, "ARTIFACTS", new=self.outdir)
            mock(owner.isolation.accepted, "GLOBUS", new=self.transfer)
            mock(owner, "preflight", return_value=state)
            mock(owner, "resolved_profile", return_value=self.profile)
            source_check = mock(owner, "source_preflight", return_value=([], IDENTITIES))
            mock(owner.isolation, "import_identity", return_value={"package_file": str(config), "pathing_file": str(config)})
            mock(owner.isolation, "owned_worktree", side_effect=owned)
            overlay = mock(owner.isolation, "prepare_runtime_overlay", return_value=({}, {"baseline_paths": baseline}))
            probes = mock(owner.isolation, "probe_paths", return_value=[{"paths": runtime_paths}])
            external = mock(owner, "external_setup", return_value={"background_simc_path": "", "passed": True, "links": []})
            link = mock(owner, "private_output", wraps=owner.private_output)
            mock(owner, "private_callers", return_value=[])
            mock(owner, "skim_inventory", return_value={"skim.root": {"exists": False}})
            copy = mock(owner, "stage_copies", wraps=owner.stage_copies)
            run = mock(owner.debug, "run_analysis", side_effect=child)
            preservation = mock(owner.isolation, "verify_preservation", return_value={"passed": True})
            originals = mock(owner, "verify_originals", wraps=owner.verify_originals)
            accepted_check = mock(owner, "verify_accepted", wraps=owner.verify_accepted)
            mock(owner.isolation, "collection_module", side_effect=collection)
            installed_check = mock(owner.isolation, "verify_ltsep_preservation", wraps=owner.isolation.verify_ltsep_preservation)
            mock(owner, "pdf_pages", return_value=6)
            verify = mock(owner, "verify_zip", wraps=owner.verify_zip)
            if fail == "source": source_check.side_effect = ValueError("source_failed")
            if fail == "overlay": overlay.side_effect = ValueError("copy_import_failed")
            if fail == "paths": probes.side_effect = ValueError("ltsep_path_mismatch")
            if fail == "external": external.side_effect = ValueError("external_link_missing")
            if fail == "output_link": link.side_effect = ValueError("private_output_collision")
            if fail == "copy": copy.side_effect = ValueError("private_copy_identity")
            if fail == "zip": verify.side_effect = ValueError("zip_failed")
            if fail == "preservation": preservation.side_effect = ValueError("ordinary_checkout_changed")
            if fail == "original": originals.side_effect = ValueError("original_input_changed")
            if fail == "accepted": accepted_check.side_effect = ValueError("accepted_artifact_changed")
            if fail == "installed": installed_check.side_effect = ValueError("installed_ltsep_changed")
            try:
                yield run, preservation, verify
            finally:
                self.assertEqual(ordinary_before, {p.name: owner.file_identity(p) for p in self.outdir.iterdir()
                                                  if p.suffix == ".root" or p.name in owner.names()})

    def test_owner_success_once_and_companions(self):
        with self.flow() as (run, preservation, verify):
            self.assertEqual(owner.execute_gate(self.repo, self.outdir, self.output, SHA), self.output)
        run.assert_called_once()
        verify.assert_called_once()
        self.assertEqual(preservation.call_count, 2)
        self.assertFalse((self.directory / "owned-runtime").exists())
        status = json.loads((self.transfer / "comparison-gate-status.json").read_text())
        self.assertEqual(status["status"], "success")
        self.assertEqual(status["diagnostic_outcomes"]["Left"]["status"], "available")
        self.assertTrue((self.transfer / "comparison.log").is_file())
        self.assertTrue((self.transfer / "comparison-run-summary.json").is_file())

    def test_each_failure_stops_delivery_preserves_receipt_and_checkout(self):
        for failure in ("source", "overlay", "paths", "external", "output_link", "copy", "child", "collector", "cleanup", "zip", "preservation", "original", "accepted", "installed"):
            with self.subTest(failure=failure):
                # Each attempt owns separate fresh destinations and target inventory.
                with tempfile.TemporaryDirectory() as fresh:
                    old = self.outdir, self.transfer, self.output
                    self.outdir = Path(fresh) / "volatile/OUTPUT/Analysis/KaonLT"; self.outdir.mkdir(parents=True)
                    self.transfer = Path(fresh) / "transfer"; self.transfer.mkdir()
                    self.output = self.transfer / "attempt.zip"
                    try:
                        with self.flow(failure) as (run, preservation, verify):
                            with self.assertRaises(ValueError):
                                owner.execute_gate(self.repo, self.outdir, self.output, SHA)
                        self.assertGreaterEqual(preservation.call_count, 1)
                        status = json.loads((self.outdir / "attempt-gate-status.json").read_text())
                        self.assertEqual(status["status"], "failed")
                        self.assertFalse((self.transfer / "attempt-gate-status.json").exists())
                        if failure not in ("zip", "preservation", "cleanup"):
                            verify.assert_not_called()
                        if failure in ("child", "cleanup", "collector", "zip", "preservation"):
                            self.assertTrue((self.outdir / "attempt.log").is_file())
                        if failure == "preservation":
                            self.assertIn("original_inputs", status["preservation"])
                            self.assertIn("accepted_artifacts", status["preservation"])
                            self.assertIn("installed_ltsep", status["preservation"])
                    finally:
                        self.outdir, self.transfer, self.output = old

    def test_collisions_and_existing_comparison_are_not_overwritten(self):
        for destination in owner.destinations(self.outdir, self.output).values():
            destination.write_bytes(b"accepted old evidence")
            with self.flow() as (run, _, _), self.assertRaisesRegex(ValueError, "output_collision"):
                owner.execute_gate(self.repo, self.outdir, self.output, SHA)
            run.assert_not_called()
            self.assertEqual(destination.read_bytes(), b"accepted old evidence")
            destination.unlink()
        target = self.outdir / next(iter(owner.names()))
        target.write_bytes(b"accepted old comparison")
        with self.flow() as (run, _, _):
            owner.execute_gate(self.repo, self.outdir, self.output, SHA)
        self.assertEqual(target.read_bytes(), b"accepted old comparison")
        run.assert_called_once()

    def test_wrong_owner_boundary_rejected_before_any_receipt(self):
        with self.assertRaisesRegex(ValueError, "owner_path_boundary_invalid"):
            owner.execute_gate(self.repo, self.outdir, self.output, SHA)
        self.assertFalse(any(self.outdir.iterdir()))


def populate_inputs(root):
    for eps in ("lowe", "highe"):
        for s in ("Center", "Left", "Right"):
            filenames = [f"{s}_kaon_Analysed_{r}_Q4p4W2p74_{eps}.root" for r in ("Data", "Dummy")]
            filenames += [f"Prod_Coin_Q4p4W2p74{s.lower()}_{eps}{suffix}"
                          for suffix in ((".root", ".hist") if eps == "lowe" else (".root",))]
            for name in filenames:
                p = root / name
                if not p.exists():
                    p.write_bytes(("fixture:" + name).encode())


class PrivateTopology(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.work = self.base / "analysis"
        self.work.mkdir()
        self.shared = self.base / "volatile/OUTPUT/Analysis/KaonLT"
        self.shared.mkdir(parents=True)
        self.paths = {"LTANAPATH": str(self.work), "OUTPATH": str(self.work / "OUTPUT/Analysis/KaonLT"),
                      "VOLATILEPATH": str(self.base / "volatile"), "ANATYPE": "Kaon"}

    def private(self):
        return owner.private_output(self.work, self.paths, self.shared)

    def test_private_writes_keep_original_bytes_links_and_independent_inodes(self):
        populate_inputs(self.shared)
        artifacts(self.shared)
        (self.shared / "ordinary.pdf").write_bytes(b"ordinary science")
        (self.shared / "accepted-alias").symlink_to(self.shared / "ordinary.pdf")
        before = {p.name: owner.file_identity(p) for p in self.shared.iterdir()}
        accepted = owner.accepted_identities(self.shared)
        private = self.private()
        plan = owner.copy_plan(self.shared, private)
        self.assertEqual(len(plan), 21)
        self.assertEqual(sum(p["bytes"] for p in plan), 1017)
        copies = {}
        owner.stage_copies(self.shared, private, plan, copies)
        for item in copies.values():
            a, b = item["original"], item["runtime_container_before"]
            self.assertNotEqual((a["device"], a["inode"]), (b["device"], b["inode"]))
        self.assertEqual(owner.effective_csv_path(ROOT, self.paths).resolve(), private)
        (Path(self.paths["OUTPATH"]) / "table.csv").write_text("private CSV")
        (private / "ordinary.pdf").write_bytes(b"new private PDF")
        (private / next(iter(copies))).write_bytes(b"derived ROOT branch stand-in")
        owner.verify_private(self.work, private)
        owner.verify_originals(copies)
        owner.verify_accepted(self.shared, accepted)
        self.assertEqual(before, {p.name: owner.file_identity(p) for p in self.shared.iterdir()})

    def test_optional_right_without_archived_simc_keeps_raw_inputs_and_diamond(self):
        populate_inputs(self.shared)
        (self.shared / "Prod_Coin_Q4p4W2p74right_lowe.root").unlink()
        # Right is unavailable downstream, but its raw files remain Diamond inputs.
        extra = self.shared / "Right_kaon_Q4p4W2p74_mide.root"
        extra.write_bytes(b"raw Diamond control")
        artifacts(self.shared)
        before = {p.name: owner.file_identity(p) for p in self.shared.iterdir()}
        private = self.private()
        plan = owner.copy_plan(self.shared, private)
        planned = {Path(p["source"]).name for p in plan}
        self.assertIn(extra.name, planned)
        self.assertNotIn("Prod_Coin_Q4p4W2p74right_lowe.hist", planned)
        for eps in ("lowe", "highe"):
            for s in ("Center", "Left", "Right"):
                for role in ("Data", "Dummy"):
                    self.assertIn(f"{s}_kaon_Analysed_{role}_Q4p4W2p74_{eps}.root", planned)
        for s in ("center", "left"):
            for suffix in (".root", ".hist"):
                self.assertIn(f"Prod_Coin_Q4p4W2p74{s}_lowe{suffix}", planned)
        copies = {}
        owner.stage_copies(self.shared, private, plan, copies)
        self.assertEqual(owner.diamond_choices(self.shared), owner.diamond_choices(private))
        # Execute only frozen shift_prep's early availability checks, then the
        # frozen main.py availability filter; no ROOT or production import.
        tree = ast.parse((ROOT / "src/cuts/shift_prep.py").read_text())
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "shift_prep")
        guard = next(i for i, n in enumerate(function.body) if isinstance(n, ast.If) and
                     ast.unparse(n.test) == "not os.path.exists(simc_filename)")
        check = ast.FunctionDef(name="available", args=deepcopy(function.args),
            body=deepcopy(function.body[:guard+1]) + [ast.Return(ast.Constant(True))], decorator_list=[])
        compiled = compile(ast.fix_missing_locations(ast.Module(body=[check], type_ignores=[])), "<frozen-availability>", "exec")
        main = ast.parse((ROOT / "src/main.py").read_text())
        selector = next(n for n in main.body if isinstance(n, ast.Assign) and
                        any(isinstance(t, ast.Name) and t.id == "available_phisetlist" for t in n.targets))
        for root in (self.shared, private):
            ns = {"OUTPATH": str(root), "os": os, "_to_float": lambda v: float(v.replace("p", ".")), "print": lambda *a: None}
            exec(compiled, ns)
            for eps in ("low", "high"):
                config = {"ParticleType": "kaon", "Q2": "4p4", "W": "2p74", "EPSSET": eps,
                          "mm_min": .7, "mm_max": 1.3, "tmin": 0., "tmax": 1.,
                          "InDATAFilename": f"Analysed_Data_Q4p4W2p74_{eps}e",
                          "InDUMMYFilename": f"Analysed_Dummy_Q4p4W2p74_{eps}e"}
                available = {s: True for s in ("Center", "Left", "Right") if ns["available"](s, config) is True}
                state = {"phisetlist": ["Center", "Left", "Right"], "mm_shift_summary": available,
                         "kaon_data_coordinate_summary": available}
                exec(compile(ast.Module(body=[selector], type_ignores=[]), "<frozen-main-filter>", "exec"), state)
                self.assertEqual(state["available_phisetlist"], ["Center", "Left"] if eps == "low" else ["Center", "Left", "Right"])
        self.assertEqual(before, {p.name: owner.file_identity(p) for p in self.shared.iterdir()})

    def test_optional_high_simc_absence_and_missing_optional_dummy_preserve_raw_scope(self):
        populate_inputs(self.shared)
        private = self.private()
        for setting in ("Center", "Left", "Right"):
            p = self.shared / f"Prod_Coin_Q4p4W2p74{setting.lower()}_highe.root"
            raw = p.read_bytes(); p.unlink()
            # Without SIMC shift_prep returns before consuming the dummy file.
            dummy = self.shared / "Right_kaon_Analysed_Dummy_Q4p4W2p74_highe.root"
            if setting == "Right":
                dummy.unlink()
            with self.subTest(setting=setting):
                plan = owner.copy_plan(self.shared, private)
                planned = {Path(p["source"]).name for p in plan}
                expected = {p.name for p in self.shared.glob("*kaon*Q4p4W2p74*.root") if "simc" not in p.name.lower()}
                self.assertTrue(expected.issubset(planned))
                for available in {"Center", "Left", "Right"} - {setting}:
                    self.assertIn(f"Prod_Coin_Q4p4W2p74{available.lower()}_highe.root", planned)
                for s in ("center", "left"):
                    for suffix in (".root", ".hist"):
                        self.assertIn(f"Prod_Coin_Q4p4W2p74{s}_lowe{suffix}", planned)
                if setting == "Right":
                    owner.stage_copies(self.shared, private, plan, {})
                    self.assertEqual(owner.diamond_choices(self.shared), owner.diamond_choices(private))
            p.write_bytes(raw)

    def test_mandatory_low_simc_and_available_optional_hist_remain_strict(self):
        populate_inputs(self.shared)
        private = self.private()
        missing = [f"Prod_Coin_Q4p4W2p74{s}_lowe{suffix}" for s in ("center", "left")
                   for suffix in (".root", ".hist")]
        missing.append("Prod_Coin_Q4p4W2p74right_lowe.hist")
        for name in missing:
            p = self.shared / name
            raw = p.read_bytes(); p.unlink()
            with self.subTest(missing=name), self.assertRaisesRegex(ValueError, "required_simc_input_missing:" + name):
                owner.copy_plan(self.shared, private)
            p.write_bytes(raw)
        # Existing invalid optional SIMC is not treated as an absent setting.
        p = self.shared / "Prod_Coin_Q4p4W2p74right_lowe.root"
        p.unlink(); p.mkdir()
        with self.assertRaisesRegex(ValueError, "required_simc_input_missing"):
            owner.copy_plan(self.shared, private)

    def test_existing_wrong_and_dangling_links_never_replaced(self):
        for leaf in (self.work / "OUTPUT", self.base / "ab-output"):
            for kind in ("file", "directory", "correct", "wrong", "dangling"):
                with self.subTest(leaf=leaf.name, kind=kind):
                    if kind == "file": leaf.write_text("preserve")
                    elif kind == "directory": leaf.mkdir()
                    else: leaf.symlink_to(self.shared if kind != "dangling" else self.base / "missing", target_is_directory=True)
                    with self.assertRaisesRegex(ValueError, "collision"):
                        self.private()
                    self.assertTrue(os.path.lexists(leaf))
                    if leaf.is_dir() and not leaf.is_symlink(): leaf.rmdir()
                    else: leaf.unlink()
        bad = {**self.paths, "OUTPATH": str(self.shared)}
        with self.assertRaisesRegex(ValueError, "authority"):
            owner.private_output(self.work, bad, self.shared)
        alias = self.base / "alias"; alias.symlink_to(self.work)
        with self.assertRaisesRegex(ValueError, "worktree_invalid"):
            owner.private_output(alias, self.paths, self.shared)

    def test_changed_private_link_and_unsafe_parent_stop(self):
        self.base.chmod(0o755)
        with self.assertRaisesRegex(ValueError, "permissions"):
            self.private()
        self.base.chmod(0o700)
        private = self.private()
        (self.work / "OUTPUT").unlink()
        (self.work / "OUTPUT").symlink_to(self.shared)
        with self.assertRaisesRegex(ValueError, "private_link_changed"):
            owner.verify_private(self.work, private)

    def test_external_role_alias_preserves_canonical_authority(self):
        alias = self.base / "volatile-alias"
        alias.symlink_to(self.base / "volatile")
        paths = {**self.paths, "VOLATILEPATH": str(alias)}
        private = owner.private_output(self.work, paths, alias / "OUTPUT/Analysis/KaonLT")
        owner.verify_private(self.work, private)
        self.assertNotEqual(private, self.shared.resolve())

    def test_missing_input_collision_path_escape_and_changed_plan_fail_before_copy(self):
        populate_inputs(self.shared)
        private = self.private()
        missing = self.shared / "Center_kaon_Analysed_Dummy_Q4p4W2p74_highe.root"
        raw = missing.read_bytes(); missing.unlink()
        with self.assertRaisesRegex(ValueError, "incomplete_data_dummy"):
            owner.copy_plan(self.shared, private)
        missing.write_bytes(raw)
        p = private / missing.name; p.symlink_to(self.base / "absent")
        with self.assertRaisesRegex(ValueError, "collision"):
            owner.copy_plan(self.shared, private)
        p.unlink()
        plan = owner.copy_plan(self.shared, private)
        plan[0]["target"] = str(self.base / "escape")
        with self.assertRaisesRegex(ValueError, "plan_changed"):
            owner.stage_copies(self.shared, private, plan, {})
        self.assertEqual(list(private.iterdir()), [])

    def test_diamond_choices_match_frozen_pure_selectors_and_ties(self):
        import glob
        populate_inputs(self.shared)
        for name in ("Center_kaon_Q4p4W2p74_lowe_B.root", "Center_kaon_Q4p4W2p74_lowe_A.root",
                     "center_kaon_Q4p4W2p74_mide.root", "Center_kaon_ShiftPrep_Q4p4W2p74_lowe.root"):
            (self.shared / name).write_bytes(b"diamond fixture")
        private = self.private()
        owner.stage_copies(self.shared, private, owner.copy_plan(self.shared, private), {})
        tree = ast.parse((ROOT / "src/cuts/diamond.py").read_text())
        keep = {"_is_shiftprep_artifact", "_find_shortest_root_any", "_find_center_fallback"}
        functions = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name in keep}
        for root in (self.shared, private):
            ns = {"glob": glob, "os": os, "OUTPATH": str(root), "ParticleType": "kaon", "FilenameOverride": "Q4p4W2p74"}
            for n in functions.values():
                exec(compile(ast.Module(body=[n], type_ignores=[]), "<pure-diamond-selector>", "exec"), ns)
            for eps in ("low", "mid", "high"):
                for setting in ("Center", "Left", "Right"):
                    chosen = ns["_find_shortest_root_any"]([setting, setting.lower()], eps)
                    self.assertEqual(owner.diamond_choices(root)[f"overlay:{setting}:{eps}"], None if chosen is None else Path(chosen).name)
                chosen = ns["_find_center_fallback"](eps)
                self.assertEqual(owner.diamond_choices(root)[f"fallback:Center:{eps}"], None if chosen is None else Path(chosen).name)

    def test_skim_predicate_is_existence_only_and_preserves_absence(self):
        skim = self.base / "skims"; skim.mkdir()
        util = self.base / "util"; lists = util / "run_list/KaonLT"; lists.mkdir(parents=True)
        (lists / "Q4p4W2p74center_lowe").write_text("1\n2\n")
        (skim / "kaon_1_-1_Raw_Data.root").mkdir()  # exists() is true even for a directory
        with patch.object(owner, "skim_directories", return_value=(str(skim), str(skim))), patch.object(
                owner, "sha256", side_effect=AssertionError("skim byte hashing forbidden")):
            result = owner.skim_inventory(ROOT, {"UTILPATH": str(util)})
        self.assertTrue(result[str(skim / "kaon_1_-1_Raw_Data.root")]["exists"])
        self.assertFalse(result[str(skim / "kaon_2_-1_Raw_Data.root")]["exists"])

    def test_external_writer_guard_and_independent_cache_locations(self):
        with patch.object(owner.isolation, "external_symlink_preflight", return_value={"background_simc_path": "", "links": []}):
            with self.assertRaisesRegex(ValueError, "background_mkdir"):
                owner.external_setup(ROOT, self.paths, self.shared)
        env = owner.private_environment({"VOLATILEPATH": "stable"}, self.work)
        self.assertEqual(env["VOLATILEPATH"], "stable")
        for key in ("TMPDIR", "XDG_CACHE_HOME", "MPLCONFIGDIR", "PYTHONPYCACHEPREFIX"):
            self.assertTrue(Path(env[key]).is_relative_to(self.base))
            self.assertTrue(Path(env[key]).is_dir())

    def test_default_plot_callers_aliases_and_wrong_writer_fail_closed(self):
        private = self.private()
        ordinary = self.base / "ordinary-checkout"
        (ordinary / "src/setup").mkdir(parents=True)
        (ordinary / "OUTPUT").symlink_to(self.base / "volatile/OUTPUT")
        shutil.copyfile(ROOT / "src/setup/findEffectiveCharge.py", ordinary / "src/setup/findEffectiveCharge.py")
        (self.work / "src/setup").mkdir(parents=True)
        shutil.copyfile(ROOT / "src/setup/findEffectiveCharge.py", self.work / "src/setup/findEffectiveCharge.py")
        shutil.copyfile(ROOT / "src/setup/getDataTable.py", self.work / "src/setup/getDataTable.py")
        baseline = {k: "stable-"+k for k in owner.isolation.PATH_FIELDS}
        baseline.update(LTANAPATH=str(ordinary), OUTPATH=str(self.shared), ANATYPE="Kaon")
        runtime = {**baseline, "LTANAPATH": str(self.work), "OUTPATH": str(self.work / "OUTPUT/Analysis/KaonLT")}
        installed = self.base / "installed/ltsep"
        copied = self.base / "runtime-python/ltsep"
        copied.mkdir(parents=True)
        alias = self.base / "package-alias"; alias.symlink_to(copied)
        calls = [("src/main.py", ("Plot_LTSep",), "src"),
                 ("src/setup/findEffectiveCharge.py", (), "src/setup"),
                 ("src/setup/getDataTable.py", (), "src/setup")]
        def probe(cwd, code, arguments, env):
            is_runtime = Path(arguments[0]).is_relative_to(self.work)
            package = alias if is_runtime else installed
            fields = dict(runtime if is_runtime else baseline)
            if not json.loads(arguments[1]):
                fields["OUTPATH"] = str(self.shared)  # actual no-mode field may remain shared
            return {"paths": fields, "package": str(package / "__init__.py"), "pathing": str(package / "pathing.py")}
        with patch.object(owner, "reached_callers", return_value=calls), patch.object(owner.isolation, "python_json", side_effect=probe):
            records = owner.private_callers(ordinary, self.work, {}, baseline, private, installed)
        self.assertEqual(len(records), 3)
        runtime["OUTPATH"] = str(self.shared)
        with patch.object(owner, "reached_callers", return_value=calls), patch.object(owner.isolation, "python_json", side_effect=probe):
            with self.assertRaisesRegex(ValueError, "private_outpath"):
                owner.private_callers(ordinary, self.work, {}, baseline, private, installed)
        self.assertEqual(len(owner.reached_callers(ROOT)), 24)

    def test_copy_failure_keeps_original_snapshot_and_source_symlink_identity(self):
        populate_inputs(self.shared)
        private = self.private()
        plan = owner.copy_plan(self.shared, private)
        copies = {}
        with patch.object(owner.shutil, "copyfileobj", side_effect=OSError("copy failure")):
            with self.assertRaisesRegex(OSError, "copy failure"):
                owner.stage_copies(self.shared, private, plan, copies)
        self.assertEqual(len(copies), 1)
        owner.verify_originals(copies)
        target = self.shared / "identity-alias"
        target.symlink_to(Path(plan[0]["source"]))
        before = owner.file_identity(target)
        target.unlink(); target.symlink_to(Path(plan[0]["source"]).relative_to(self.shared))
        self.assertNotEqual(before, owner.file_identity(target))


if __name__ == "__main__":
    unittest.main()
