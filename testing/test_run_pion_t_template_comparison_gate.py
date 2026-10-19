"""Synthetic standard-library fixtures; no production launch, ROOT or network."""
from contextlib import contextmanager, ExitStack
from copy import deepcopy
import hashlib
import io
import json
import os
from pathlib import Path
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
        self.outdir = self.directory / "artifacts"
        self.outdir.mkdir()
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
        self.assertEqual(tracked["source_identity"]["required_analysis_commit"], owner.BASE_HEAD)
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
        with patch.object(owner.isolation, "snapshot", return_value=state), patch.object(owner, "git", side_effect=["", "src/main.py"]):
            with self.assertRaisesRegex(ValueError, "unexpected_committed_source"):
                owner.preflight(self.repo, SHA)
        for sha in ("short", "A"*40, owner.BASE_HEAD):
            with patch.object(owner.isolation, "snapshot", return_value={**state, "head": sha, "origin_test": sha}):
                with self.assertRaises(ValueError):
                    owner.preflight(self.repo, sha)

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

    def collect(self, output=None, profile=None):
        profile_path = self.directory / "effective.json"
        profile_path.write_text(json.dumps(profile or self.profile))
        module = owner.collector
        with patch.object(module, "collect_source_state", return_value=("source fixture", SHA)), patch.object(module, "collect_source_checks", return_value=("checks fixture", [])), patch.object(module, "_committed_identity", return_value=(True, [], [])):
            return module.collect_validation_bundle(outdir=self.outdir, kinematic="Q4p4W2p74",
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
        worktree = self.directory / "owned-analysis"
        worktree.mkdir(exist_ok=True)
        config = self.directory / "installed.path"
        config.write_bytes(b"installed unchanged")
        installed = {"original_source_sha256": {str(config): owner.sha256(config)}}
        state = {"branch": "test", "head": SHA, "origin_test": SHA, "porcelain": ""}
        @contextmanager
        def owned(repo, sha, leaf):
            yield worktree
            if fail == "cleanup":
                raise ValueError("cleanup_failed")
        @contextmanager
        def collection(source):
            module = Mock()
            module.collect_validation_bundle.side_effect = lambda **kw: self.collect()
            yield module
        def child(tree, log, env):
            self.assertEqual(tree, worktree)
            self.assertEqual(env[owner.ENABLE_ENV], "1")
            log.write_text(MARKERS)
            artifacts(self.outdir)
            # Explicit fresh fixture timestamps avoid filesystem clock rounding.
            fixture_ns = time.time_ns() + 10**9
            for name in owner.names():
                os.utime(self.outdir / name, ns=(fixture_ns, fixture_ns))
            return 1 if fail == "child" else 0
        with ExitStack() as stack:
            def mock(obj, name, **kwargs):
                return stack.enter_context(patch.object(obj, name, **kwargs))
            mock(owner.isolation.accepted, "REPO", new=self.repo)
            mock(owner.isolation.accepted, "ARTIFACTS", new=self.outdir)
            mock(owner.isolation.accepted, "GLOBUS", new=self.transfer)
            mock(owner, "preflight", return_value=state)
            mock(owner, "resolved_profile", return_value=self.profile)
            source = mock(owner, "source_preflight", return_value=([], IDENTITIES))
            mock(owner.isolation, "import_identity", return_value={"package_file": str(config), "pathing_file": str(config)})
            mock(owner.isolation, "owned_worktree", side_effect=owned)
            overlay = mock(owner.isolation, "prepare_runtime_overlay", return_value=({}, {"baseline_paths": {}}))
            probes = mock(owner.isolation, "probe_paths", return_value=[{"paths": {}}])
            external = mock(owner.isolation, "external_symlink_preflight", return_value={"background_simc_path": "", "passed": True})
            link = mock(owner.debug, "prepare_debug_output_link", return_value={"passed": True})
            run = mock(owner.debug, "run_analysis", side_effect=child)
            preservation = mock(owner.isolation, "verify_preservation", return_value={"passed": True})
            mock(owner.isolation, "collection_module", side_effect=collection)
            mock(owner, "pdf_pages", return_value=6)
            verify = mock(owner, "verify_zip", wraps=owner.verify_zip)
            if fail == "source": source.side_effect = ValueError("source_failed")
            if fail == "overlay": overlay.side_effect = ValueError("copy_import_failed")
            if fail == "paths": probes.side_effect = ValueError("ltsep_path_mismatch")
            if fail == "external": external.side_effect = ValueError("external_link_missing")
            if fail == "output_link": link.side_effect = ValueError("debug_output_path_already_exists")
            if fail == "zip": verify.side_effect = ValueError("zip_failed")
            if fail == "preservation": preservation.side_effect = ValueError("ordinary_checkout_changed")
            yield run, preservation, verify

    def test_owner_success_once_and_companions(self):
        with self.flow() as (run, preservation, verify):
            self.assertEqual(owner.execute_gate(self.repo, self.outdir, self.output, SHA), self.output)
        run.assert_called_once()
        verify.assert_called_once()
        self.assertEqual(preservation.call_count, 2)
        status = json.loads((self.transfer / "comparison-gate-status.json").read_text())
        self.assertEqual(status["status"], "success")
        self.assertEqual(status["diagnostic_outcomes"]["Left"]["status"], "available")
        self.assertTrue((self.transfer / "comparison.log").is_file())
        self.assertTrue((self.transfer / "comparison-run-summary.json").is_file())

    def test_each_failure_stops_delivery_preserves_receipt_and_checkout(self):
        for failure in ("source", "overlay", "paths", "external", "output_link", "child", "cleanup", "zip", "preservation"):
            with self.subTest(failure=failure):
                # Each attempt owns separate fresh destinations and target inventory.
                with tempfile.TemporaryDirectory() as fresh:
                    old = self.outdir, self.transfer, self.output
                    self.outdir = Path(fresh) / "artifacts"; self.outdir.mkdir()
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
                        if failure not in ("zip", "preservation"):
                            verify.assert_not_called()
                        if failure in ("child", "cleanup", "zip", "preservation"):
                            self.assertTrue((self.outdir / "attempt.log").is_file())
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
        with self.flow() as (run, _, _), self.assertRaisesRegex(ValueError, "comparison_artifact_preexists"):
            owner.execute_gate(self.repo, self.outdir, self.output, SHA)
        self.assertEqual(target.read_bytes(), b"accepted old comparison")
        run.assert_not_called()

    def test_wrong_owner_boundary_rejected_before_any_receipt(self):
        with self.assertRaisesRegex(ValueError, "owner_path_boundary_invalid"):
            owner.execute_gate(self.repo, self.outdir, self.output, SHA)
        self.assertFalse(any(self.outdir.iterdir()))


if __name__ == "__main__":
    unittest.main()
