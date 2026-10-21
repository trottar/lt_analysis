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

from testing import run_pion_local_t_event_validation_gate as owner

ROOT = Path(__file__).resolve().parents[1]
SHA = "a" * 40
FP = "b" * 64
MARKERS = ("All available low-epsilon settings completed in debug mode.\n"
           "Full high-epsilon processing is intentionally skipped in -d debug mode.\n")
IDENTITIES = {name: {"sha256": FP, "bytes": 12} for name in owner.SOURCE_FILES}


from testing.test_pion_local_t_event_validation import fixture
from testing import pion_local_t_event_validation_capture as capture


def artifacts(outdir):
    data, dummy = outdir / "source-data.bin", outdir / "source-dummy.bin"
    data.write_bytes(b"synthetic data source"); dummy.write_bytes(b"synthetic dummy source")
    identities = {"data": capture.file_identity(data), "dummy": capture.file_identity(dummy)}
    for setting in ("Center", "Left"):
        packet, _, _ = fixture(setting)
        packet["input_identity"] = identities
        packet["source_input_paths"] = {"data": str(data), "dummy": str(dummy)}
        packet["source_files"] = deepcopy(IDENTITIES)
        packet["source_HEAD"] = SHA
        for population, rows in packet["events"].items():
            for row in rows:
                role = "dummy" if row["source_label"].startswith("dummy") else "data"
                row["source_sha256"] = identities[role]["sha256"]
                row["fold"] = capture.fold(row)
            packet["event_hashes"][population] = capture.digest(rows)
        (outdir / (setting+"_kaon_Q4p4W2p74_lowe_pion-local-t-events.json")).write_bytes(capture.encoded(packet))


class GateTests(unittest.TestCase):
    def test_different_but_self_consistent_source_identity_rejected_before_collection(self):
        artifacts(self.outdir)
        path=self.outdir / "Center_kaon_Q4p4W2p74_lowe_pion-local-t-events.json"
        packet=json.loads(path.read_bytes())
        packet["input_identity"]["data"]={"sha256":"e"*64,"bytes":100}
        for population,rows in packet["events"].items():
            for row in rows:
                if row["source_label"] in ("prompt","rand"):
                    row["source_sha256"]="e"*64
                    row["fold"]=capture.fold(row)
            packet["event_hashes"][population]=capture.digest(rows)
        path.write_bytes(capture.encoded(packet))
        capture.validate_packet(packet)
        with patch.object(capture,"file_identity",side_effect=AssertionError("late source-file check")):
            with self.assertRaisesRegex(ValueError,"accepted_AB_input_identity_mismatch:Center:data"):
                self.verify()

    def test_environment_only_one_opt_in_and_completion_gate(self):
        clean = owner.child_environment({"LT_ANALYSIS_DEBUG_LEFT_LOW":"1", "KEEP":"value"})
        self.assertEqual(clean, {owner.ENABLE_ENV:"1", "KEEP":"value"})
        for env in ({owner.ENABLE_ENV:"0"}, {"LT_ANALYSIS_PION_T_TEMPLATE_COMPARISON":"1"}, {"LT_ANALYSIS_UNREVIEWED":"1"}):
            with self.assertRaises(ValueError):owner.child_environment(env)
        log = self.directory / "completion.log"
        log.write_text(MARKERS)
        owner.verify_completion(log,0)
        with self.assertRaisesRegex(ValueError,"analysis_failed"):owner.verify_completion(log,1)
        for value in ("",MARKERS+"High Epsilon Completed!"):
            log.write_text(value)
            with self.assertRaises(ValueError):owner.verify_completion(log,0)

    def test_packet_provenance_completeness_and_source_data_checks(self):
        artifacts(self.outdir)
        path = self.outdir / "Center_kaon_Q4p4W2p74_lowe_pion-local-t-events.json"
        good = json.loads(path.read_bytes())
        cases = [(lambda p:p.update(source_HEAD="c"*40),"source_or_setting"),
                 (lambda p:p["source_files"][owner.SOURCE_FILES[0]].update(sha256="c"*64),"committed_source_hash"),
                 (lambda p:p["parents"].pop(),"parent_inventory"),
                 (lambda p:p["source_counts"]["pion"].update(prompt=0),"missing_source_rows"),
                 (lambda p:p["parents"][0]["pion_input"]["values"].__setitem__(1,100.),"histogram_fingerprint")]
        for mutate,reason in cases:
            packet=deepcopy(good);mutate(packet);path.write_bytes(capture.encoded(packet))
            with self.subTest(reason=reason),self.assertRaisesRegex(ValueError,reason):self.verify()
        path.write_bytes(capture.encoded(good))
        (self.outdir / "source-data.bin").write_bytes(b"changed data")
        with self.assertRaisesRegex(ValueError,"source_data_changed"):self.verify()

    def test_packet_freshness_and_missing_setting(self):
        artifacts(self.outdir)
        with self.assertRaisesRegex(ValueError,"packet_stale"):
            owner.verify_artifacts(self.outdir,SHA,IDENTITIES,time.time_ns()+10**9,self.before)
        next(self.outdir.glob("Left*.json")).unlink()
        with self.assertRaisesRegex(ValueError,"packet_missing_or_empty"):self.verify()

    def setUp(self):
        # Synthetic source bytes stand in for the independently pinned farm inputs.
        synthetic_inputs = {role: {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
                            for role,raw in (("data",b"synthetic data source"),("dummy",b"synthetic dummy source"))}
        pins = patch.object(capture, "ACCEPTED_INPUT_IDENTITIES",
                            {s:deepcopy(synthetic_inputs) for s in ("Center","Left")})
        pins.start()
        self.addCleanup(pins.stop)
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
        return owner.verify_artifacts(self.outdir, SHA, IDENTITIES, 0, self.before)

    def test_exact_profile_and_unchanged_collector_acceptance(self):
        tracked = owner.collector.load_validation_profile(ROOT / owner.PROFILE)
        self.assertEqual(tracked["artifacts"], owner.ARTIFACTS)
        self.assertEqual(tracked["settings"], owner.SETTINGS)
        self.assertEqual(tracked["artifacts"]["settings"], [])
        self.assertEqual(len(owner.names()), 2)
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
            self.assertEqual(len(archive.namelist()), 11)
            self.assertEqual(sum(not name.endswith("/") and name.startswith("global/") for name in archive.namelist()), 2)
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
        self.assertEqual(status["diagnostic_outcomes"]["Left"]["status"], "captured_pending_independent_validation")
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
