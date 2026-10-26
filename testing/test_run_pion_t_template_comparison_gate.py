"""Synthetic standard-library fixtures; no production launch, ROOT or network."""
from contextlib import contextmanager, ExitStack, redirect_stderr
import ast
from copy import deepcopy
import hashlib
import io
import json
import os
import stat
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


@contextmanager
def infrastructure_aliases(base, infra):
    """Real logical /volatile -> /lustre/enp/volatile -> expphy chain."""
    physical = base.parents[2]
    infra.mkdir()
    (infra / "volatile").symlink_to("lustre/enp/volatile")
    (infra / "lustre").symlink_to("lustre24")
    (infra / "lustre24").symlink_to(physical)
    (physical / "enp").symlink_to("expphy")
    modes = {p: stat.S_IMODE(p.stat().st_mode) for p in (infra, physical)}
    for p in modes:
        p.chmod(0o555)
    try:
        yield infra / "volatile", infra / "volatile" / base.name
    finally:
        for p, mode in modes.items():
            p.chmod(mode)
        (physical / "enp").unlink()
        shutil.rmtree(infra)


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
        self.volatile_temp = tempfile.TemporaryDirectory(dir="/dev/shm")
        self.addCleanup(self.volatile_temp.cleanup)
        self.repo = self.directory / "group/repo"
        self.repo.mkdir(parents=True)
        self.outdir = Path(self.volatile_temp.name) / "expphy/volatile/user/OUTPUT/Analysis/KaonLT"
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
        with patch.object(owner.isolation, "snapshot", return_value=state), patch.object(owner, "git", side_effect=lambda repo, *args: owner.CLEANUP_REPAIR_BLOBS[args[1].split(":", 1)[1]] if args[0] == "rev-parse" else "") as run:
            self.assertEqual(owner.preflight(ROOT, SHA), state)
            self.assertIn((ROOT, "merge-base", "--is-ancestor", owner.BASE_HEAD, owner.REVIEWED_HEAD),
                          [call.args for call in run.call_args_list])
            self.assertIn((ROOT, "diff", "--name-only", owner.REVIEWED_HEAD, SHA), [c.args for c in run.call_args_list])
        for path in extra + ["unrelated.py"]:
            with self.subTest(path=path), patch.object(owner.isolation, "snapshot", return_value=state), patch.object(
                    owner, "git", side_effect=["", "", "", path]):
                with self.assertRaisesRegex(ValueError, "unexpected_committed_source"):
                    owner.preflight(ROOT, SHA)

    def test_cleanup_successor_source_guard_is_exact_and_pinned(self):
        state = {"branch": "test", "head": SHA, "origin_test": SHA, "porcelain": ""}
        for defect in (None, "helper", "test", "unrelated", "ancestor", "dirty"):
            def git(repo, *args):
                if args[0] == "diff":
                    return "unrelated.py" if defect == "unrelated" else "\n".join(owner.CLEANUP_REPAIR_BLOBS)
                if args[0] == "rev-parse":
                    path = args[1].split(":", 1)[1]
                    if defect == ("test" if "/test_" in path else "helper"): return "0" * 40
                    return owner.CLEANUP_REPAIR_BLOBS[path]
                if args[:3] == ("merge-base", "--is-ancestor", owner.CLEANUP_REPAIR_BASE) and defect == "ancestor":
                    raise ValueError("not_cleanup_successor")
                return ""
            observed = {**state, "porcelain": " M testing/run_e8_4_fix5_canonical_five_plot_gate.py"} if defect == "dirty" else state
            with self.subTest(defect=defect), patch.object(owner.isolation, "snapshot", return_value=observed), patch.object(owner, "git", side_effect=git):
                if defect:
                    with self.assertRaisesRegex(ValueError, "blob_changed|unexpected_committed_source|not_cleanup_successor|dirty_gate_source"):
                        owner.preflight(self.repo, SHA)
                else:
                    self.assertEqual(owner.preflight(self.repo, SHA), state)
        for path, blob in owner.CLEANUP_REPAIR_BLOBS.items():
            raw = (ROOT / path).read_bytes().replace(b"\r\n", b"\n")
            self.assertEqual(hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest(), blob)

    def test_split_entry_with_real_group_selected_package_and_private_data(self):
        with VolatileRuntime().fixture(aliases=True) as (repo, base, shared, package, config, calls), ExitStack() as stack:
            old = self.repo, self.outdir, self.transfer, self.output
            self.repo, self.outdir = repo, shared
            self.transfer = base / "transfer"; self.transfer.mkdir()
            self.output = self.transfer / "split.zip"
            populate_inputs(shared); artifacts(shared)
            previous = [base / ("kaonlt-ab-analysis-" + x) for x in ("c16r9oxs", "0o8v_67p")]
            for path in previous: path.mkdir(); (path / "KEEP").write_bytes(b"old attempt")
            caller_paths = [("src/main.py", ("Plot_LTSep",), "src"),
                            ("src/setup/findEffectiveCharge.py", (), "src/setup"),
                            ("src/setup/getDataTable.py", (), "src/setup")]
            def copy_callers(tree):
                for caller, _, _ in caller_paths:
                    target = tree / caller; target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes((ROOT / caller).read_bytes())
            copy_callers(repo)
            real_git = owner.isolation.git
            def git(tree, *args):
                result = real_git(tree, *args)
                if args[:2] == ("worktree", "add"): copy_callers(Path(args[3]))
                return result
            def mock(obj, name, **kwargs): return stack.enter_context(patch.object(obj, name, **kwargs))
            mock(owner.isolation, "git", side_effect=git)
            mock(owner.isolation.accepted, "REPO", new=repo)
            mock(owner.isolation.accepted, "ARTIFACTS", new=shared)
            mock(owner.isolation.accepted, "GLOBUS", new=self.transfer)
            source = mock(owner, "source_preflight", return_value=([], IDENTITIES))
            mock(owner, "preflight", return_value={"branch": "test", "head": SHA, "origin_test": SHA, "porcelain": ""})
            mock(owner, "resolved_profile", return_value=self.profile)
            mock(owner, "reached_callers", return_value=caller_paths)
            mock(owner, "external_setup", return_value={"links": [], "directories": []})
            mock(owner, "skim_inventory", return_value={"fixture": {"exists": False}})
            mock(owner.isolation, "verify_preservation", return_value={"passed": True})
            mock(owner, "pdf_pages", return_value=6)
            def child(tree, log, env):
                self.assertTrue(source.called)
                private = (tree / "OUTPUT/Analysis/KaonLT").resolve()
                self.assertTrue(tree.is_relative_to(repo.parent))
                self.assertTrue(private.is_relative_to(base))
                self.assertNotEqual(tree.stat().st_dev, private.stat().st_dev)
                for key in ("TMPDIR", "XDG_CACHE_HOME", "MPLCONFIGDIR", "PYTHONPYCACHEPREFIX"):
                    self.assertTrue(Path(env[key]).is_relative_to(base))
                self.assertFalse(list(tree.parent.rglob("*.root")))
                self.assertFalse(list(tree.parent.rglob("*.hist")))
                self.assertFalse((tree.parent / "ab-output").exists())
                for setting in ("Center", "Left", "Right"):
                    self.assertEqual(len(open_frozen_pion_inputs(private, setting)), 2)
                log.write_text(MARKERS); time.sleep(.02); artifacts(private)
                (private / "Center_kaon_Analysed_Data_Q4p4W2p74_lowe.root").write_bytes(b"derived private branches")
                return 0
            run = mock(owner.debug, "run_analysis", side_effect=child)
            @contextmanager
            def collection(clean):
                module = Mock()
                def collect(**kw):
                    self.assertEqual(kw["repo_root"], clean)
                    self.assertFalse((clean / "OUTPUT").exists())
                    self.assertTrue(Path(kw["outdir"]).is_relative_to(base))
                    self.assertTrue(Path(kw["outdir"]).is_dir())
                    return self.collect(outdir=kw["outdir"])
                module.collect_validation_bundle.side_effect = collect
                yield module
            mock(owner.isolation, "collection_module", side_effect=collection)
            stack.enter_context(redirect_stderr(io.StringIO()))
            try:
                owner.execute_gate(repo, shared, self.output, SHA)
                run.assert_called_once()
                summary = json.loads((self.transfer / "split-run-summary.json").read_text())
                placement = summary["runtime_placement"]
                self.assertFalse(Path(placement["parent"]["path"]).exists())
                self.assertFalse(Path(placement["code_parent"]["path"]).exists())
                self.assertEqual(len(summary["copy_plan"]), 27)
                for path in previous: self.assertEqual((path / "KEEP").read_bytes(), b"old attempt")
            finally:
                self.repo, self.outdir, self.transfer, self.output = old

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
    def flow(self, fail=None, aliases=False):
        """Keep actual status, child mock, artifacts, collector and delivery logic."""
        worktree = None
        volatile = self.outdir.parents[2]
        (self.repo / "OUTPUT").symlink_to(volatile / "OUTPUT") if not (self.repo / "OUTPUT").is_symlink() else None
        source = self.directory / "clean-source"
        config = self.directory / "installed.path"
        config.write_bytes(b"installed unchanged")
        installed = {"original_source_sha256": {str(config): owner.sha256(config)}}
        state = {"branch": "test", "head": SHA, "origin_test": SHA, "porcelain": ""}
        baseline = {"LTANAPATH": str(self.repo), "OUTPATH": str(self.outdir),
                    "VOLATILEPATH": str(self.outdir.parents[2]), "ANATYPE": "Kaon"}
        populate_inputs(self.outdir)
        if fail == "pion_missing":
            (self.outdir / "Center_pion_Analysed_Dummy_Q4p4W2p74_lowe.root").unlink()
        ordinary_before = {p.name: owner.file_identity(p) for p in self.outdir.iterdir()
                           if p.suffix == ".root" or p.name in owner.names()}
        @contextmanager
        def owned(repo, sha, leaf):
            nonlocal worktree
            if leaf == "analysis":
                with tempfile.TemporaryDirectory(prefix="kaonlt-e8-4-canonical-five-analysis-", dir=repo.parent) as d:
                    worktree = Path(d) / "analysis"; worktree.mkdir()
                    yield worktree
                return
            self.assertEqual(leaf, "source")
            target = source
            target.mkdir()
            try:
                yield target
            finally:
                shutil.rmtree(source)
                if fail == "cleanup":
                    raise ValueError("cleanup_failed")
        def fake_git(repo, *args):
            nonlocal worktree
            if args[:2] == ("worktree", "add"):
                worktree = Path(args[3]); worktree.mkdir()
            elif args[:2] == ("worktree", "remove"):
                shutil.rmtree(Path(args[3]))
            elif args == ("rev-parse", "HEAD"):
                return SHA
            return ""
        def overlay_copy(repo, work, placement):
            copied = work.parent / "runtime-python/ltsep"
            copied.mkdir(parents=True)
            return {}, {"baseline_paths": baseline}
        def runtime_probes(*args):
            return [{"paths": {**baseline, "LTANAPATH": str(worktree),
                               "OUTPATH": str(worktree / "OUTPUT/Analysis/KaonLT")}}]
        @contextmanager
        def collection(clean):
            self.assertEqual(clean, source)
            self.assertFalse((clean / "OUTPUT").exists())
            module = Mock()
            def collect(**kw):
                self.assertEqual(kw["repo_root"], source)
                private = worktree.joinpath("OUTPUT/Analysis/KaonLT").resolve()
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
            private = tree.joinpath("OUTPUT/Analysis/KaonLT").resolve()
            for setting in ("Center", "Left", "Right"):
                self.assertEqual(len(open_frozen_pion_inputs(private, setting)), 2)
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
            mock(owner, "git", side_effect=fake_git)
            anchor = volatile.parent
            if aliases:
                anchor, raw = stack.enter_context(infrastructure_aliases(volatile, self.directory / "infrastructure"))
                baseline["VOLATILEPATH"] = str(raw)
            mock(owner, "VOLATILE_MOUNT", new=anchor)
            mock(owner.isolation, "path_fields", return_value=baseline)
            overlay = mock(owner, "prepare_runtime_overlay", side_effect=overlay_copy)
            probes = mock(owner.isolation, "probe_paths", side_effect=runtime_probes)
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

    def test_missing_pion_is_reported_before_copy_or_analysis_child(self):
        with self.flow("pion_missing") as (run, _, verify), patch.object(owner, "stage_copies") as stage:
            with self.assertRaisesRegex(ValueError,
                    "required_pion_input_missing:Center_pion_Analysed_Dummy_Q4p4W2p74_lowe.root"):
                owner.execute_gate(self.repo, self.outdir, self.output, SHA)
            stage.assert_not_called()
            run.assert_not_called()
            verify.assert_not_called()
        status = json.loads((self.outdir / "comparison-gate-status.json").read_text())
        self.assertEqual(status["status"], "failed")
        self.assertFalse(status.get("analysis_started", False))
        self.assertIn("required_pion_input_missing:", status["failure_reason"])
        self.assertFalse(self.output.exists())

    def test_owner_success_once_and_companions(self):
        with self.flow() as (run, preservation, verify):
            self.assertEqual(owner.execute_gate(self.repo, self.outdir, self.output, SHA), self.output)
        run.assert_called_once()
        verify.assert_called_once()
        self.assertEqual(preservation.call_count, 2)
        self.assertEqual(list(self.outdir.parents[2].glob("kaonlt-ab-analysis-*")), [])
        status = json.loads((self.transfer / "comparison-gate-status.json").read_text())
        self.assertEqual(status["status"], "success")
        self.assertEqual(status["diagnostic_outcomes"]["Left"]["status"], "available")
        self.assertTrue((self.transfer / "comparison.log").is_file())
        self.assertTrue((self.transfer / "comparison-run-summary.json").is_file())

    def test_each_failure_stops_delivery_preserves_receipt_and_checkout(self):
        for failure in ("source", "overlay", "paths", "external", "output_link", "copy", "child", "collector", "cleanup", "zip", "preservation", "original", "accepted", "installed"):
            with self.subTest(failure=failure):
                # Each attempt owns separate fresh destinations and target inventory.
                with tempfile.TemporaryDirectory(dir="/dev/shm") as fresh:
                    old = self.outdir, self.transfer, self.output
                    self.outdir = Path(fresh) / "expphy/volatile/user/OUTPUT/Analysis/KaonLT"; self.outdir.mkdir(parents=True)
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

    def test_group_placement_rejected_before_copy_or_child(self):
        with self.flow() as (run, _, _), patch.object(owner.isolation, "path_fields",
                return_value={"VOLATILEPATH": str(self.repo.parent)}), patch.object(owner, "stage_copies") as copy:
            with self.assertRaisesRegex(ValueError, "volatile_base_invalid"):
                owner.execute_gate(self.repo, self.outdir, self.output, SHA)
        copy.assert_not_called()
        run.assert_not_called()
        self.assertEqual(list(self.repo.parent.glob("kaonlt-ab-analysis-*")), [])

    def test_owner_alias_success_once_preserves_raw_resolver_and_source_collection(self):
        with self.flow(aliases=True) as (run, _, verify):
            owner.execute_gate(self.repo, self.outdir, self.output, SHA)
        run.assert_called_once(); verify.assert_called_once()
        summary = json.loads((self.transfer / "comparison-run-summary.json").read_text())
        placement = summary["runtime_placement"]
        self.assertEqual(placement["base"]["canonical"], str(self.outdir.parents[2]))
        self.assertNotEqual(placement["resolver_VOLATILEPATH"], placement["base"]["canonical"])
        self.assertTrue(placement["aliases"])
        self.assertFalse(Path(placement["parent"]["path"]).exists())

    def test_owner_mutable_alias_rejects_before_copy_or_child(self):
        with self.flow(aliases=True) as (run, _, _), patch.object(owner, "stage_copies") as copy:
            parent = owner.VOLATILE_MOUNT.parent
            parent.chmod(0o755)
            try:
                with self.assertRaisesRegex(ValueError, "alias_parent_writable"):
                    owner.execute_gate(self.repo, self.outdir, self.output, SHA)
            finally:
                parent.chmod(0o555)
        copy.assert_not_called(); run.assert_not_called()
        self.assertFalse((self.transfer / "comparison-gate-status.json").exists())


class VolatileRuntime(unittest.TestCase):
    """Real Linux filesystem/probe fixtures; Git and the analysis child are mocks."""

    def test_group_selector_reproduces_volatile_caller_failure_and_group_success(self):
        """Diagnostic only: a working caller does not validate split cleanup."""
        with self.fixture() as (repo, base, shared, package, config, calls):
            fields = owner.isolation.read_path_config(config.read_bytes())
            fields["ANALYSISPATH"] = str(repo.parent)
            config.write_text("".join(k + "=" + v + "\n" for k, v in fields.items()))
            pathing = package / "pathing.py"
            pathing.write_text('class InvalidPath(Exception): pass\n' + pathing.read_text().replace(
                'raise ValueError("no matching caller")', 'raise InvalidPath("PATHFILE not defined")'))
            original = {p: p.read_bytes() for p in package.rglob("*") if p.is_file()}
            previous = [base / ("kaonlt-ab-analysis-" + suffix) for suffix in ("c16r9oxs", "0o8v_67p")]
            for parent in previous:
                parent.mkdir(mode=0o700); (parent / "KEEP").write_bytes(b"previous failed attempt")
            real_glob = Path.glob
            def no_previous_scan(path, pattern):
                self.assertNotEqual(path, base, "must not enumerate prior attempts")
                return real_glob(path, pattern)
            old_text = subprocess.check_output(["git", "show", "f77b73786960f8a4136267372efa2b607f064706:testing/run_pion_t_template_comparison_gate.py"], cwd=ROOT).decode()
            selected = [n for n in ast.parse(old_text).body if isinstance(n, ast.FunctionDef) and n.name in {"analysis_worktree", "prepare_runtime_overlay"}]
            old = dict(vars(owner))
            exec(compile(ast.Module(body=selected, type_ignores=[]), "<starting-volatile-owner>", "exec"), old)
            baseline = owner.isolation.path_fields(repo, repo)
            record = owner.volatile_base(repo, baseline)
            with patch.object(Path, "glob", no_previous_scan):
                with old["analysis_worktree"](repo, SHA, record) as work:
                    record.update(baseline_paths=baseline, shared_output=str(shared))
                    env, overlay = old["prepare_runtime_overlay"](repo, work, record)
                    self.assertTrue(work.is_relative_to(base))
                    self.assertEqual(work.name, "analysis")
                    self.assertTrue(work.parent.name.startswith("kaonlt-ab-analysis-"))
                    with self.assertRaisesRegex(ValueError, "(?s)ltsep_probe_failed:.*PATHFILE not defined"):
                        owner.isolation.probe_paths(work, env, baseline)
                # The actual frozen helper and real create_ltsep_root/realpath
                # probes succeed under the historical group-sibling geometry.
                with patch.object(owner.isolation, "git", side_effect=owner.git):
                    with owner.isolation.owned_worktree(repo, SHA, "analysis") as work:
                        env, overlay = owner.isolation.prepare_runtime_overlay(repo, work)
                        probes = owner.isolation.probe_paths(work, env, baseline)
                        self.assertEqual(len(probes), 3)
                        self.assertTrue(work.is_relative_to(repo.parent))
                        self.assertNotEqual(work.stat().st_dev, base.stat().st_dev)
                        self.assertTrue(all(Path(p["paths"]["LTANAPATH"]) == work for p in probes))
                        self.assertTrue(all(p["paths"]["VOLATILEPATH"] == baseline["VOLATILEPATH"] for p in probes))
                        self.assertEqual(owner.isolation.read_path_config(
                            Path(overlay["selected_copied_config"]).read_bytes())["LTANAPATH"], str(work))
            self.assertEqual({p: p.read_bytes() for p in original}, original)
            for parent in previous:
                self.assertEqual((parent / "KEEP").read_bytes(), b"previous failed attempt")

    def test_frozen_group_parent_swap_refuses_cleanup_outside_owned_area(self):
        """Retained safety regression: reject the swap without any removal.

        Git is mocked exactly as in the existing Linux fixture. Only disposable
        fixture bytes can be removed; no real Git mutation or analysis runs.
        """
        with self.fixture() as (repo, base, shared, package, config, calls):
            unrelated = repo.parent / "unrelated-existing"
            (unrelated / "analysis").mkdir(parents=True)
            sentinel = unrelated / "analysis/KEEP"; sentinel.write_bytes(b"unrelated bytes")
            with patch.object(owner.isolation, "git", side_effect=owner.git):
                with self.assertRaisesRegex(ValueError, "worktree_cleanup_path_invalid"):
                    with owner.isolation.owned_worktree(repo, SHA, "analysis") as work:
                        parent = work.parent
                        saved = parent.with_name(parent.name + "-saved")
                        parent.rename(saved)
                        parent.symlink_to(unrelated)
                        # Reproduce guard rejection followed by context exit;
                        # neither cleanup route may touch the replacement.
                        with self.assertRaisesRegex(ValueError, "alias_or_invalid"):
                            owner.directory_identity(parent, private=True)
                        raise ValueError("owner rejected swapped code parent")
            removals = [args for _, args in calls if args[:2] == ("worktree", "remove")]
            self.assertTrue((saved / "analysis/science.py").is_file())
            self.assertEqual({"removals": len(removals), "sentinel_preserved": sentinel.exists()},
                             {"removals": 0, "sentinel_preserved": True},
                             "frozen cleanup issued removal through a swapped parent")

    def test_00_starting_owner_rejects_real_2700_and_leaves_parent(self):
        with self.fixture() as (repo, base, shared, package, config, calls):
            base.chmod(0o2700)
            text = subprocess.check_output(["git", "show",
                "e628f727d780a4d87c2ce928fedfcdc62958eb00:testing/run_pion_t_template_comparison_gate.py"], cwd=ROOT).decode()
            functions = [n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef)
                         and n.name in {"directory_identity", "analysis_worktree"}]
            namespace = dict(vars(owner))
            exec(compile(ast.Module(body=functions, type_ignores=[]), "<starting-mode-guard>", "exec"), namespace)
            record = owner.volatile_base(repo, owner.isolation.path_fields(repo, repo))
            with self.assertRaisesRegex(ValueError, "runtime_directory_owner_or_mode"):
                with namespace["analysis_worktree"](repo, SHA, record):
                    self.fail("old guard accepted inherited setgid")
            parents = list(base.glob("kaonlt-ab-analysis-*"))
            self.assertEqual(len(parents), 1)
            parent = parents[0]
            self.assertEqual(parent.stat().st_uid, os.geteuid())
            self.assertEqual(stat.S_IMODE(parent.stat().st_mode), 0o2700)
            self.assertEqual(owner.directory_identity(parent, private=True)["mode"], 0o2700)
            self.assertEqual(calls, [])
            parent.rmdir()  # Only this test's reproduction of the old leak.

    def test_private_modes_allow_only_0700_and_2700(self):
        with tempfile.TemporaryDirectory(dir="/dev/shm") as directory:
            path = Path(directory).resolve() / "private"; path.mkdir()
            denied = {0o770, 0o750, 0o755, 0o777, 0o1700, 0o4700, 0o3770}
            denied.update(0o700 | access | special for access in range(0o100)
                          for special in (0, 0o1000, 0o2000, 0o3000, 0o4000, 0o5000, 0o6000, 0o7000)
                          if (0o700 | access | special) not in (0o700, 0o2700))
            denied.update((0o600, 0o2600))
            for mode in (0o700, 0o2700):
                path.chmod(mode)
                self.assertEqual(stat.S_IMODE(path.stat().st_mode), mode)
                self.assertEqual(owner.directory_identity(path, private=True)["mode"], mode)
                with patch.object(os, "geteuid", return_value=os.geteuid()+1), self.assertRaisesRegex(ValueError, "owner_or_mode"):
                    owner.directory_identity(path, private=True)
            for mode in sorted(denied):
                with self.subTest(mode=oct(mode)):
                    path.chmod(mode)
                    self.assertEqual(stat.S_IMODE(path.stat().st_mode), mode)
                    with self.assertRaisesRegex(ValueError, "owner_or_mode"):
                        owner.directory_identity(path, private=True)
            path.chmod(0o700)

    def test_inherited_setgid_passes_actual_parent_and_layout_paths(self):
        for mode in (0o700, 0o2700):
            with self.subTest(mode=oct(mode)), self.fixture() as (repo, base, shared, package, config, calls):
                base.chmod(mode)
                baseline = owner.isolation.path_fields(repo, repo)
                record = owner.volatile_base(repo, baseline)
                with owner.analysis_worktree(repo, SHA, record) as work:
                    self.assertEqual(record["parent"]["mode"], mode)
                    self.assertEqual(owner.verify_runtime_parent(repo, record), Path(record["parent"]["path"]))
                    record.update(baseline_paths=baseline, shared_output=str(shared))
                    env, _ = owner.prepare_runtime_overlay(repo, work, record)
                    env = owner.private_environment(env, Path(record["parent"]["path"]))
                    paths = owner.isolation.probe_paths(work, env, baseline)[0]["paths"]
                    private = owner.private_output(work, paths, shared, Path(record["parent"]["path"]))
                    owner.verify_runtime_layout(repo, work, record, private, env)
                    roots = [Path(record["parent"]["path"]) / "ab-output"]
                    roots += [Path(env[k]) for k in ("TMPDIR", "XDG_CACHE_HOME", "MPLCONFIGDIR", "PYTHONPYCACHEPREFIX")]
                    for path in roots:
                        self.assertEqual(stat.S_IMODE(path.stat().st_mode), mode)
                        for denied in (0o770, 0o750, 0o755, 0o777, 0o1700, 0o4700, 0o3770):
                            path.chmod(denied)
                            try:
                                with self.subTest(path=path.name, mode=oct(denied)), self.assertRaisesRegex(ValueError, "owner_or_mode"):
                                    owner.runtime_directory(repo, record, path, private=True)
                                if path != work.parent / "runtime-python":
                                    with self.assertRaisesRegex(ValueError, "owner_or_mode"):
                                        owner.verify_runtime_layout(repo, work, record, private, env)
                            finally:
                                path.chmod(mode)
                    changed = roots[0]; changed.chmod(0o2700 if mode == 0o700 else 0o700)
                    try:
                        with self.assertRaisesRegex(ValueError, "location_identity_changed"):
                            owner.verify_runtime_layout(repo, work, record, private, env)
                    finally:
                        changed.chmod(mode)
                    Path(record["parent"]["path"]).chmod(0o2700 if mode == 0o700 else 0o700)
                    try:
                        with self.assertRaisesRegex(ValueError, "parent_identity_changed"):
                            owner.verify_runtime_parent(repo, record)
                    finally:
                        Path(record["parent"]["path"]).chmod(mode)
                self.assertFalse(work.parent.exists())

    def test_first_identity_failure_removes_only_new_empty_parent(self):
        for mode in (0o700, 0o2700):
            with self.subTest(mode=oct(mode)), self.fixture() as (repo, base, shared, package, config, calls):
                base.chmod(mode)
                artifacts(shared); before = owner.accepted_identities(shared)
                unrelated = base / "unrelated"; unrelated.mkdir(); (unrelated / "keep").write_bytes(b"keep")
                record = owner.volatile_base(repo, owner.isolation.path_fields(repo, repo))
                real_identity = owner.directory_identity
                def fail_first(path, private=False, owned=True):
                    if private and Path(path).name.startswith("kaonlt-ab-analysis-"):
                        raise ValueError("injected_first_identity_failure")
                    return real_identity(path, private=private, owned=owned)
                with patch.object(owner, "directory_identity", side_effect=fail_first), patch.object(shutil, "rmtree") as recursive:
                    with self.assertRaisesRegex(ValueError, "injected_first_identity_failure"):
                        with owner.analysis_worktree(repo, SHA, record):
                            self.fail("analysis reached after identity failure")
                    recursive.assert_not_called()
                self.assertEqual(list(base.glob("kaonlt-ab-analysis-*")), [])
                self.assertEqual((unrelated / "keep").read_bytes(), b"keep")
                self.assertEqual(calls, [])
                owner.verify_accepted(shared, before)

    def test_first_identity_failure_refuses_ambiguous_or_nonempty_cleanup(self):
        for defect in ("symlink", "replacement", "mode", "owner", "nonempty"):
            with self.subTest(defect=defect), self.fixture() as (repo, base, shared, package, config, calls):
                record = owner.volatile_base(repo, owner.isolation.path_fields(repo, repo))
                real_identity = owner.directory_identity
                created = []
                def fail_first(path, private=False, owned=True):
                    path = Path(path)
                    if private and path.name.startswith("kaonlt-ab-analysis-"):
                        created.append(path)
                        if defect in ("symlink", "replacement"):
                            path.rename(path.with_name(path.name + "-saved"))
                            if defect == "symlink":
                                path.symlink_to(repo.parent)
                            else:
                                path.mkdir(mode=0o700)
                        elif defect == "mode":
                            path.chmod(0o2700)
                        elif defect == "nonempty":
                            (path / "keep").write_bytes(b"keep")
                        raise ValueError("injected_first_identity_failure")
                    if defect == "owner" and created and path == created[0]:
                        with patch.object(os, "geteuid", return_value=os.geteuid()+1):
                            return real_identity(path, private=private, owned=owned)
                    return real_identity(path, private=private, owned=owned)
                with patch.object(owner, "directory_identity", side_effect=fail_first), patch.object(shutil, "rmtree") as recursive:
                    with self.assertRaisesRegex((ValueError, OSError), "alias_or_invalid|cleanup_unsafe|owner_or_mode|Directory not empty"):
                        with owner.analysis_worktree(repo, SHA, record):
                            self.fail("analysis reached after identity failure")
                    recursive.assert_not_called()
                self.assertTrue(os.path.lexists(created[0]))
                self.assertTrue((repo / "OUTPUT").is_symlink())
                self.assertEqual(calls, [])

    def test_00_pushed_owner_rejects_real_alias_chain(self):
        with self.fixture(aliases=True) as (repo, base, shared, package, config, calls):
            text = subprocess.check_output(["git", "show",
                "df07ecedf01b2d299c6fc61ac293c1fb4101ffd2:testing/run_pion_t_template_comparison_gate.py"], cwd=ROOT).decode()
            functions = [n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef)
                         and n.name in {"volatile_base", "directory_identity"}]
            namespace = {**vars(owner), "VOLATILE_MOUNT": owner.VOLATILE_MOUNT}
            exec(compile(ast.Module(body=functions, type_ignores=[]), "<pushed-volatile-guard>", "exec"), namespace)
            baseline = owner.isolation.path_fields(repo, repo)
            self.assertTrue(owner.VOLATILE_MOUNT.is_symlink())
            self.assertFalse(owner.VOLATILE_MOUNT.resolve(strict=True).is_mount())
            with self.assertRaisesRegex(ValueError, "volatile_mount_invalid"):
                namespace["volatile_base"](repo, baseline)
            record = owner.volatile_base(repo, baseline)
            self.assertEqual(record["base"]["path"], str(base))
            self.assertEqual(record["anchor"]["canonical"], str(base.parent))
            self.assertGreaterEqual(len(record["aliases"]), 3)
            self.assertEqual(calls, [])

    def test_user_writable_aliases_and_repository_device_are_rejected(self):
        with self.fixture(aliases=True) as (repo, base, shared, package, config, calls):
            baseline = owner.isolation.path_fields(repo, repo)
            anchor_parent = owner.VOLATILE_MOUNT.parent
            anchor_parent.chmod(0o755)
            try:
                with self.assertRaisesRegex(ValueError, "alias_parent_writable"):
                    owner.volatile_base(repo, baseline)
            finally:
                anchor_parent.chmod(0o555)
            for target in (repo, repo.parent, base, package.parent):
                link = base / "untrusted"; link.symlink_to(target)
                try:
                    with self.subTest(target=target), self.assertRaisesRegex(ValueError, "alias_parent_writable"):
                        owner.volatile_base(repo, {**baseline, "VOLATILEPATH": str(link)})
                finally:
                    link.unlink()
            with self.assertRaisesRegex(ValueError, "on_repository_device"):
                owner.volatile_base(package.parent, baseline)
            read_only = base / "read-only"; read_only.mkdir()
            (read_only / "alias").symlink_to(base)
            read_only.chmod(0o555)
            try:
                with self.assertRaisesRegex(ValueError, "untrusted_base_alias"):
                    owner.volatile_base(repo, {**baseline, "VOLATILEPATH": str(read_only / "alias")})
            finally:
                read_only.chmod(0o755)
            dangling = base / "dangling"; dangling.symlink_to(base / "missing")
            with self.assertRaises(FileNotFoundError):
                owner.volatile_base(repo, {**baseline, "VOLATILEPATH": str(dangling)})
            self.assertEqual(calls, [])
            self.assertEqual(list(base.glob("kaonlt-ab-analysis-*")), [])

    def test_trusted_anchor_alias_to_group_or_repository_is_rejected(self):
        with self.fixture(aliases=True) as (repo, base, shared, package, config, calls):
            outside = repo.parent / "user"; outside.mkdir()
            for target, raw in ((repo.parent, outside), (repo, repo / "user")):
                raw.mkdir(exist_ok=True)
                anchor = owner.VOLATILE_MOUNT
                anchor.parent.chmod(0o755)
                anchor.unlink(); anchor.symlink_to(target)
                anchor.parent.chmod(0o555)
                with self.subTest(target=target), self.assertRaisesRegex(ValueError, "base_overlap|on_repository_device"):
                    owner.volatile_base(repo, {"VOLATILEPATH": str(anchor / "user")})
            self.assertEqual(calls, [])

    def test_same_referent_alias_rebinding_stops_cleanup_without_traversal(self):
        with self.fixture(aliases=True) as (repo, base, shared, package, config, calls):
            artifacts(shared)
            before = owner.accepted_identities(shared)
            record = owner.volatile_base(repo, owner.isolation.path_fields(repo, repo))
            with self.assertRaisesRegex(ValueError, "base_identity_changed"):
                with owner.analysis_worktree(repo, SHA, record) as work:
                    anchor = owner.VOLATILE_MOUNT
                    anchor.parent.chmod(0o755)
                    anchor.unlink(); anchor.symlink_to(base.parent)
                    anchor.parent.chmod(0o555)
            self.assertTrue(Path(record["parent"]["path"]).is_dir())
            self.assertTrue(all(Path(args[3]).is_relative_to(repo.parent) for _, args in calls if args[:2] == ("worktree", "remove")))
            owner.verify_accepted(shared, before)

    @contextmanager
    def fixture(self, ambiguous=False, aliases=False):
        # Extract only the frozen harmless package fixture. Importing that test
        # module would also import plotting/science dependencies unnecessarily.
        source = ast.parse((ROOT / "testing/test_run_e8_4_fix5_canonical_five_plot_gate.py").read_text())
        builder = next(n for n in ast.walk(source) if isinstance(n, ast.FunctionDef) and n.name == "runtime_fixture")
        namespace = {"ROOT": ROOT, "Path": Path, "owner": owner.isolation}
        exec(compile(ast.Module(body=[builder], type_ignores=[]), "<installed-package-fixture>", "exec"), namespace)
        with tempfile.TemporaryDirectory(dir="/dev/shm") as directory, tempfile.TemporaryDirectory() as group, ExitStack() as stack:
            root = Path(directory).resolve()
            old_repo, old_work, package, config = namespace["runtime_fixture"](self, directory)
            repo = Path(group) / "repo"
            shutil.move(str(old_repo), repo)
            shutil.rmtree(old_work.parent)
            mount = root / "expphy/volatile"; base = mount / "user"; base.mkdir(parents=True)
            shared = base / "OUTPUT/Analysis/KaonLT"; shared.mkdir(parents=True)
            (repo / "OUTPUT").symlink_to(base / "OUTPUT")
            fields = owner.isolation.read_path_config(config.read_bytes())
            if aliases:
                mount, raw = stack.enter_context(infrastructure_aliases(base, root / "infrastructure"))
            else:
                raw = base
            fields.update(LTANAPATH=str(repo), VOLATILEPATH=str(raw), ANALYSISPATH=str(repo.parent))
            config.write_text("".join(k + "=" + v + "\n" for k, v in fields.items()))
            if ambiguous:
                config.with_name("duplicate.path").write_bytes(config.read_bytes())
            calls = []
            def fake_git(tree, *args):
                calls.append((Path(tree), args))
                if args[:2] == ("worktree", "add"):
                    work = Path(args[3]); work.mkdir()
                    shutil.copytree(repo / "farm_env", work / "farm_env")
                    (work / "science.py").write_text("# detached fixture\n")
                elif args[:2] == ("worktree", "remove"):
                    shutil.rmtree(Path(args[3]))
                elif args == ("rev-parse", "HEAD"):
                    return SHA
                return ""
            stack.enter_context(patch.object(owner, "git", side_effect=fake_git))
            stack.enter_context(patch.object(owner.isolation, "git", side_effect=fake_git))
            stack.enter_context(patch.object(owner, "VOLATILE_MOUNT", mount))
            stack.enter_context(patch.dict(os.environ, {"PYTHONPATH": str(package.parent),
                "PYTHONDONTWRITEBYTECODE": "1"}))
            yield repo, base, shared, package, config, calls

    def test_real_copied_overlay_and_all_writable_locations_are_volatile(self):
        with self.fixture(aliases=True) as (repo, base, shared, package, config, calls):
            original = {p: p.read_bytes() for p in package.rglob("*") if p.is_file()}
            populate_inputs(shared); artifacts(shared)
            accepted = owner.accepted_identities(shared)
            baseline = owner.isolation.path_fields(repo, repo)
            record = owner.volatile_base(repo, baseline)
            with owner.analysis_worktree(repo, SHA, record) as work:
                record.update(baseline_paths=baseline, shared_output=str(shared))
                env, overlay = owner.prepare_runtime_overlay(repo, work, record)
                env = owner.private_environment(env, Path(record["parent"]["path"]))
                probes = owner.isolation.probe_paths(work, env, baseline)
                private = owner.private_output(work, probes[0]["paths"], shared, Path(record["parent"]["path"]))
                owner.verify_runtime_layout(repo, work, record, private, env)
                copies = {}; plan = owner.copy_plan(shared, private)
                owner.stage_copies(shared, private, plan, copies)
                for p in (private,
                          *[Path(env[k]) for k in ("TMPDIR", "XDG_CACHE_HOME", "MPLCONFIGDIR", "PYTHONPYCACHEPREFIX")]):
                    self.assertTrue(p.is_relative_to(base))
                    self.assertFalse(p.is_relative_to(repo.parent))
                self.assertTrue(work.is_relative_to(repo.parent))
                self.assertTrue(Path(overlay["overlay_parent"]).is_relative_to(work.parent))
                self.assertNotEqual(work.stat().st_dev, private.stat().st_dev)
                self.assertEqual(record["parent"]["mode"], 0o700)
                for copy in copies.values():
                    self.assertTrue(Path(copy["runtime_container_before"]["resolved"]).is_relative_to(base))
                    self.assertNotEqual(copy["original"]["inode"], copy["runtime_container_before"]["inode"])
                patched = Path(overlay["selected_copied_config"]).read_bytes()
                self.assertEqual(patched, original[config].replace(str(repo).encode(), str(work).encode()))
                self.assertEqual({p: p.read_bytes() for p in original}, original)
                owner.isolation.verify_ltsep_preservation(overlay)
                owner.verify_accepted(shared, accepted); owner.verify_originals(copies)
                self.assertTrue(all(p["paths"]["VOLATILEPATH"] == baseline["VOLATILEPATH"] for p in probes))
                self.assertEqual(len(plan), 27)
                self.assertNotEqual(baseline["VOLATILEPATH"], str(base))
            self.assertFalse(work.parent.exists())
            self.assertTrue((repo / "OUTPUT").is_symlink())
            self.assertEqual(calls[0][1][:3], ("worktree", "add", "--detach"))
            self.assertEqual(calls[-1][1][:3], ("worktree", "remove", "--force"))

    def test_invalid_base_and_mount_fail_before_creation(self):
        with self.fixture() as (repo, base, shared, package, config, calls):
            alias = base.parent / "alias"; alias.symlink_to(repo.parent)
            for path in ("", "relative", str(repo.parent), str(alias), str(base / "missing"), str(base.parent)):
                with self.subTest(path=path), self.assertRaises((ValueError, OSError)):
                    owner.volatile_base(repo, {"VOLATILEPATH": path})
            self.assertFalse(base.parent.is_mount())
            with patch.object(owner, "VOLATILE_MOUNT", base / "missing"), self.assertRaises(FileNotFoundError):
                owner.volatile_base(repo, {"VOLATILEPATH": str(base)})
            with patch.object(os, "access", return_value=False), self.assertRaisesRegex(ValueError, "base_invalid"):
                owner.volatile_base(repo, {"VOLATILEPATH": str(base)})
            # Synthetic directories share a device: an alleged volatile mount
            # on the real group device must be refused even with free space.
            with patch.object(owner, "GROUP_MOUNT", base.parent / ".."), self.assertRaisesRegex(ValueError, "base_overlap"):
                owner.volatile_base(repo, {"VOLATILEPATH": str(base)})
            with patch.object(owner, "GROUP_MOUNT", package.parent), self.assertRaisesRegex(ValueError, "on_group_device"):
                owner.volatile_base(repo, {"VOLATILEPATH": str(base)})
            with patch.object(os, "geteuid", return_value=os.geteuid()+1), self.assertRaisesRegex(ValueError, "owner_or_mode"):
                owner.volatile_base(repo, {"VOLATILEPATH": str(base)})
            self.assertEqual(calls, [])
            self.assertEqual(list(base.glob("kaonlt-ab-analysis-*")), [])

    def test_injected_body_and_overlay_failures_clean_owned_parent(self):
        for defect in ("body", "ambiguous"):
            with self.subTest(defect=defect), self.fixture(ambiguous=defect == "ambiguous") as (repo, base, shared, package, config, calls):
                baseline = owner.isolation.path_fields(repo, repo)
                record = owner.volatile_base(repo, baseline)
                with self.assertRaisesRegex(ValueError, "injected|matching_config_not_unique"):
                    with owner.analysis_worktree(repo, SHA, record) as work:
                        record.update(baseline_paths=baseline, shared_output=str(shared))
                        owner.prepare_runtime_overlay(repo, work, record)
                        raise ValueError("injected")
                self.assertFalse(work.parent.exists())
                self.assertTrue(config.exists())

    def test_swapped_parent_refuses_cleanup_and_never_touches_target(self):
        with self.fixture() as (repo, base, shared, package, config, calls):
            record = owner.volatile_base(repo, owner.isolation.path_fields(repo, repo))
            with self.assertRaisesRegex(ValueError, "alias_or_invalid"):
                with owner.analysis_worktree(repo, SHA, record) as work:
                    parent = Path(record["parent"]["path"])
                    saved = parent.with_name(parent.name + "-saved")
                    parent.rename(saved)
                    parent.symlink_to(repo.parent)
            self.assertTrue(saved.is_dir())
            self.assertTrue((repo / "OUTPUT").is_symlink())
            self.assertTrue(all(Path(args[3]).is_relative_to(repo.parent) for _, args in calls if args[:2] == ("worktree", "remove")))

    def test_swapped_cache_escape_and_group_output_are_rejected(self):
        with self.fixture() as (repo, base, shared, package, config, calls):
            baseline = owner.isolation.path_fields(repo, repo)
            record = owner.volatile_base(repo, baseline)
            with owner.analysis_worktree(repo, SHA, record) as work:
                record.update(baseline_paths=baseline, shared_output=str(shared))
                env, _ = owner.prepare_runtime_overlay(repo, work, record)
                env = owner.private_environment(env, Path(record["parent"]["path"]))
                paths = owner.isolation.probe_paths(work, env, baseline)[0]["paths"]
                private = owner.private_output(work, paths, shared, Path(record["parent"]["path"]))
                owner.verify_runtime_layout(repo, work, record, private, env)
                with self.assertRaisesRegex(ValueError, "private_output_path_changed"):
                    owner.verify_runtime_layout(repo, work, record, repo.parent / "ab-output", env)
                with self.assertRaisesRegex(ValueError, "outside_volatile_parent"):
                    owner.verify_runtime_layout(repo, work, record, private, {**env, "TMPDIR": str(repo.parent)})
                cache = Path(env["TMPDIR"]); saved = cache.with_name("saved-tmp")
                cache.rename(saved); cache.symlink_to(repo.parent)
                try:
                    with self.assertRaisesRegex(ValueError, "alias_or_invalid"):
                        owner.verify_runtime_layout(repo, work, record, private, env)
                finally:
                    cache.unlink(); saved.rename(cache)  # Fixture repairs only.
                link = work / "OUTPUT"; link.unlink(); link.symlink_to(shared)
                try:
                    with self.assertRaisesRegex(ValueError, "private_link_changed"):
                        owner.verify_runtime_layout(repo, work, record, private, env)
                finally:
                    link.unlink(); link.symlink_to(Path(record["parent"]["path"]) / "ab-output")

    def test_cleanup_capability_and_changed_base_reject_before_worktree(self):
        with self.fixture() as (repo, base, shared, package, config, calls):
            record = owner.volatile_base(repo, owner.isolation.path_fields(repo, repo))
            with patch.object(shutil.rmtree, "avoids_symlink_attacks", False), self.assertRaisesRegex(ValueError, "cleanup_not_symlink_safe"):
                with owner.analysis_worktree(repo, SHA, record):
                    self.fail("unsafe cleanup accepted")
            saved = base.with_name("old-base"); base.rename(saved); base.mkdir()
            with self.assertRaisesRegex(ValueError, "base_identity_changed"):
                with owner.analysis_worktree(repo, SHA, record):
                    self.fail("replaced base accepted")
            self.assertEqual(calls, [])

    def test_changed_worktree_identity_and_parent_mode_refuse_cleanup(self):
        for defect in ("worktree", "mode"):
            with self.subTest(defect=defect), self.fixture() as (repo, base, shared, package, config, calls):
                record = owner.volatile_base(repo, owner.isolation.path_fields(repo, repo))
                with self.assertRaisesRegex(ValueError, "identity_changed|owner_or_mode"):
                    with owner.analysis_worktree(repo, SHA, record) as work:
                        if defect == "worktree":
                            work.rename(work.with_name("saved-analysis")); work.mkdir()
                        else:
                            Path(record["parent"]["path"]).chmod(0o755)
                self.assertTrue(Path(record["parent"]["path"]).exists())
                self.assertTrue((repo / "OUTPUT").is_symlink())
                if defect == "worktree":
                    self.assertFalse(any(args[:2] == ("worktree", "remove") for _, args in calls))
                else:
                    self.assertFalse(work.parent.exists())  # Safe code role cleaned independently.


def populate_inputs(root):
    for eps in ("lowe", "highe"):
        for s in ("Center", "Left", "Right"):
            filenames = [f"{s}_kaon_Analysed_{r}_Q4p4W2p74_{eps}.root" for r in ("Data", "Dummy")]
            if eps == "lowe":
                filenames += [f"{s}_pion_Analysed_{r}_Q4p4W2p74_{eps}.root" for r in ("Data", "Dummy")]
            filenames += [f"Prod_Coin_Q4p4W2p74{s.lower()}_{eps}{suffix}"
                          for suffix in ((".root", ".hist") if eps == "lowe" else (".root",))]
            for name in filenames:
                p = root / name
                if not p.exists():
                    p.write_bytes(("fixture:" + name).encode())


def open_frozen_pion_inputs(root, setting):
    """Exercise production filename/first-use semantics without importing ROOT."""
    tree = ast.parse((ROOT / "src/cuts/rand_sub.py").read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and
                    n.name == "_open_subtracted_particle_tree_bundle")
    opened = []
    def open_file(path):
        opened.append((Path(path), Path(path).read_bytes()))
        return Mock()
    ns = {"os": os, "open_root_file": open_file,
          "_print_rand_debug": lambda *a, **kw: None,
          "_debug_tree_status": lambda *a: None,
          "get_prompt_tree_name": lambda *a, **kw: "fixture_prompt",
          "get_rand_tree_name": lambda *a, **kw: "fixture_rand"}
    exec(compile(ast.Module(body=[function], type_ignores=[]), "<frozen-pion-opener>", "exec"), ns)
    ns[function.name](str(root), setting, "pion", "Analysed_Data_Q4p4W2p74_lowe",
                      "Analysed_Dummy_Q4p4W2p74_lowe", "low")
    return opened


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
        return owner.private_output(self.work, self.paths, self.shared, self.base)

    def test_private_writes_keep_original_bytes_links_and_independent_inodes(self):
        populate_inputs(self.shared)
        artifacts(self.shared)
        (self.shared / "ordinary.pdf").write_bytes(b"ordinary science")
        (self.shared / "accepted-alias").symlink_to(self.shared / "ordinary.pdf")
        before = {p.name: owner.file_identity(p) for p in self.shared.iterdir()}
        accepted = owner.accepted_identities(self.shared)
        private = self.private()
        plan = owner.copy_plan(self.shared, private)
        self.assertEqual(len(plan), 27)
        self.assertEqual(sum(p["bytes"] for p in plan), 1332)
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

    def test_pion_pairs_reach_frozen_first_use_with_private_bytes(self):
        populate_inputs(self.shared)
        private = self.private()
        copies = {}
        owner.stage_copies(self.shared, private, owner.copy_plan(self.shared, private), copies)
        for setting in ("Center", "Left", "Right"):
            for path, raw in open_frozen_pion_inputs(private, setting):
                self.assertEqual(raw, (self.shared / path.name).read_bytes())
                self.assertEqual(path.parent, private)
                self.assertIn(path.name, copies)
                self.assertNotEqual(path.stat().st_ino, (self.shared / path.name).stat().st_ino)
        self.assertFalse(list(private.glob("*pion*highe*")))
        owner.verify_originals(copies)

    def test_old_inventory_reproduces_exact_center_first_use_failure(self):
        populate_inputs(self.shared)
        private = self.private()
        source = subprocess.check_output(["git", "show",
            "b5769e51c01655d9756481c91b258f029e219211:testing/run_pion_t_template_comparison_gate.py"], cwd=ROOT)
        function = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == "copy_plan")
        ns = dict(vars(owner))
        exec(compile(ast.Module(body=[function], type_ignores=[]), "<old-copy-plan>", "exec"), ns)
        old_plan = ns["copy_plan"](self.shared, private)
        self.assertFalse(any("_pion_" in p["source"] for p in old_plan))
        for entry in old_plan:
            shutil.copyfile(entry["source"], entry["target"])
        with self.assertRaisesRegex(FileNotFoundError, "Center_pion_Analysed_Data_Q4p4W2p74_lowe.root"):
            open_frozen_pion_inputs(private, "Center")
        self.assertFalse((private / "Center_pion_Analysed_Dummy_Q4p4W2p74_lowe.root").exists())

    def test_missing_consumed_pion_pair_fails_before_any_staging(self):
        populate_inputs(self.shared)
        private = self.private()
        for setting in ("Center", "Left", "Right"):
            pair = [self.shared / f"{setting}_pion_Analysed_{r}_Q4p4W2p74_lowe.root" for r in ("Data", "Dummy")]
            for missing in ([pair[0]], [pair[1]], pair):
                raw = {p: p.read_bytes() for p in missing}
                for p in missing: p.unlink()
                with self.subTest(setting=setting, missing=missing):
                    with self.assertRaisesRegex(ValueError, "required_pion_input_missing:" + missing[0].name):
                        owner.copy_plan(self.shared, private)
                    self.assertEqual(list(private.iterdir()), [])
                for p, value in raw.items(): p.write_bytes(value)

    def test_invalid_pion_file_and_plan_change_fail_closed(self):
        populate_inputs(self.shared)
        private = self.private()
        p = self.shared / "Center_pion_Analysed_Dummy_Q4p4W2p74_lowe.root"
        raw = p.read_bytes()
        for kind in ("directory", "dangling"):
            p.unlink()
            if kind == "directory": p.mkdir()
            else: p.symlink_to(self.base / "missing-pion")
            with self.subTest(kind=kind), self.assertRaisesRegex(ValueError, "required_pion_input_missing:" + p.name):
                owner.copy_plan(self.shared, private)
            if kind == "directory": p.rmdir()
            else: p.unlink()
            p.write_bytes(raw)
        plan = owner.copy_plan(self.shared, private)
        p.unlink()
        with self.assertRaisesRegex(ValueError, "required_pion_input_missing:" + p.name):
            owner.stage_copies(self.shared, private, plan, {})
        self.assertEqual(list(private.iterdir()), [])

    def test_unavailable_right_needs_no_pion_and_high_prepass_needs_no_pion(self):
        populate_inputs(self.shared)
        (self.shared / "Prod_Coin_Q4p4W2p74right_lowe.root").unlink()
        for p in self.shared.glob("Right_pion_*.root"): p.unlink()
        # No high-epsilon pion files are present, including available high settings.
        private = self.private()
        plan = owner.copy_plan(self.shared, private)
        self.assertEqual(sum("_pion_" in p["source"] for p in plan), 4)
        owner.stage_copies(self.shared, private, plan, {})
        self.assertEqual(owner.diamond_choices(self.shared), owner.diamond_choices(private))
        for setting in ("Center", "Left"):
            self.assertEqual(len(open_frozen_pion_inputs(private, setting)), 2)

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
            owner.private_output(self.work, bad, self.shared, self.base)
        alias = self.base / "alias"; alias.symlink_to(self.work)
        with self.assertRaisesRegex(ValueError, "worktree_invalid"):
            owner.private_output(alias, self.paths, self.shared, self.base)

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
        private = owner.private_output(self.work, paths, alias / "OUTPUT/Analysis/KaonLT", self.base)
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
        env = owner.private_environment({"VOLATILEPATH": "stable"}, self.base)
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
