"""Synthetic isolated-runtime gate checks: no farm, ROOT or real launcher."""
from contextlib import contextmanager, ExitStack, redirect_stdout, redirect_stderr
from copy import deepcopy
import io
import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
import tempfile
import time
from types import SimpleNamespace
import unittest
from unittest.mock import patch, Mock
import zipfile

from testing import run_e8_4_fix5_canonical_five_plot_gate as owner
from testing.test_run_e8_4_fix5_left_lowe_plot_gate import pages as accepted_pages

ROOT = Path(__file__).resolve().parents[1]
SHA = "a" * 40


def pages(phi="Left", epsilon="lowe"):
    payload = accepted_pages()
    payload["setting"].update(phi_setting=phi, epsilon_filename_token=epsilon,
                              epsilon_setting={"lowe": "low", "highe": "high"}[epsilon])
    payload["pdf_basename"] = owner.pdf_name(phi, epsilon)
    return payload


class OwnerTests(unittest.TestCase):
    def runtime_fixture(self, directory):
        root = Path(directory).resolve()
        repo = root / "ordinary"; repo.mkdir()
        parent = root / "kaonlt-owner"; parent.mkdir()
        worktree = parent / "analysis"; worktree.mkdir()
        for tree in (repo, worktree):
            (tree / "farm_env").mkdir()
            for name in ("ltsep_paths.py", "print_ltsep_path_fields.py"):
                (tree / "farm_env" / name).write_bytes((ROOT / "farm_env" / name).read_bytes())
        package = root / "installed/ltsep"; package.mkdir(parents=True)
        (package / "PATH_TO_DIR").mkdir()
        (package / "__init__.py").write_text(
            'from .pathing import SetPath\n'
            'class Root:\n'
            ' def __init__(self, caller, *args):\n'
            '  p=SetPath(caller)\n'
            '  for k in p.values: setattr(self,k,p.getPath(k))\n'
            '  self.OUTPATH = str(__import__("pathlib").Path(self.LTANAPATH)/"OUTPUT/Analysis"/(self.ANATYPE+"LT"))\n', encoding="utf-8")
        (package / "pathing.py").write_text(
            'from pathlib import Path\n'
            'class SetPath:\n'
            ' def __init__(self, caller):\n'
            '  self.values={}\n'
            '  for f in (Path(__file__).parent/"PATH_TO_DIR").glob("*.path"):\n'
            '   d=dict(l.split("=",1) for l in f.read_text().splitlines())\n'
            '   if Path(caller).resolve().is_relative_to(Path(d["ANALYSISPATH"])):\n'
            '    self.values=d\n'
            '  if not self.values: raise ValueError("no matching caller")\n'
            ' def getPath(self,key):\n'
            '  return self.values[key].replace("${USER}",self.values["USER"])\n', encoding="utf-8")
        fields = {k: str(root / (k.lower() + "-fixture")) for k in owner.PATH_FIELDS}
        fields.update(LTANAPATH=str(repo), ANALYSISPATH=str(root), USER="resolved-user",
                      HOST="synthetic-host", ANATYPE="Kaon")
        config = package / "PATH_TO_DIR/user.path"
        config.write_text("".join(k + "=" + v + "\n" for k, v in fields.items()), encoding="utf-8")
        return repo, worktree, package, config

    def local_python(self, command, **kwargs):
        # Execute the actual probes with the workstation interpreter; change
        # only the farm executable name, never their imports or returned fields.
        command = list(command)
        if command[0] == "python3": command[0] = sys.executable
        return self.real_subprocess_run(command, **kwargs)

    def test_real_fixed_config_overlay_imports_and_three_callers(self):
        with tempfile.TemporaryDirectory() as directory:
            repo, worktree, package, config = self.runtime_fixture(directory)
            before = {str(p): p.read_bytes() for p in package.rglob("*") if p.is_file()}
            ambient = {**os.environ, "PYTHONPATH": str(package.parent),
                "LT_ANALYSIS_DEBUG_LEFT_LOW": "1", "LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE": "bad.npz",
                "LT_ANALYSIS_ALLOW_UNPAIRED_CANONICAL_BINNING": "1", "LT_BG_SIGMA0_LEFT_LOW_ROOT": "external.root"}
            self.real_subprocess_run = subprocess.run
            with patch.dict(os.environ, ambient, clear=True), patch.object(owner.subprocess, "run", side_effect=self.local_python):
                baseline = owner.path_fields(repo, repo)
                self.assertEqual(Path(baseline["OUTPATH"]), repo / "OUTPUT/Analysis/KaonLT")
                # A changed caller alone still reads the ordinary fixed config.
                self.assertEqual(owner.path_fields(worktree, worktree)["LTANAPATH"], str(repo))
                env, record = owner.prepare_runtime_overlay(repo, worktree)
                for key in ("LT_ANALYSIS_DEBUG_LEFT_LOW", "LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE",
                            "LT_ANALYSIS_ALLOW_UNPAIRED_CANONICAL_BINNING"):
                    self.assertNotIn(key, env)
                self.assertEqual(env["LT_BG_SIGMA0_LEFT_LOW_ROOT"], "external.root")
                probes = owner.probe_paths(worktree, env, baseline)
                self.assertEqual(len(probes), 3)
                self.assertTrue(all(p["paths"]["LTANAPATH"] == str(worktree) for p in probes))
                self.assertTrue(all(p["paths"] == probes[0]["paths"] for p in probes))
                self.assertEqual(Path(probes[0]["paths"]["OUTPATH"]), worktree / "OUTPUT/Analysis/KaonLT")
                self.assertNotEqual(probes[0]["paths"]["OUTPATH"], baseline["OUTPATH"])
                for key in set(owner.PATH_FIELDS) - {"OUTPATH", "LTANAPATH"}:
                    self.assertEqual(probes[0]["paths"][key], baseline[key])
                child = owner.import_identity(worktree, env)
                self.assertEqual(Path(child["package_file"]),
                    Path(record["overlay_parent"]) / "ltsep/__init__.py")
                patched = Path(record["selected_copied_config"]).read_bytes()
                original = before[str(config)]
                self.assertEqual(patched, original.replace(b"LTANAPATH=" + str(repo).encode(),
                                                          b"LTANAPATH=" + str(worktree).encode()))
                self.assertNotEqual(record["original_copied_file_sha256"], record["patched_copied_file_sha256"])
                self.assertTrue(owner.verify_ltsep_preservation(record)["passed"])
                self.assertEqual({str(p): p.read_bytes() for p in package.rglob("*") if p.is_file()}, before)
                # Exercise the unchanged run function's environment propagation.
                process = Mock(stdout=io.BytesIO(b"synthetic\n")); process.wait.return_value = 0
                with patch.object(owner.subprocess, "Popen", return_value=process) as child_run, redirect_stderr(io.StringIO()):
                    owner.run_analysis(worktree, Path(directory) / "child.log", env)
                self.assertEqual(child_run.call_args.kwargs["env"], env)
                config.write_bytes(original + b"unexpected")
                with self.assertRaisesRegex(ValueError, "installed_ltsep_changed"):
                    owner.verify_ltsep_preservation(record)
                self.assertEqual(config.read_bytes(), original + b"unexpected")  # Never restored.

    def test_real_overlay_rejects_ambiguous_configuration(self):
        with tempfile.TemporaryDirectory() as directory:
            repo, worktree, package, config = self.runtime_fixture(directory)
            config.with_name("duplicate.path").write_bytes(config.read_bytes())
            self.real_subprocess_run = subprocess.run
            with patch.dict(os.environ, {"PYTHONPATH": str(package.parent)}), patch.object(owner.subprocess, "run", side_effect=self.local_python):
                with self.assertRaisesRegex(ValueError, "matching_config_not_unique"):
                    owner.prepare_runtime_overlay(repo, worktree)

    def test_plot_ltsep_gate_rejects_stable_field_drift_and_wrong_outpath(self):
        with tempfile.TemporaryDirectory() as directory:
            repo, worktree, package, config = self.runtime_fixture(directory)
            self.real_subprocess_run = subprocess.run
            with patch.dict(os.environ, {"PYTHONPATH": str(package.parent)}), patch.object(owner.subprocess, "run", side_effect=self.local_python):
                env, record = owner.prepare_runtime_overlay(repo, worktree)
                baseline = record["baseline_paths"]
                actual = owner.probe_paths(worktree, env, baseline)[0]["paths"]
            for key in set(owner.PATH_FIELDS) - {"LTANAPATH", "OUTPATH"}:
                with self.subTest(key=key), patch.object(owner, "path_fields", return_value={**actual, key: "wrong"}):
                    with self.assertRaisesRegex(ValueError, "ltsep_stable_fields_changed"):
                        owner.probe_paths(worktree, env, baseline)
            for wrong in (baseline["OUTPATH"], "relative/output"):
                with self.subTest(outpath=wrong), patch.object(owner, "path_fields", return_value={**actual, "OUTPATH": wrong}):
                    with self.assertRaisesRegex(ValueError, "ltsep_plot_ltsep_outpath_invalid"):
                        owner.probe_paths(worktree, env, baseline)

    def test_background_resolution_sources_config_with_ltsep_user(self):
        with tempfile.TemporaryDirectory() as directory:
            tree = Path(directory).resolve()
            (tree / "background_samples").mkdir()
            (tree / "background_samples/background_samples.conf").write_text(
                f'BACKGROUND_SIMC_PATH="{tree.as_posix()}/nonexistent/${{USER}}/background-test"\n', encoding="utf-8")
            simc = tree / "simc"
            paths = {"SIMCPATH": str(simc), "VOLATILEPATH": str(tree), "USER": "ltsep-user"}
            real_run = subprocess.run
            def shell(command, **kwargs):
                self.assertEqual(kwargs["env"]["USER"], "ltsep-user")
                command = list(command)
                if os.name == "nt":
                    command[0] = str(Path(shutil.which("git")).resolve().parents[1] / "bin/bash.exe")
                    command[-1] = Path(command[-1]).as_posix()
                return real_run(command, **kwargs)
            with patch.dict(os.environ, {"USER": "unrelated-user"}), patch.object(owner.subprocess, "run", side_effect=shell), patch.object(Path, "is_dir", lambda p: p == simc), patch.object(Path, "is_symlink", return_value=True), patch.object(Path, "exists", return_value=True), patch.object(owner.os, "readlink", return_value="existing-target"):
                result = owner.external_symlink_preflight(tree, paths)
            self.assertEqual(result["background_simc_path"], f"{tree.as_posix()}/nonexistent/ltsep-user/background-test")

    def test_first_atomic_status_publication_is_canonical(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            original_link = os.link
            first = []
            def publish(source, destination):
                first.append(json.loads(Path(source).read_text()))
                original_link(source, destination)
            with patch.object(owner.accepted.os, "link", side_effect=publish):
                status = owner.GateStatus(root, root / "fresh.zip", SHA)
            self.assertEqual(len(first), 1)
            self.assertEqual(first[0]["schema_version"], "e8_4_fix5_canonical_five_owner_gate_status/v1")
            self.assertEqual(first[0]["settings"], list(owner.CANONICAL_SETTINGS))
            self.assertNotIn("phi", first[0]); self.assertNotIn("epsilon", first[0])
            self.assertEqual(json.loads(status.path.read_text()), first[0])

    def test_launcher_command_and_detached_cwd(self):
        with tempfile.TemporaryDirectory() as directory:
            worktree = Path(directory) / "analysis"
            process = Mock(stdout=io.BytesIO(b"synthetic\n"))
            process.wait.return_value = 0
            with patch.object(owner.subprocess, "Popen", return_value=process) as launch, redirect_stderr(io.StringIO()):
                self.assertEqual(owner.run_analysis(worktree, Path(directory) / "analysis.log"), 0)
            self.assertEqual(launch.call_args.args[0], ["./run_Prod_Analysis.sh", "4p4", "2p74"])
            self.assertEqual(launch.call_args.kwargs["cwd"], worktree)

    def test_path_probes_require_exact_analysis_root_in_all_caller_contexts(self):
        worktree = Path(tempfile.gettempdir()).resolve() / "isolated-analysis"
        callers = (worktree, worktree / "src/setup/set_sig_fortran.py", worktree / "set_SymLinks.sh")
        for wrong in (None, 0, 1, 2):
            calls = []
            def run(command, **kwargs):
                index = len(calls); calls.append(command)
                fields = {name: "value" for name in owner.PATH_FIELDS}
                fields["LTANAPATH"] = str(ROOT if index == wrong else worktree)
                return SimpleNamespace(returncode=0, stdout=",".join(fields[name] for name in owner.PATH_FIELDS), stderr="")
            with self.subTest(wrong=wrong), patch.object(owner.subprocess, "run", side_effect=run):
                if wrong is None:
                    self.assertEqual(len(owner.probe_paths(worktree)), 3)
                else:
                    with self.assertRaisesRegex(ValueError, "ltanapath_not_isolated"):
                        owner.probe_paths(worktree)
            self.assertEqual([c[-1] for c in calls], [str(c) for c in callers[:len(calls)]])

    def test_external_links_fail_closed_without_repairs(self):
        base = Path(tempfile.gettempdir()).resolve()
        simc, background, volatile = base / "simc-fixture", base / "background-fixture", base / "volatile-fixture"
        paths = {"SIMCPATH": str(simc), "VOLATILEPATH": str(volatile), "USER": "resolved-user"}
        links = {simc / leaf: "existing-target" for leaf in ("OUTPUTS", "input", "worksim")}
        links[background / "worksim"] = str(volatile) + "/worksim/"
        for defect in (None, "OUTPUTS", "input", "worksim", "broken", "background"):
            def is_link(path):
                return path in links and path != simc / str(defect)
            def exists(path):
                return not (defect == "broken" and path == simc / "input")
            def readlink(path):
                return "wrong" if defect == "background" and Path(path) == background / "worksim" else links[Path(path)]
            with self.subTest(defect=defect), patch.object(Path, "is_symlink", is_link), patch.object(Path, "exists", exists), patch.object(Path, "is_dir", return_value=True), patch.object(owner.os, "readlink", side_effect=readlink), patch.object(owner.subprocess, "run", return_value=SimpleNamespace(returncode=0, stdout=str(background), stderr="")), patch.object(owner.os, "symlink") as repair, patch.object(Path, "unlink") as remove:
                if defect is None:
                    self.assertTrue(owner.external_symlink_preflight(ROOT, paths)["passed"])
                else:
                    with self.assertRaisesRegex(ValueError, "external_simc_link|background_simc_worksim"):
                        owner.external_symlink_preflight(ROOT, paths)
                repair.assert_not_called(); remove.assert_not_called()

    def test_worktree_identity_and_cleanup_bounds(self):
        # Synthetic ordinary repository; owned sibling directories never use
        # the real workstation checkout's parent.
        fixture = tempfile.TemporaryDirectory()
        self.addCleanup(fixture.cleanup)
        ordinary = Path(fixture.name) / "ordinary"
        ordinary.mkdir()
        for defect in (None, "head", "branch", "status"):
            calls = []
            def git(repo, *args):
                calls.append((repo, args))
                if args[0] == "rev-parse": return "b" * 40 if defect == "head" else SHA
                if args[0] == "branch": return "test" if defect == "branch" else ""
                if args[0] == "status": return " M src/main.py" if defect == "status" else ""
                return ""
            with self.subTest(defect=defect), patch.object(owner, "git", side_effect=git):
                if defect is None:
                    with owner.owned_worktree(ordinary, SHA, "analysis") as worktree:
                        self.assertEqual(worktree.name, "analysis")
                        self.assertEqual(worktree.parent.parent, ordinary.parent)
                else:
                    with self.assertRaisesRegex(ValueError, "worktree_"), owner.owned_worktree(ordinary, SHA, "analysis"):
                        self.fail("invalid worktree yielded")
            self.assertEqual(calls[0][1][:3], ("worktree", "add", "--detach"))
            self.assertEqual(calls[-1], (ordinary, ("worktree", "remove", "--force", calls[0][1][-2])))
        for parent, tree in ((ROOT, ROOT / "analysis"), (ROOT / "nested", ROOT / "nested/analysis"),
                             (ROOT.parent, ROOT)):
            with self.assertRaisesRegex(ValueError, "cleanup_path_invalid"):
                owner.cleanup_path(ROOT, parent, tree, "analysis")

    def test_page_gates_all_five_and_negative_cases(self):
        for setting in owner.CANONICAL_SETTINGS:
            payload = pages(**setting)
            original = deepcopy(payload)
            self.assertEqual(owner.verify_pages(payload, **setting), len(payload["pages"]))
            self.assertEqual(payload, original)
        for change, reason in (
            (lambda p: p["setting"].update(phi_setting="Right"), "setting_invalid"),
            (lambda p: p.update(renderer_failures=["failure"]), "renderer_failures"),
            (lambda p: p["pages"].pop(), "new_page_missing"),
            (lambda p: p["pages"][-2].update(represented_phi_inventory=[]), "child_inventory"),
            (lambda p: p["pages"][-2].update(t_index=8), "t_identity"),
            (lambda p: p["pages"].pop(0), "critical_page")):
            payload = pages(); change(payload)
            with self.subTest(reason=reason), self.assertRaisesRegex(ValueError, reason):
                owner.verify_pages(payload, "Left", "lowe")
        with self.assertRaisesRegex(ValueError, "not_canonical"):
            owner.verify_pages(pages(), "Right", "lowe")

    def exercise(self, directory, defect=None, cli=False, lineage_override=None, preflight_only=False):
        root = Path(directory).resolve()
        repo = root / "ordinary"; repo.mkdir()
        outdir = root / "volatile/OUTPUT/Analysis/KaonLT"; outdir.mkdir(parents=True)
        output = root / "fresh.zip"
        companions = {"log": root / "fresh.log", "run_summary": root / "fresh-run-summary.json",
                      "gate_status": root / "fresh-gate-status.json"}
        if defect and defect.startswith("existing_"):
            companions[defect.removeprefix("existing_")].write_bytes(b"earlier attempt")
        (repo / "testing").mkdir()
        (repo / owner.PROFILE).write_bytes((ROOT / owner.PROFILE).read_bytes())
        for name in owner.FARM_OUTPUTS:
            path = repo / name; path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"accepted farm-local bytes")
        primary_bytes = {name: (repo / name).read_bytes() for name in owner.FARM_OUTPUTS}
        dirty = " M src/models/xmodel_kaon_pl.f\n?? src/kaon/functions/Q4p4W2p74.model"
        profile = owner.resolved_profile(repo, SHA)
        candidates = {}
        for name in owner.CANDIDATES:
            (outdir / name).write_text('{"synthetic": true}')
            candidates[name] = owner.sha256(outdir / name)
        events, calls = [], []
        ran = False
        def materialization(directory):
            if defect == "materialization": raise ValueError("synthetic_materialization_failure")
            events.append("materialization_verify")
            return {"passed": True, "files": {}}
        def staging(*args, **kwargs):
            if defect == "staging": raise ValueError("synthetic_staging_failure")
            events.append("candidate_staging")
            return {"passed": True, "files": []}
        def lineage(*args):
            events.append("lineage_preflight")
            if lineage_override is not None: return lineage_override(*args)
            if defect and defect.startswith("lineage_"): raise ValueError("synthetic_" + defect)
            return {"passed": True, "f4_shared_reproduction_passed": True,
                    "settings": [{"setting_id": k} for k in owner.F1_SHA256]}
        def git(cwd, *args):
            calls.append((cwd, args))
            if args[:2] == ("worktree", "add"):
                worktree = Path(args[-2]); worktree.mkdir()
                (worktree / "testing").mkdir()
                (worktree / owner.PROFILE).write_bytes((ROOT / owner.PROFILE).read_bytes())
                path = worktree / "src/models/xmodel_kaon_pl.f"
                path.parent.mkdir(parents=True); path.write_bytes(b"committed model")
            if args[0] == "branch": return "other" if defect == "branch" and cwd == repo else "test" if cwd == repo else ""
            if args[0] == "rev-parse":
                if cwd == repo and ((defect == "head" and args[1] == "HEAD") or
                                    (defect == "origin" and args[1] != "HEAD")): return "b" * 40
                return SHA
            if args[0] == "status":
                if cwd != repo: return ""
                if defect == "dirty": return dirty + "\n?? testing/unreviewed.py"
                return dirty + ("\n M src/main.py" if defect == "status_changed" and ran else "")
            return ""
        def probes(worktree, env=None, baseline=None):
            if defect == "path": raise ValueError("ltanapath_not_isolated")
            return [{"caller": str(worktree), "paths": {"VOLATILEPATH": str(root / "volatile")}, "isolation_passed": True}]
        def external(*args):
            if defect == "external": raise ValueError("background_simc_worksim_would_mutate")
            return {"passed": True}
        def run(worktree, log, env=None):
            nonlocal ran
            ran = True; events.append("isolated_analysis")
            self.assertNotEqual(worktree, repo)
            self.assertEqual(worktree.name, "analysis")
            (worktree / "src/models/xmodel_kaon_pl.f").write_bytes(b"runtime rewritten")
            self.assertEqual((repo / "src/models/xmodel_kaon_pl.f").read_bytes(), primary_bytes["src/models/xmodel_kaon_pl.f"])
            if defect == "hash_changed": (repo / "src/models/xmodel_kaon_pl.f").write_bytes(b"unexpected primary mutation")
            text = "Low Epsilon Completed!\nHigh Epsilon Completed!\n"
            if defect == "markers": text = "Left / lowe debug analysis completed.\n"
            log.write_text(text)
            for name in owner.artifact_names(profile) - {owner.SUMMARY} - set(candidates):
                if name.endswith("-manifest.json"):
                    setting = next(s for s in owner.CANONICAL_SETTINGS if owner.pdf_name(**s).replace('.pdf', '-manifest.json') == name)
                    content = json.dumps(pages(**setting)).encode()
                elif name.endswith(".pdf"): content = b"%PDF-fixture"
                else:
                    epsilon = "lowe" if "lowe" in name else "highe"
                    phis = [s["phi"] for s in owner.CANONICAL_SETTINGS if s["epsilon"] == epsilon]
                    if "ledger" not in name:
                        content = json.dumps({"histlist": [{"phi_setting": phi} for phi in phis],
                            "inpDict": {"ParticleType": "kaon", "EPSSET": "low" if epsilon == "lowe" else "high",
                                "Q2": "4p4", "W": "2p74", "OutFilename": f"FullAnalysis_{owner.KINEMATIC}_{epsilon}"}}).encode()
                    elif name.endswith(".json"):
                        content = json.dumps({"active_profile": "no_empirical_residual", "particle_type": "kaon",
                            "epsset": "low" if epsilon == "lowe" else "high", "q2": "4p4", "w": "2p74",
                            "settings": [{"phi_setting": phi} for phi in phis]}).encode()
                    else:
                        content = ("row_kind,phi_setting\n" + "".join("setting_total," + phi + "\n" for phi in phis)).encode()
                (outdir / name).write_bytes(content)
                # Synthetic freshness is test-owned; avoid filesystem clock skew.
                fresh_ns = time.time_ns() + 1_000_000_000
                os.utime(outdir / name, ns=(fresh_ns, fresh_ns))
            target = outdir / owner.pdf_name("Center", "highe")
            if defect == "stale": os.utime(target, (1, 1))
            if defect == "missing": target.unlink()
            if defect == "candidate": (outdir / next(iter(candidates))).write_bytes(b"wrong candidate")
            return 7 if defect == "analysis" else 0
        real_preserve, real_verify, real_zip = owner.verify_preservation, owner.verify_artifacts, owner.verify_zip
        def preserve(*args):
            events.append("preservation_check")
            if defect == "final_preservation": raise ValueError("synthetic_final_preservation_failure")
            return real_preserve(*args)
        def verify(*args):
            events.append("five_setting_verify")
            return real_verify(*args)
        checks = [{"name": "required_analysis_commit_ancestor", "returncode": 0, "stdout": "", "stderr": ""},
                  {"name": "committed_files_after_required_analysis_commit", "returncode": 0, "stdout": "", "stderr": ""}]
        def source_checks(*args):
            events.append("source_recheck" if ran else "source_preflight")
            if defect == "source_preflight": raise ValueError("synthetic_source_failure")
            return checks
        real_collect = owner.collector.collect_validation_bundle
        def collect(**kwargs):
            events.append("collection")
            self.assertNotEqual(kwargs["repo_root"], repo)
            self.assertNotIn("phi", kwargs); self.assertNotIn("epsilon", kwargs)
            if defect == "collector": return {"returncode": 1}
            def command_runner(command, cwd=None):
                return {"command": command, "returncode": 0, "stdout": SHA + "\n" if command[:3] == ["git", "rev-parse", "HEAD"] else "", "stderr": ""}
            return real_collect(**kwargs, command_runner=command_runner)
        @contextmanager
        def module(worktree):
            yield SimpleNamespace(collect_validation_bundle=collect)
        def verify_zip(*args):
            events.append("zip_verify")
            if defect == "zip": raise ValueError("synthetic_zip_failure")
            return real_zip(*args)
        real_copy = owner.copy_companion
        def copy(source, destination):
            name = next(k for k, p in companions.items() if p == destination)
            events.append("deliver_" + name)
            self.assertEqual(sum(args[:2] == ("worktree", "add") for _, args in calls),
                             sum(args[:2] == ("worktree", "remove") for _, args in calls))
            self.assertIn("zip_verify", events)
            self.assertEqual(events[events.index("zip_verify") + 1], "preservation_check")
            zip_hash = owner.sha256(output)
            if defect == "copy_" + name:
                destination.write_bytes(b"partial failed copy")
                raise OSError("synthetic_companion_copy_failure:" + name)
            real_copy(source, destination)
            self.assertEqual(owner.sha256(output), zip_hash)
        real_ltsep_preserve = owner.verify_ltsep_preservation
        def ltsep_preserve(*args):
            if preflight_only: events.append("ltsep_preservation_check")
            if defect == "final_ltsep": raise ValueError("synthetic_final_ltsep_failure")
            return real_ltsep_preserve(*args)
        with ExitStack() as stack:
            for name, replacement in {"git": git, "CANDIDATES": candidates, "probe_paths": probes,
                    "verify_candidate_materialization": materialization, "stage_candidates": staging,
                    "lineage_preflight": lineage,
                    "prepare_runtime_overlay": lambda *args: ({}, {"baseline_paths": {}, "original_source_sha256": {}}),
                    "external_symlink_preflight": external, "run_analysis": run,
                    "verify_preservation": preserve, "verify_artifacts": verify,
                    "collector_source_preflight": source_checks, "collection_module": module,
                    "verify_zip": verify_zip, "copy_companion": copy,
                    "verify_ltsep_preservation": ltsep_preserve}.items():
                stack.enter_context(patch.object(owner, name, replacement))
            if preflight_only:
                # Explicit tripwires: even an ignored return value is forbidden.
                for name in ("run_analysis", "verify_completion", "verify_artifacts",
                             "collector_source_preflight", "collection_module",
                             "verify_zip", "deliver_evidence", "copy_companion", "artifact_names"):
                    stack.enter_context(patch.object(owner, name,
                        Mock(side_effect=AssertionError("preflight-only called " + name))))
            try:
                if cli:
                    rc = owner.main(["--source-commit", SHA, "--repo", str(repo),
                                     "--outdir", str(outdir), "--output", str(output),
                                     "--candidate-materialization-dir", str(root / "materialization")] +
                                    (["--lineage-preflight-only"] if preflight_only else []))
                    result = (owner.accepted.gate_status_path(outdir, output) if preflight_only else output) if rc == 0 else None
                else:
                    result = owner.execute_gate(repo, outdir, output, SHA, root / "materialization",
                                                lineage_preflight_only=preflight_only)
            finally:
                if preflight_only:
                    self.assertFalse(ran)
                    self.assertFalse(output.exists())
                    adds = [args for _, args in calls if args[:2] == ("worktree", "add")]
                    removes = [args for _, args in calls if args[:2] == ("worktree", "remove")]
                    self.assertEqual(len(adds), len(removes))
                    self.assertTrue(all(not Path(args[-2]).exists() for args in adds))
                for cwd, args in calls:
                    self.assertNotIn(args[0], {"clean", "reset", "stash", "checkout"})
                    if args[:2] == ("worktree", "remove"):
                        self.assertEqual(args[2], "--force")
                        self.assertFalse(Path(args[-1]).is_relative_to(repo))
                if defect not in {"hash_changed"}:
                    self.assertEqual({name: (repo / name).read_bytes() for name in primary_bytes}, primary_bytes)
                if defect in {"branch", "head", "origin", "dirty"}:
                    self.assertFalse(any(args[0] == "worktree" for _, args in calls))
                    self.assertFalse(ran)
                if defect in {"analysis", "markers", "hash_changed", "status_changed", "stale", "missing", "candidate", "source_preflight", "path", "external"}:
                    self.assertNotIn("collection", events)
                if defect in {"path", "external"}: self.assertFalse(ran)
                if defect in {"materialization", "staging"} or (defect and defect.startswith("lineage_")):
                    self.assertFalse(ran)
                    self.assertNotIn("collection", events)
        return result, events, outdir


    def materialization_fixture(self, root):
        directory = root / "materialization"; directory.mkdir()
        f1 = [{"alias": k, "setting_id": k, "source_file_sha256": v,
               "setting": {"kinematic_token": owner.KINEMATIC}} for k,v in owner.F1_SHA256.items()]
        comparison = {"schema_version": "method_a_current_baseline_authority_comparison/v1",
            "non_authoritative": True, "accepted_file_sha256": owner.HISTORICAL_INPUT_SHA256,
            "candidate_serialized_sha256": {}, "diagnostic_f3_authority_override": deepcopy(owner.F3_RECONSTRUCTION),
            "f1_inputs": f1, "summary": dict(owner.SCIENTIFIC_GATE)}
        payloads, pins, outputs = {}, {}, {}
        for stage, body_name, fpkey in (("f2","representation","representation_fingerprint"),
                ("f3","acceptance_map","map_fingerprint"),("f4","correction","correction_fingerprint")):
            fp = owner.CANDIDATE_FINGERPRINTS[stage]
            body = {"fingerprint": fp[fpkey]}
            if stage == "f3": body["algorithm_fingerprint"] = fp["algorithm_fingerprint"]
            if stage == "f4":
                body.update({"f3_source_file_sha256": "pending", **{"f3_"+k:v for k,v in owner.CANDIDATE_FINGERPRINTS["f3"].items()}})
            payloads[stage] = {"non_authoritative": True, "artifact_fingerprint": fp["artifact_fingerprint"], body_name: body}
            if stage == "f4": body["f3_source_file_sha256"] = pins["f3"]
            path = directory / owner.MATERIALIZATION_NAMES[stage]
            path.write_text(json.dumps(payloads[stage]),encoding="utf-8"); pins[stage] = owner.sha256(path)
            outputs[stage] = {"basename": path.name, "raw_sha256": pins[stage], "non_authoritative": True, **fp}
            if stage == "f4": outputs[stage].update({k:v for k,v in body.items() if k.startswith("f3_")})
        comparison["candidate_serialized_sha256"] = {k:pins[k] for k in ("f2","f3")}
        comparison["diagnostic_f3_authority_override"][owner.KINEMATIC]["source_file_sha256"] = pins["f3"]
        path = directory / owner.MATERIALIZATION_NAMES["comparison"]
        path.write_text(json.dumps(comparison),encoding="utf-8"); pins["comparison"] = owner.sha256(path)
        manifest = {"schema_version": "method_a_current_baseline_authority_materialization/v1",
            "source_head": owner.MATERIALIZATION_HEAD, "kinematic_token": owner.KINEMATIC,
            "complete": True, "errors": [], "non_authoritative": True,
            "accepted_authority_mutated": False, "production_objects_mutated": False,
            "production_application_performed": False, "method_a_promoted": False,
            "accepted_inputs": {k:{"raw_sha256":v} for k,v in owner.HISTORICAL_INPUT_SHA256.items()},
            "comparison_input": {"raw_sha256": pins["comparison"],"copied_output_raw_sha256": pins["comparison"],
                "copied_output_basename": path.name,"schema_version": comparison["schema_version"]},
            "current_f1_inputs": deepcopy(f1), "scientific_gate": dict(owner.SCIENTIFIC_GATE),
            "candidate_outputs": outputs, "diagnostic_f3_authority_override": {"accepted_authority": False,
                "purpose": "diagnostic_candidate_construction_only", "record": comparison["diagnostic_f3_authority_override"]}}
        path = directory / owner.MATERIALIZATION_NAMES["manifest"]
        path.write_text(json.dumps(manifest),encoding="utf-8"); pins["manifest"] = owner.sha256(path)
        return directory, pins, manifest, comparison

    def materialization_patches(self, pins):
        stack = ExitStack()
        stack.enter_context(patch.object(owner,"MATERIALIZATION_SHA256",pins))
        stack.enter_context(patch.object(owner,"CANDIDATES",{owner.MATERIALIZATION_NAMES[k]:pins[k] for k in ("f2","f3","f4")}))
        authority = deepcopy(owner.F3_RECONSTRUCTION)
        authority[owner.KINEMATIC]["source_file_sha256"] = pins["f3"]
        stack.enter_context(patch.object(owner,"F3_RECONSTRUCTION",authority))
        return stack

    def test_fresh_candidate_mapping_is_independent_of_historical_owner(self):
        self.assertEqual(owner.CANDIDATES, {
            owner.MATERIALIZATION_NAMES["f2"]: owner.MATERIALIZATION_SHA256["f2"],
            owner.MATERIALIZATION_PREFIX+"acceptance-map-current-baseline-candidate.json": "c5b86452b790ecbaf5b8f0df05da67efa2fa92aab157b12b153ed7e491a38228",
            owner.MATERIALIZATION_PREFIX+"parent-preserving-correction-current-baseline-candidate.json": "79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7"})
        self.assertEqual(owner.accepted.CANDIDATES, {
            owner.MATERIALIZATION_NAMES["f3"]: "eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d",
            owner.MATERIALIZATION_NAMES["f4"]: "1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902"})
        self.assertIsNot(owner.CANDIDATES, owner.accepted.CANDIDATES)

    def test_reviewed_materialization_verifies_and_stages_exact_f2_f3_f4(self):
        for previous in ("absent","old","fresh","unknown"):
            with self.subTest(previous=previous), tempfile.TemporaryDirectory() as d:
                root=Path(d); directory,pins,manifest,comparison=self.materialization_fixture(root)
                outdir=root/"output";outdir.mkdir()
                original={p.name:p.read_bytes() for p in directory.iterdir()}
                old={}
                for k in ("f2","f3","f4"):
                    name=owner.MATERIALIZATION_NAMES[k]; target=outdir/name
                    if previous in {"old","unknown"}: target.write_bytes(b"known-old" if previous=="old" else b"unknown")
                    if previous=="fresh" or (previous=="old" and k=="f2"): target.write_bytes((directory/name).read_bytes())
                    if k != "f2": old[name]=__import__("hashlib").sha256(b"known-old").hexdigest()
                with self.materialization_patches(pins), patch.object(owner.accepted,"CANDIDATES",old):
                    verification=owner.verify_candidate_materialization(directory)
                    self.assertTrue(verification["passed"])
                    if previous=="unknown":
                        with self.assertRaisesRegex(ValueError,"unknown_candidate_target"):
                            owner.stage_candidates(directory,outdir,verification)
                        self.assertTrue(all(p.read_bytes()==b"unknown" for p in outdir.iterdir()))
                    else:
                        with patch.object(owner.os,"replace", wraps=os.replace) as replace:
                            installed=owner.stage_candidates(directory,outdir,verification)
                        self.assertEqual(replace.call_count,0 if previous=="fresh" else 2 if previous=="old" else 3)
                        self.assertEqual({p.name for p in outdir.iterdir()},set(owner.CANDIDATES))
                        for row in installed["files"]:
                            target=Path(row["path"])
                            self.assertEqual(target.read_bytes(),original[target.name])
                            self.assertEqual(row["action"], "no-op" if previous=="old" and target.name==owner.MATERIALIZATION_NAMES["f2"] else {"absent":"installed","old":"replaced","fresh":"no-op"}[previous])
                self.assertEqual({p.name:p.read_bytes() for p in directory.iterdir()},original)

    def test_materialization_rejects_bad_hashes_and_semantics(self):
        mutations=(
            ("source",lambda m:m.update(source_head="b"*40)),
            ("gate",lambda m:m["scientific_gate"].update(f3_scientific_payload_match=False)),
            ("f1",lambda m:m["current_f1_inputs"][0].update(source_file_sha256="b"*64)),
            ("f3",lambda m:m["candidate_outputs"]["f3"].update(map_fingerprint="b"*64)),
            ("f4",lambda m:m["candidate_outputs"]["f4"].update(artifact_fingerprint="b"*64)),
            ("history",lambda m:m["accepted_inputs"]["f2"].update(raw_sha256="b"*64)),
            ("sentinel",lambda m:m["diagnostic_f3_authority_override"]["record"][owner.KINEMATIC].update(farm_source_head="b"*40)))
        for key in owner.MATERIALIZATION_NAMES:
            with self.subTest(hash=key), tempfile.TemporaryDirectory() as d:
                directory,pins,manifest,comparison=self.materialization_fixture(Path(d))
                (directory/owner.MATERIALIZATION_NAMES[key]).write_bytes(b"wrong")
                with self.materialization_patches(pins), self.assertRaisesRegex(ValueError,"hash_mismatch"):
                    owner.verify_candidate_materialization(directory)
        for label,change in mutations:
            with self.subTest(semantic=label), tempfile.TemporaryDirectory() as d:
                directory,pins,manifest,comparison=self.materialization_fixture(Path(d))
                change(manifest);path=directory/owner.MATERIALIZATION_NAMES["manifest"]
                path.write_text(json.dumps(manifest),encoding="utf-8");pins["manifest"]=owner.sha256(path)
                with self.materialization_patches(pins), self.assertRaises(ValueError):
                    owner.verify_candidate_materialization(directory)

    def test_f2_predecessor_policy_is_one_exact_staging_only_identity(self):
        predecessor = "182433ccd2d13d3f7af80b963e69190ad8f2eb4e295c8216b36585e8ee7409b2"
        name = owner.MATERIALIZATION_NAMES["f2"]
        self.assertEqual(owner.RECOGNIZED_STAGING_PREDECESSORS, {name: predecessor})
        self.assertEqual(owner.CANDIDATES[name],
                         "2fa715b2b5c2d5077e38416fb9eceb814e013e72f308c43526fddff36ddd962e")
        self.assertNotIn(predecessor, owner.MATERIALIZATION_SHA256.values())
        self.assertNotIn(predecessor, owner.CANDIDATES.values())
        self.assertNotIn(predecessor, owner.HISTORICAL_INPUT_SHA256.values())
        self.assertNotIn(predecessor, owner.accepted.CANDIDATES.values())

    def test_observed_f2_predecessor_replaced_with_real_copy_and_final_hash(self):
        # Observed farm bytes are not local. Mock only the pre-existing target's
        # planning/recheck hashes; synthetic reviewed sources, temporary copies
        # and final destinations use the real hash reader and fixture pins.
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            directory, pins, _, _ = self.materialization_fixture(root)
            outdir = root / "output"; outdir.mkdir()
            name = owner.MATERIALIZATION_NAMES["f2"]
            target = outdir / name
            marker = b"synthetic predecessor; not the observed farm bytes"
            target.write_bytes(marker)
            predecessor = owner.RECOGNIZED_STAGING_PREDECESSORS[name]
            real_sha = owner.sha256
            predecessor_reads = []
            original = {p.name: p.read_bytes() for p in directory.iterdir()}
            def observed_hash(path):
                if Path(path) == target and target.read_bytes() == marker:
                    predecessor_reads.append(str(path))
                    return predecessor
                return real_sha(path)
            status = Mock()
            with self.materialization_patches(pins), patch.object(owner, "sha256", side_effect=observed_hash), \
                    patch.object(owner.os, "replace", wraps=os.replace) as replace:
                verification = owner.verify_candidate_materialization(directory)
                result = owner.stage_candidates(directory, outdir, verification, status=status)
            self.assertTrue(result["passed"])
            self.assertEqual(len(predecessor_reads), 2)
            self.assertEqual(replace.call_count, 3)
            self.assertEqual(result["files"][0], {"path": str(target),
                "before_sha256": predecessor, "after_sha256": pins["f2"], "action": "replaced"})
            self.assertEqual([r["action"] for r in result["files"]], ["replaced", "installed", "installed"])
            self.assertEqual(status.update.call_args.kwargs["candidate_installation"]["files"], result["files"])
            for stage in ("f2", "f3", "f4"):
                destination = outdir / owner.MATERIALIZATION_NAMES[stage]
                self.assertEqual(destination.read_bytes(), original[destination.name])
                self.assertEqual(real_sha(destination), pins[stage])
            self.assertEqual({p.name: p.read_bytes() for p in directory.iterdir()}, original)

    def test_other_f2_hashes_still_fail_before_any_replacement(self):
        predecessor = owner.RECOGNIZED_STAGING_PREDECESSORS[owner.MATERIALIZATION_NAMES["f2"]]
        for unknown in (predecessor[:-1] + "3", "0" * 64,
                        owner.HISTORICAL_INPUT_SHA256["f2"],
                        owner.accepted.CANDIDATES[owner.MATERIALIZATION_NAMES["f3"]]):
            with self.subTest(hash=unknown), tempfile.TemporaryDirectory() as d:
                root = Path(d)
                directory, pins, _, _ = self.materialization_fixture(root)
                outdir = root / "output"; outdir.mkdir()
                target = outdir / owner.MATERIALIZATION_NAMES["f2"]
                target.write_bytes(b"unknown synthetic target")
                real_sha = owner.sha256
                def unknown_hash(path):
                    return unknown if Path(path) == target else real_sha(path)
                with self.materialization_patches(pins), patch.object(owner, "sha256", side_effect=unknown_hash), \
                        patch.object(owner.os, "replace") as replace:
                    verification = owner.verify_candidate_materialization(directory)
                    with self.assertRaisesRegex(ValueError, "unknown_candidate_target:" + target.name):
                        owner.stage_candidates(directory, outdir, verification)
                replace.assert_not_called()
                self.assertEqual(list(outdir.iterdir()), [target])
                self.assertEqual(target.read_bytes(), b"unknown synthetic target")

    def test_staging_rejects_changed_source_without_replacement(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            directory, pins, _, _ = self.materialization_fixture(root)
            outdir = root / "output"; outdir.mkdir()
            with self.materialization_patches(pins), patch.object(owner.os, "replace") as replace:
                verification = owner.verify_candidate_materialization(directory)
                source = directory / owner.MATERIALIZATION_NAMES["f2"]
                source.write_bytes(source.read_bytes() + b"changed")
                with self.assertRaisesRegex(ValueError, "materialization_source_changed:" + source.name):
                    owner.stage_candidates(directory, outdir, verification)
            replace.assert_not_called()
            self.assertEqual(list(outdir.iterdir()), [])

    def test_predecessor_target_change_during_copy_prevents_replacement(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            directory, pins, _, _ = self.materialization_fixture(root)
            outdir = root / "output"; outdir.mkdir()
            target = outdir / owner.MATERIALIZATION_NAMES["f2"]
            target.write_bytes(b"synthetic predecessor")
            real_sha = owner.sha256
            reads = []
            def observed_then_changed(path):
                if Path(path) == target:
                    reads.append(str(path))
                    if len(reads) == 1:
                        return owner.RECOGNIZED_STAGING_PREDECESSORS[target.name]
                    target.write_bytes(b"concurrent change")
                return real_sha(path)
            with self.materialization_patches(pins), patch.object(owner, "sha256", side_effect=observed_then_changed), \
                    patch.object(owner.os, "replace") as replace:
                verification = owner.verify_candidate_materialization(directory)
                with self.assertRaisesRegex(ValueError, "candidate_target_changed_during_staging"):
                    owner.stage_candidates(directory, outdir, verification)
            replace.assert_not_called()
            self.assertEqual(len(reads), 2)
            self.assertEqual(list(outdir.iterdir()), [target])
            self.assertEqual(target.read_bytes(), b"concurrent change")

    def test_staging_rejects_symlink_nonfile_and_source_target_alias(self):
        for defect in ("symlink", "nonfile", "alias"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as d:
                root = Path(d)
                directory, pins, _, _ = self.materialization_fixture(root)
                original = {p.name: p.read_bytes() for p in directory.iterdir()}
                outdir = root / "output"; outdir.mkdir()
                target = outdir / owner.MATERIALIZATION_NAMES["f2"]
                if defect == "symlink": target.symlink_to(directory / target.name)
                if defect == "nonfile": target.mkdir()
                if defect == "alias": outdir = directory
                reason = {"symlink": "candidate_source_target_alias", "nonfile": "candidate_target_not_file",
                          "alias": "candidate_destination_inside_materialization"}[defect]
                with self.materialization_patches(pins), patch.object(owner.os, "replace") as replace:
                    verification = owner.verify_candidate_materialization(directory)
                    with self.assertRaisesRegex(ValueError, reason):
                        owner.stage_candidates(directory, outdir, verification)
                replace.assert_not_called()
                self.assertEqual({p.name: p.read_bytes() for p in directory.iterdir()}, original)

    def test_real_five_setting_lineage_preflight_loads_files_and_reconstructs_f4(self):
        from testing import test_f6_3_parallel_full_procedure_method_a as scientific
        fixture=type("LocalLineageFixture",(scientific.CandidateLineageTests,),{})
        fixture.setUpClass()
        module=scientific.f63
        with tempfile.TemporaryDirectory() as d:
            outdir=Path(d); paths=module.accepted_f6_3_artifact_paths(outdir,owner.KINEMATIC)
            hashes={}
            for artifact in fixture.artifacts:
                setting=artifact["setting"]
                sid=setting["phi_setting"]+"-"+setting["epsilon_filename_token"]
                path=Path(paths["f1"][sid]);path.write_text(json.dumps(artifact,sort_keys=True),encoding="utf-8")
                hashes[sid]=owner.sha256(path)
            f2=scientific.f2.build_pion_hgcer_method_a_acceptance_representation_artifact(fixture.artifacts,input_file_hashes=hashes,input_paths={})
            Path(paths["f2"]).write_bytes(module._writer_bytes(f2))
            f2sha=owner.sha256(Path(paths["f2"]))
            f2authority={"source_file_sha256":f2sha,"representation_fingerprint":f2["representation"]["fingerprint"],"artifact_fingerprint":f2["artifact_fingerprint"]}
            f3=scientific.f3_builder.build_pion_hgcer_method_a_acceptance_map_artifact(fixture.artifacts,f2,f1_input_file_hashes=hashes,f2_input_file_sha256=f2sha,input_paths={})
            Path(paths["f3"]).write_text(json.dumps(f3,sort_keys=True),encoding="utf-8")
            f3sha=owner.sha256(Path(paths["f3"]))
            f3authority=scientific.f5_fixtures.f4_fixtures._authority(f3,f3sha)
            f3authority[owner.KINEMATIC]["farm_source_head"]="0"*40
            f4=scientific.f4.build_pion_hgcer_method_a_parent_preserving_correction_artifact(
                fixture.artifacts,f3,f1_input_file_hashes=hashes,f3_input_file_sha256=f3sha,
                accepted_f3_runtime_authority_by_kinematic=f3authority,input_paths={})
            Path(paths["f4"]).write_text(json.dumps(f4,sort_keys=True),encoding="utf-8")
            f4sha=owner.sha256(Path(paths["f4"]))
            f4authority={owner.KINEMATIC:{"source_file_sha256":f4sha,
                "correction_fingerprint":f4["correction"]["fingerprint"],
                "artifact_fingerprint":f4["artifact_fingerprint"],"farm_source_head":owner.MATERIALIZATION_HEAD,
                "f1_source_file_sha256":hashes,**{k:f4["correction"][k] for k in
                    ("f3_source_file_sha256","f3_map_fingerprint","f3_algorithm_fingerprint","f3_artifact_fingerprint")}}}
            before={p.name:p.read_bytes() for p in outdir.iterdir()}
            def real_lineage(*args):
                # The outer orchestration fixture owns different synthetic
                # candidate bytes; restore this calculator fixture's identities.
                with patch.object(owner,"CANDIDATES",{
                        Path(paths["f2"]).name:f2sha,Path(paths["f3"]).name:f3sha,Path(paths["f4"]).name:f4sha}):
                    return owner.validate_candidate_lineage(outdir,module)
            with patch.object(owner,"F1_SHA256",hashes), patch.object(owner,"CANDIDATES",{
                    Path(paths["f2"]).name:f2sha,Path(paths["f3"]).name:f3sha,Path(paths["f4"]).name:f4sha}), \
                 patch.object(module,"F6_3_CANDIDATE_F1_SOURCE_FILE_SHA256",hashes), \
                 patch.object(module,"F6_3_CANDIDATE_F2_AUTHORITY",f2authority), \
                 patch.object(module,"F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC",f3authority), \
                 patch.object(module,"F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC",f4authority), \
                 patch.object(module,"reconstruct_transient_factor_map",wraps=module.reconstruct_transient_factor_map) as reconstruct:
                result=owner.validate_candidate_lineage(outdir,module)
                self.assertTrue(result["passed"]);self.assertEqual(reconstruct.call_count,5)
                self.assertEqual([r["setting_id"] for r in result["settings"]],list(hashes))
                self.assertTrue(all(r["factor_count"]>0 and r["f4_shared_reproduction_passed"] for r in result["settings"]))
                self.assertEqual(result["observed_sha256"],{**hashes,"f2":f2sha,"f3":f3sha,"f4":f4sha})
                with tempfile.TemporaryDirectory() as flow:
                    receipt, events, _ = self.exercise(flow, lineage_override=real_lineage,
                                                       preflight_only=True)
                    transferred = json.loads(receipt.read_text())
                    self.assertEqual(reconstruct.call_count, 10)
                    self.assertEqual(transferred["f6_3_lineage_preflight"], result)
                    self.assertFalse(transferred["analysis_started"])
                    self.assertEqual(transferred["status"], "success")
                    self.assertEqual(events, ["materialization_verify", "candidate_staging",
                        "lineage_preflight", "preservation_check", "ltsep_preservation_check"])
                self.assertTrue(all(r["scientific_equivalence_mode"] == "exact-lineage"
                                    for r in result["settings"]))
                original_f1 = {sid: Path(path).read_bytes() for sid, path in paths["f1"].items()}
                first_sid = next(iter(paths["f1"]))
                Path(paths["f1"][first_sid]).write_bytes(original_f1[first_sid] + b" ")
                mixed = owner.validate_candidate_lineage(outdir, module)
                self.assertTrue(all(r["scientific_equivalence_mode"] == "equivalent-new-lineage"
                                    for r in mixed["settings"]))
                self.assertEqual({sid for sid in hashes if
                    mixed["observed_f1_source_file_sha256"][sid] != hashes[sid]}, {first_sid})
                # Serialization alone changes every raw hash without changing
                # any scientific object; the real bridge rebuilds all stages.
                for sid, path in paths["f1"].items():
                    Path(path).write_bytes(original_f1[sid] + b"\n ")
                observed = {sid: owner.sha256(Path(path)) for sid, path in paths["f1"].items()}
                regenerated = owner.validate_candidate_lineage(outdir, module)
                self.assertTrue(all(observed[sid] != hashes[sid] for sid in hashes))
                self.assertEqual(regenerated["reviewed_f1_source_file_sha256"], hashes)
                self.assertEqual(regenerated["observed_f1_source_file_sha256"], observed)
                for record in regenerated["settings"]:
                    self.assertEqual(record["scientific_equivalence_mode"], "equivalent-new-lineage")
                    provenance = record["provenance"]
                    self.assertEqual(provenance["accepted_f1_source_file_sha256"], observed)
                    self.assertEqual(provenance["current_runtime_lineage"]["f1_source_file_sha256"], observed)
                    self.assertEqual(provenance["reviewed_candidate"]["f1_source_file_sha256"], hashes)
                    self.assertTrue(all(provenance["scientific_equivalence"][stage] ==
                        {"passed": True, "first_mismatch_path": None} for stage in ("f2", "f3", "f4")))
                for preflight_only in (False, True):
                    with tempfile.TemporaryDirectory() as flow:
                        receipt, _, flow_outdir = self.exercise(flow, lineage_override=real_lineage,
                                                               preflight_only=preflight_only)
                        status_path = receipt if preflight_only else flow_outdir / "fresh-gate-status.json"
                        self.assertEqual(json.loads(status_path.read_text())["f6_3_lineage_preflight"], regenerated)
                real = module.reconstruct_transient_factor_map
                def reject_in_both_modes(reason):
                    with self.assertRaisesRegex(Exception, reason):
                        owner.validate_candidate_lineage(outdir, module)
                    for preflight_only in (False, True):
                        with tempfile.TemporaryDirectory() as flow:
                            with self.assertRaisesRegex(Exception, reason):
                                self.exercise(flow, lineage_override=real_lineage, preflight_only=preflight_only)
                            status = json.loads((Path(flow) / "volatile/OUTPUT/Analysis/KaonLT/fresh-gate-status.json").read_text())
                            self.assertFalse(status["analysis_started"])
                # Tamper only with returned evidence after real reconstruction;
                # never patch scientific acceptance or builder results.
                def altered(decision):
                    def reconstruct_result(*args, **kwargs):
                        factors, provenance, rows = real(*args, **kwargs)
                        decision(factors, provenance, rows)
                        return factors, provenance, rows
                    return reconstruct_result
                manipulations = (
                    ("factor_parent_inventory_invalid", lambda f,p,r: r[next(iter(r))].update(t_index=99)),
                    ("setting_provenance_invalid", lambda f,p,r: p.update(transient_factor_identity_fingerprint="a"*64)),
                    ("setting_provenance_invalid", lambda f,p,r: p.update(accepted_f1_source_file_sha256=hashes)),
                    ("roles_provenance_invalid", lambda f,p,r: p["reviewed_candidate"].update(f1_source_file_sha256=observed)),
                    ("roles_provenance_invalid", lambda f,p,r: p["current_runtime_lineage"].update(f1_source_file_sha256=hashes)),
                    ("roles_provenance_invalid", lambda f,p,r: p["current_runtime_lineage"].pop("f3")),
                    ("scientific_equivalence_invalid", lambda f,p,r: p["scientific_equivalence"].update(mode="exact-lineage")),
                    ("scientific_equivalence_invalid", lambda f,p,r: p["scientific_equivalence"].update(schema_version="stale")),
                    ("scientific_equivalence_invalid", lambda f,p,r: p["scientific_equivalence"].update(provenance_exclusions={})),
                    ("scientific_equivalence_invalid", lambda f,p,r: p["scientific_equivalence"]["f2"].update(first_mismatch_path="science")),
                    ("scientific_equivalence_invalid", lambda f,p,r: p["scientific_equivalence"]["f4"].update(passed=False)),
                    ("scientific_equivalence_invalid", lambda f,p,r: p["scientific_equivalence"]["f3"].update(passed=1)),
                )
                for reason, manipulation in manipulations:
                    with self.subTest(reason=reason), patch.object(module, "reconstruct_transient_factor_map",
                            side_effect=altered(manipulation)):
                        reject_in_both_modes(reason)
                for sid, path in paths["f1"].items():
                    Path(path).write_bytes(original_f1[sid])
                first_sid = next(iter(paths["f1"]))
                first_path = Path(paths["f1"][first_sid])
                # A fully resealed scientific change is distinct from a raw
                # serialization change and must fail the first exact stage.
                changed = json.loads(original_f1[first_sid])
                changed["contract"]["method_a_training_records"][0]["SHMS_delta"] += 0.5
                scientific.f5_fixtures.f1_fixtures._seal_f1_artifact(changed)
                scientific.f3_builder._validate_f1_artifacts([changed] + fixture.artifacts[1:])
                first_path.write_bytes(module._writer_bytes(changed))
                reject_in_both_modes("f6_3_scientific_equivalence_mismatch:f2:")
                for malformed in (b'{"x":1,"x":2}', b'{"x":NaN}', b'{}'):
                    first_path.write_bytes(malformed)
                    reject_in_both_modes("f6_3_")
                first_path.unlink()
                reject_in_both_modes("accepted_artifact_load_failed")
                first_path.write_bytes(original_f1[first_sid])
                duplicate_sid = list(paths["f1"])[1]
                first_path.write_bytes(original_f1[duplicate_sid])
                reject_in_both_modes("f1_setting_duplicate")
                first_path.write_bytes(original_f1[first_sid])
                for inventory in (dict(list(paths["f1"].items())[1:]),
                                  {**paths["f1"], "Right-lowe": str(first_path)}):
                    with patch.object(module, "accepted_f6_3_artifact_paths",
                                      return_value={**paths, "f1": inventory}):
                        reject_in_both_modes("f1_path_inventory_invalid")
                loader = module.load_accepted_f6_3_authority
                for bad in ("", "g"*64, None):
                    def bad_hashes(*args, **kwargs):
                        loaded = loader(*args, **kwargs)
                        return (*loaded[:-1], {**loaded[-1], first_sid: bad})
                    with patch.object(module, "load_accepted_f6_3_authority", side_effect=bad_hashes):
                        reject_in_both_modes("lineage_f1_hash_inventory_invalid")
                for key in ("f2","f3","f4"):
                    path=Path(paths[key] if key in {"f2","f3","f4"} else paths["f1"][key])
                    original=path.read_bytes();path.write_bytes(original+b" ")
                    with self.subTest(hash=key),self.assertRaisesRegex(ValueError,"lineage_.*identity_mismatch"):
                        owner.validate_candidate_lineage(outdir,module)
                    for preflight_only in (False, True):
                        with tempfile.TemporaryDirectory() as flow, self.assertRaisesRegex(ValueError,"lineage_.*identity_mismatch"):
                            self.exercise(flow,"lineage_"+key, lineage_override=real_lineage, preflight_only=preflight_only)
                    path.write_bytes(original)
                invalid=deepcopy(f4authority);invalid[owner.KINEMATIC]["correction_fingerprint"]="b"*64
                with patch.object(module,"F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC",invalid),self.assertRaises(Exception):
                    owner.validate_candidate_lineage(outdir,module)
                for preflight_only in (False, True):
                    with patch.object(module,"F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC",invalid), tempfile.TemporaryDirectory() as flow, self.assertRaises(Exception):
                        self.exercise(flow,"lineage_authority",lineage_override=real_lineage, preflight_only=preflight_only)
                for preflight_only in (False, True):
                    with patch.object(scientific.f4,"build_pion_hgcer_method_a_parent_preserving_correction_with_review_data",side_effect=ValueError("synthetic_reconstruction_failure")), tempfile.TemporaryDirectory() as flow, self.assertRaisesRegex(Exception,"synthetic_reconstruction_failure"):
                        self.exercise(flow,"lineage_reconstruction",lineage_override=real_lineage, preflight_only=preflight_only)
                real=module.reconstruct_transient_factor_map
                def wrong_setting(*args,**kwargs):
                    factors,provenance,rows=real(*args,**kwargs)
                    provenance["selected_setting_id"]="Right-lowe"
                    return factors,provenance,rows
                with patch.object(module,"reconstruct_transient_factor_map",side_effect=wrong_setting), \
                        self.assertRaisesRegex(ValueError,"setting_provenance_invalid"):
                    owner.validate_candidate_lineage(outdir,module)
                for preflight_only in (False, True):
                    with patch.object(module,"reconstruct_transient_factor_map",side_effect=wrong_setting), tempfile.TemporaryDirectory() as flow, self.assertRaisesRegex(ValueError,"setting_provenance_invalid"):
                        self.exercise(flow,"lineage_setting",lineage_override=real_lineage, preflight_only=preflight_only)
                for value in (None, 0.0, -1.0, float("nan")):
                    def invalid_factors(*args, **kwargs):
                        factors, provenance, rows = real(*args, **kwargs)
                        if value is None: factors = {}
                        else: factors[next(iter(factors))] = value
                        return factors, provenance, rows
                    for preflight_only in (False, True):
                        with self.subTest(factor=value), patch.object(module,"reconstruct_transient_factor_map",side_effect=invalid_factors), tempfile.TemporaryDirectory() as flow, self.assertRaisesRegex(ValueError,"factor_population_invalid"):
                            self.exercise(flow,"lineage_setting",lineage_override=real_lineage, preflight_only=preflight_only)
                for key in ("accepted_f4_correction_fingerprint", "accepted_f4_artifact_fingerprint"):
                    def incomplete_provenance(*args, **kwargs):
                        factors, provenance, rows = real(*args, **kwargs)
                        provenance.pop(key)
                        return factors, provenance, rows
                    for preflight_only in (False, True):
                        with self.subTest(missing=key), patch.object(module,"reconstruct_transient_factor_map",side_effect=incomplete_provenance), tempfile.TemporaryDirectory() as flow, self.assertRaisesRegex(ValueError,"setting_provenance_invalid"):
                            self.exercise(flow,"lineage_setting",lineage_override=real_lineage, preflight_only=preflight_only)
            self.assertEqual({p.name:p.read_bytes() for p in outdir.iterdir()},before)


    def test_materialization_staging_and_all_lineage_failures_prevent_analysis(self):
        for defect in ("materialization","staging","lineage_f1","lineage_f3","lineage_f4",
                       "lineage_reconstruction","lineage_authority","lineage_setting"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as d:
                with self.assertRaisesRegex(ValueError,"synthetic_"):
                    self.exercise(d,defect)
                status=json.loads((Path(d)/"volatile/OUTPUT/Analysis/KaonLT/fresh-gate-status.json").read_text())
                self.assertEqual(status["status"],"failed")
                self.assertFalse(status["analysis_started"])

    def test_preflight_only_cli_returns_final_receipt_without_full_run_operations(self):
        with tempfile.TemporaryDirectory() as directory, redirect_stdout(io.StringIO()) as stdout:
            receipt, events, outdir = self.exercise(directory, cli=True, preflight_only=True)
            self.assertEqual(stdout.getvalue(), receipt.as_posix() + "\n")
            status = json.loads(receipt.read_text())
            self.assertEqual(status["source_commit"], SHA)
            self.assertEqual(status["mode"], "lineage-preflight-only")
            self.assertIn("not canonical-five runtime validation", status["role"])
            self.assertFalse(status["canonical_five_runtime_validation"])
            self.assertFalse(status["analysis_started"])
            self.assertTrue(status["candidate_materialization"]["passed"])
            self.assertTrue(status["candidate_installation"]["passed"])
            self.assertTrue(status["f6_3_lineage_preflight"]["f4_shared_reproduction_passed"])
            self.assertTrue(status["worktree_cleanup_completed"])
            self.assertTrue(status["final_ordinary_checkout_preservation"]["passed"])
            self.assertTrue(status["final_installed_ltsep_preservation"]["passed"])
            self.assertEqual(events, ["materialization_verify", "candidate_staging",
                "lineage_preflight", "preservation_check", "ltsep_preservation_check"])
            self.assertEqual({p.name for p in outdir.iterdir()},
                             set(owner.CANDIDATES) | {receipt.name})
            for flag in ("analysis_completed", "artifact_verification_completed",
                         "collection_completed", "zip_verification_completed"):
                self.assertFalse(status[flag])

    def test_preflight_only_fails_closed_including_final_preservation(self):
        for defect in ("materialization", "staging", "path", "external", "lineage_f1",
                       "lineage_f3", "lineage_f4", "lineage_reconstruction",
                       "lineage_authority", "lineage_setting",
                       "final_preservation", "final_ltsep"):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError):
                    self.exercise(directory, defect, preflight_only=True)
                status = json.loads((Path(directory) /
                    "volatile/OUTPUT/Analysis/KaonLT/fresh-gate-status.json").read_text())
                self.assertEqual(status["status"], "failed")
                self.assertFalse(status["analysis_started"])
                self.assertEqual(status["mode"], "lineage-preflight-only")

    def test_success_real_generic_collection_order_and_zip(self):
        with tempfile.TemporaryDirectory() as directory:
            output, events, outdir = self.exercise(directory)
            self.assertEqual(events, ["materialization_verify", "source_preflight", "candidate_staging",
                "lineage_preflight", "isolated_analysis", "preservation_check",
                "five_setting_verify", "collection", "source_recheck", "preservation_check",
                "zip_verify", "preservation_check", "deliver_log", "deliver_run_summary", "deliver_gate_status"])
            with zipfile.ZipFile(output) as archive:
                manifest = json.loads(archive.read("manifest.json"))
                self.assertTrue(manifest["complete"])
                self.assertEqual(len(manifest["settings"]), 5)
                self.assertEqual(manifest["requested_settings"], list(owner.CANONICAL_SETTINGS))
                self.assertEqual(set(manifest["global_artifacts"]),
                                 {"candidate_f2", "candidate_f3", "candidate_f4", "run_summary"})
            status = json.loads((outdir / "fresh-gate-status.json").read_text())
            self.assertEqual(status["status"], "success")
            self.assertEqual(status["zip_identity"]["sha256"], owner.sha256(output))
            self.assertTrue(status["worktree_cleanup_completed"])
            self.assertTrue(status["final_ordinary_checkout_preservation"]["passed"])
            self.assertTrue(status["final_installed_ltsep_preservation"]["passed"])
            summary = json.loads((outdir / owner.SUMMARY).read_text())
            self.assertEqual(len(summary["ordinary_checkout_before"]["farm_local_sha256"]), 2)
            self.assertTrue(summary["ordinary_checkout_before"]["porcelain"])
            self.assertTrue(summary["ordinary_checkout_preservation"]["passed"])
            self.assertTrue(summary["candidate_materialization"]["passed"])
            self.assertTrue(summary["candidate_installation"]["passed"])
            self.assertEqual(set(summary["candidate_sha256"]), set(owner.CANDIDATES))
            self.assertIn(owner.MATERIALIZATION_NAMES["f2"], summary["candidate_sha256"])
            self.assertTrue(summary["f6_3_lineage_preflight"]["passed"])
            self.assertEqual((output.parent / "fresh-gate-status.json").read_bytes(),
                             (outdir / "fresh-gate-status.json").read_bytes())
            for name, source in (("log", outdir / "fresh.log"), ("run_summary", outdir / owner.SUMMARY)):
                identity = status["companion_evidence"][name]
                delivered = Path(identity["path"])
                self.assertEqual(delivered.parent, output.parent)
                self.assertEqual(delivered.read_bytes(), source.read_bytes())
                self.assertEqual(identity["sha256"], owner.sha256(delivered))
                self.assertEqual(identity["bytes"], delivered.stat().st_size)
            self.assertEqual(status["delivered_gate_status_path"], str(output.parent / "fresh-gate-status.json"))

    def test_one_command_delivers_complete_evidence_without_followup(self):
        with tempfile.TemporaryDirectory() as directory, redirect_stdout(io.StringIO()) as stdout:
            output, events, outdir = self.exercise(directory, cli=True)
            self.assertEqual(stdout.getvalue(), output.as_posix() + "\n")
            self.assertEqual({p.name for p in output.parent.iterdir() if p.is_file()},
                             {"fresh.zip", "fresh.log", "fresh-run-summary.json", "fresh-gate-status.json"})
            status = json.loads((output.parent / "fresh-gate-status.json").read_text())
            self.assertEqual(status["status"], "success")
            self.assertEqual(status["zip_identity"]["sha256"], owner.sha256(output))
            self.assertEqual(events[-1], "deliver_gate_status")

    def test_existing_companions_fail_closed(self):
        for name in ("log", "run_summary", "gate_status"):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(ValueError, "companion_destination_exists"):
                    self.exercise(directory, "existing_" + name)
                root = Path(directory)
                status = json.loads((root / "volatile/OUTPUT/Analysis/KaonLT/fresh-gate-status.json").read_text())
                self.assertEqual(status["status"], "failed")
                existing = {"log": "fresh.log", "run_summary": "fresh-run-summary.json",
                            "gate_status": "fresh-gate-status.json"}[name]
                self.assertEqual((root / existing).read_bytes(), b"earlier attempt")

    def test_companion_copy_failures_return_failure_with_failed_source_status(self):
        for name in ("log", "run_summary", "gate_status"):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory, \
                    redirect_stdout(io.StringIO()) as stdout, redirect_stderr(io.StringIO()):
                result, events, outdir = self.exercise(directory, "copy_" + name, cli=True)
                self.assertIsNone(result)
                self.assertEqual(stdout.getvalue(), "")
                status = json.loads((outdir / "fresh-gate-status.json").read_text())
                self.assertEqual(status["status"], "failed")
                self.assertIn("synthetic_companion_copy_failure", status["failure_reason"])

    def test_companion_copy_exclusively_refuses_existing_destination(self):
        with tempfile.TemporaryDirectory() as directory:
            source, destination = Path(directory) / "source", Path(directory) / "destination"
            source.write_bytes(b"new attempt"); destination.write_bytes(b"old attempt")
            with self.assertRaises(FileExistsError): owner.copy_companion(source, destination)
            self.assertEqual(destination.read_bytes(), b"old attempt")

    def test_failure_gates_stop_and_persist(self):
        for defect, reason in (("branch", "wrong_branch"), ("head", "wrong_head"), ("origin", "not_observed_pushed"),
            ("dirty", "dirty_gate_source"), ("analysis", "analysis_failed"), ("markers", "completion_missing"),
            ("hash_changed", "ordinary_checkout_changed:farm_local_sha256"),
            ("status_changed", "ordinary_checkout_changed:porcelain"), ("stale", "artifact_stale"),
            ("missing", "artifact_missing"), ("candidate", "candidate_identity_mismatch"),
            ("source_preflight", "synthetic_source_failure"), ("collector", "collection_failed"), ("zip", "synthetic_zip_failure")):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(ValueError, reason): self.exercise(directory, defect)
                status = json.loads((Path(directory) / "volatile/OUTPUT/Analysis/KaonLT/fresh-gate-status.json").read_text())
                self.assertEqual(status["status"], "failed")
                self.assertIn(reason.split(":")[0], status["failure_reason"])

    def test_isolation_and_external_failures_block_before_analysis(self):
        for defect, reason in (("path", "ltanapath_not_isolated"),
                              ("external", "background_simc_worksim_would_mutate")):
            with self.subTest(defect=defect), tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(ValueError, reason): self.exercise(directory, defect)

    def test_completion_requires_both_markers_and_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "analysis.log"
            for text, code in (("Low Epsilon Completed!", 0), ("High Epsilon Completed!", 0),
                               ("Low Epsilon Completed!\nHigh Epsilon Completed!", 1)):
                log.write_text(text)
                with self.assertRaises(ValueError): owner.verify_completion(log, code)

    def test_main_stdout_is_only_zip_on_success(self):
        for failure in (False, True):
            stdout, stderr = io.StringIO(), io.StringIO()
            with patch.object(owner, "execute_gate", side_effect=ValueError("failed") if failure else None,
                              return_value=Path("/fresh.zip")), redirect_stdout(stdout), redirect_stderr(stderr):
                rc = owner.main(["--source-commit", SHA, "--outdir", "artifacts", "--output", "fresh.zip",
                                 "--candidate-materialization-dir", "materialization"])
            self.assertEqual(rc, 1 if failure else 0)
            self.assertEqual(stdout.getvalue(), "" if failure else "/fresh.zip\n")


if __name__ == "__main__":
    unittest.main()
