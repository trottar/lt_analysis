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

    def exercise(self, directory, defect=None, cli=False):
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
            target = outdir / owner.pdf_name("Center", "highe")
            if defect == "stale": os.utime(target, (1, 1))
            if defect == "missing": target.unlink()
            if defect == "candidate": (outdir / next(iter(candidates))).write_bytes(b"wrong candidate")
            return 7 if defect == "analysis" else 0
        real_preserve, real_verify, real_zip = owner.verify_preservation, owner.verify_artifacts, owner.verify_zip
        def preserve(*args):
            events.append("preservation_check")
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
        with ExitStack() as stack:
            for name, replacement in {"git": git, "CANDIDATES": candidates, "probe_paths": probes,
                    "prepare_runtime_overlay": lambda *args: ({}, {"baseline_paths": {}, "original_source_sha256": {}}),
                    "external_symlink_preflight": external, "run_analysis": run,
                    "verify_preservation": preserve, "verify_artifacts": verify,
                    "collector_source_preflight": source_checks, "collection_module": module,
                    "verify_zip": verify_zip, "copy_companion": copy}.items():
                stack.enter_context(patch.object(owner, name, replacement))
            try:
                if cli:
                    rc = owner.main(["--source-commit", SHA, "--repo", str(repo),
                                     "--outdir", str(outdir), "--output", str(output)])
                    result = output if rc == 0 else None
                else:
                    result = owner.execute_gate(repo, outdir, output, SHA)
            finally:
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
        return result, events, outdir

    def test_success_real_generic_collection_order_and_zip(self):
        with tempfile.TemporaryDirectory() as directory:
            output, events, outdir = self.exercise(directory)
            self.assertEqual(events, ["source_preflight", "isolated_analysis", "preservation_check",
                "five_setting_verify", "collection", "source_recheck", "preservation_check",
                "zip_verify", "preservation_check", "deliver_log", "deliver_run_summary", "deliver_gate_status"])
            with zipfile.ZipFile(output) as archive:
                manifest = json.loads(archive.read("manifest.json"))
                self.assertTrue(manifest["complete"])
                self.assertEqual(len(manifest["settings"]), 5)
                self.assertEqual(manifest["requested_settings"], list(owner.CANONICAL_SETTINGS))
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
                rc = owner.main(["--source-commit", SHA, "--outdir", "artifacts", "--output", "fresh.zip"])
            self.assertEqual(rc, 1 if failure else 0)
            self.assertEqual(stdout.getvalue(), "" if failure else "/fresh.zip\n")


if __name__ == "__main__":
    unittest.main()
