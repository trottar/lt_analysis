"""Synthetic packages/repositories only; never attest the real JLab package."""
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from copy import deepcopy
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from testing import attest_pion_chain_ltsep_default_root as tool

ROOT = Path(__file__).resolve().parents[1]


def run(repo, *args):
    return subprocess.run(args, cwd=repo, check=True, capture_output=True).stdout


@contextmanager
def fixture(behavior="", pathing_extra="", duplicate=False):
    with tempfile.TemporaryDirectory(prefix="kaonlt-ltsep-synthetic-") as name, tempfile.TemporaryDirectory(prefix="kaonlt-ltsep-installed-") as installed_name:
        root = Path(name).resolve()
        repo = root / "ordinary"
        repo.mkdir()
        for source in (*tool.FROZEN, tool.TOOL, tool.TEST):
            target = repo / source
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes((ROOT / source).read_bytes())
        run(repo, "git", "init", "-q", "-b", "test")
        run(repo, "git", "add", "--all")
        run(repo, "git", "-c", "user.name=fixture", "-c", "user.email=fixture@example.invalid",
            "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null", "commit", "-qm", "synthetic source")
        sha = tool.git(repo, "rev-parse", "HEAD").decode().strip()
        run(repo, "git", "update-ref", "refs/remotes/origin/test", sha)
        # Keep the fake installed package outside ANALYSISPATH, so the package
        # exclusion is independently tested rather than hidden by a parent ban.
        installed = Path(installed_name).resolve()
        package = installed / "ltsep"
        (package / "PATH_TO_DIR").mkdir(parents=True)
        (package / "__init__.py").write_text(
            'from pathlib import Path\nfrom .pathing import SetPath\n'
            'class Root:\n'
            ' def __init__(self,caller,*args):\n'
            '  p=SetPath(caller)\n'
            '  for k in p.values: setattr(self,k,p.getPath(k))\n'
            '  self.OUTPATH=str(Path(self.LTANAPATH)/"OUTPUT/Analysis/KaonLT")\n'
            + behavior, encoding="utf-8")
        (package / "pathing.py").write_text(
            'from pathlib import Path\n'
            'class SetPath:\n'
            ' def __init__(self,caller):\n'
            '  self.values={}\n'
            '  for f in (Path(__file__).parent/"PATH_TO_DIR").glob("*.path"):\n'
            '   d=dict(line.split("=",1) for line in f.read_text().splitlines())\n'
            '   if Path(caller).resolve().is_relative_to(Path(d["ANALYSISPATH"])): self.values=d\n'
            '  if not self.values: raise ValueError("no caller match")\n'
            ' def getPath(self,key): return self.values[key].replace("${USER}",self.values["USER"])\n'
            + pathing_extra, encoding="utf-8")
        fields = {k: str(root / "inputs" / k) for k in tool.FIELDS}
        fields.update(LTANAPATH=str(repo), ANALYSISPATH=str(root), USER="fixture-user",
                      HOST="fixture-host", ANATYPE="Kaon")
        for key in fields:
            if key not in {"LTANAPATH", "OUTPATH", "USER", "HOST", "ANATYPE"}:
                Path(fields[key]).mkdir(parents=True, exist_ok=True)
        (root / "inputs/sentinel").write_bytes(b"synthetic scientific input: do not read or modify")
        config = package / "PATH_TO_DIR/user.path"
        config.write_text("".join(k + "=" + v + "\n" for k, v in fields.items()), encoding="utf-8")
        if duplicate:
            shutil.copyfile(config, config.with_name("second.path"))
        before = tool.tree_bytes(package)
        ordinary = tool.repository_state(repo)
        env = dict(os.environ, PYTHONPATH=str(installed), PYTHONDONTWRITEBYTECODE="1", GIT_OPTIONAL_LOCKS="0")
        with patch.dict(os.environ, env, clear=True), patch.object(tool, "BASELINE", sha):
            yield root, repo, package, sha
        # Assertions after both success/failure paths, with no repair/restoration.
        assert tool.tree_bytes(package) == before
        assert tool.repository_state(repo) == ordinary
        assert (root / "inputs/sentinel").read_bytes() == b"synthetic scientific input: do not read or modify"
        assert not list(root.glob("kaonlt-e8-4-canonical-five-*/"))


class AttestationTests(unittest.TestCase):
    def test_equal_default_and_plot_resolver_only_pass(self):
        with fixture() as (_, repo, _, sha):
            receipt = tool.attest(repo, sha)
            self.assertTrue(receipt["passed"], receipt)
            self.assertFalse(receipt["farm_execution_authorized"])
            self.assertEqual(len(receipt["callers"]), 7)
            self.assertEqual([r["mode"] for r in receipt["callers"][:3]], [[], [], []])
            self.assertEqual(receipt["preservation"], {"copied_package": True,
                "installed_package": True, "ordinary_checkout": True})
            for row in receipt["callers"]:
                self.assertEqual(row["baseline"]["paths"]["UTILPATH"], row["overlay"]["paths"]["UTILPATH"])
                self.assertTrue(row["input_path_identity"]["UTILPATH"]["read_only"])
            self.assertTrue(all(r["owner"] == "attempt" for r in receipt["resolved_targets"]))
            self.assertIn("complete launcher write confinement", receipt["unknown_operations"])
            self.assertFalse(list(repo.rglob("table_*.csv")))

    def test_exact_directory_caller_and_inherited_cwd_for_each_role(self):
        behavior = (
            '  expected=Path(self.LTANAPATH)\n'
            '  if caller.endswith("main.py"): expected=expected/"src"\n'
            '  elif caller.endswith(".py"): expected=expected/"src/setup"\n'
            '  if Path.cwd()!=expected: raise ValueError("wrong inherited cwd")\n'
            '  if caller.endswith("set_SymLinks.sh"): self.OUTPATH=str(Path(self.VOLATILEPATH)/"wrong-file-caller")\n')
        with fixture(behavior) as (_, repo, _, sha):
            receipt = tool.attest(repo, sha)
            self.assertTrue(receipt["passed"], receipt)
            self.assertIn("set_SymLinks.sh", receipt["source_sha256"])
            by_role = {row["role"]: row for row in receipt["callers"]}
            for role, caller, mode, cwd, source in tool.CALLERS:
                row = by_role[role]
                self.assertEqual(row["source_file"], source)
                for phase in ("baseline", "overlay"):
                    tree = Path(row[phase]["paths"]["LTANAPATH"])
                    self.assertEqual(row[phase]["cwd"], str(tree / cwd))
                    self.assertEqual(row[phase]["caller"], str(tree / caller))
            self.assertEqual(by_role["symlink_setup"]["baseline"]["caller"], str(repo))
            self.assertNotEqual(by_role["symlink_setup"]["role"], by_role["launcher"]["role"])
            wrong = tool.probe(repo, "set_SymLinks.sh", ("Plot_LTSep",), dict(os.environ))
            self.assertNotEqual(wrong["paths"]["OUTPATH"], by_role["symlink_setup"]["baseline"]["paths"]["OUTPATH"])
            # Reintroducing the reviewed wrong file argument must fail rather
            # than certify that alternate resolution as the symlink-setup role.
            callers = (*tool.CALLERS[:-1], ("symlink_setup", "set_SymLinks.sh", ("Plot_LTSep",), ".", "set_SymLinks.sh"))
            with patch.object(tool, "CALLERS", callers):
                result = tool.attest(repo, sha)
            self.assertFalse(result["passed"])
            self.assertEqual(result["failed_invariant"], "plot_root_outpath_external")

    def test_default_ignores_ltanapath_despite_good_plot(self):
        behavior = '  if not args: self.OUTPATH=str(Path(self.VOLATILEPATH)/"shared")\n'
        with fixture(behavior) as (_, repo, _, sha):
            result = tool.attest(repo, sha)
            self.assertEqual(result["failed_invariant"], "default_root_outpath_external")
            self.assertFalse(result["passed"])
            self.assertEqual(result["preservation"]["ordinary_checkout"], True)

    def test_default_differs_inside_attempt(self):
        with fixture('  if not args: self.OUTPATH=str(Path(self.LTANAPATH)/"different")\n') as (_, repo, _, sha):
            self.assertEqual(tool.attest(repo, sha)["failed_invariant"], "default_root_outpath_differs_from_plot")

    def test_default_callers_disagree(self):
        behavior = '  if not args and caller.endswith("getDataTable.py"): self.OUTPATH=str(Path(self.LTANAPATH)/"different")\n'
        with fixture(behavior) as (_, repo, _, sha):
            self.assertEqual(tool.attest(repo, sha)["failed_invariant"], "default_root_caller_mismatch")

    def test_no_mode_missing_outpath(self):
        with fixture('  if not args: del self.OUTPATH\n') as (_, repo, _, sha):
            self.assertIn("default_root_unresolved:src/setup/findEffectiveCharge.py:OUTPATH",
                          tool.attest(repo, sha)["failed_invariant"])

    def test_multiple_config_matches(self):
        with fixture(duplicate=True) as (_, repo, _, sha):
            self.assertEqual(tool.attest(repo, sha)["failed_invariant"], "ltsep_matching_config_not_unique")

    def test_changed_utilpath_authority_is_blocked(self):
        behavior = ('  if not args and "ordinary" not in self.LTANAPATH and caller.endswith("getDataTable.py"):\n'
                    '   self.UTILPATH=str(Path(self.LTANAPATH)/"wrong-input")\n')
        with fixture(behavior) as (_, repo, _, sha):
            self.assertEqual(tool.attest(repo, sha)["failed_invariant"],
                             "scientific_input_path_changed:src/setup/getDataTable.py:UTILPATH")

    def test_input_caller_disagreement_is_not_hidden_by_stability(self):
        behavior = '  if caller.endswith("getEfficiencyValue.py"): self.UTILPATH=str(Path(self.UTILPATH)/"different")\n'
        with fixture(behavior) as (_, repo, _, sha):
            self.assertEqual(tool.attest(repo, sha)["failed_invariant"], "default_root_input_caller_mismatch:UTILPATH")

    def test_readonly_import_and_resolver_reject_mutation(self):
        behaviors = [
            '  Path(__file__).write_text("corrupt installed package")\n',
            '  Path(__file__).with_name("PATH_TO_DIR").joinpath("user.path").write_text("corrupt config")\n',
            '  Path(self.UTILPATH,"forbidden.csv").write_text("not allowed")\n',
            '  Path(self.UTILPATH,"new-dir").mkdir()\n',
            '  __import__("subprocess").run(["echo","forbidden"])\n',
        ]
        for behavior in behaviors:
            with self.subTest(behavior=behavior), fixture(behavior) as (_, repo, _, sha):
                result = tool.attest(repo, sha)
                self.assertFalse(result["passed"])
                self.assertIn("resolver_subprocess_failed:", result["failed_invariant"])
                self.assertIn("resolver_", result["failed_invariant"])

    def test_resolver_cannot_read_scientific_data_or_import_science(self):
        for behavior in ('  Path(self.UTILPATH).parent.joinpath("sentinel").read_bytes()\n',
                         '  __import__("numpy")\n', '  __import__("ROOT")\n'):
            with self.subTest(behavior=behavior), fixture(behavior) as (_, repo, _, sha):
                result = tool.attest(repo, sha)
                self.assertFalse(result["passed"])
                self.assertIn("resolver_subprocess_failed:", result["failed_invariant"])

    def test_installed_fingerprint_mutation_is_detected(self):
        with fixture() as (_, repo, package, sha):
            original = tool.tree_bytes
            count = 0
            def changed(path):
                nonlocal count
                result = original(path)
                if path == package:
                    count += 1
                    if count > 1:
                        result["PATH_TO_DIR/user.path"] = "synthetic concurrent mismatch"
                return result
            with patch.object(tool, "tree_bytes", side_effect=changed):
                result = tool.attest(repo, sha)
            self.assertEqual(result["failed_invariant"], "installed_package_changed")
            self.assertFalse(result["preservation"]["installed_package"])

    def test_ordinary_state_mutation_is_detected_without_restoration(self):
        with fixture() as (_, repo, _, sha):
            original = tool.repository_state
            count = 0
            def changed(path):
                nonlocal count
                result = original(path)
                count += 1
                if count > 1:
                    result["head"] = "synthetic concurrent mismatch"
                return result
            with patch.object(tool, "repository_state", side_effect=changed):
                result = tool.attest(repo, sha)
            self.assertEqual(result["failed_invariant"], "ordinary_checkout_changed")
            self.assertFalse(result["preservation"]["ordinary_checkout"])

    def test_import_failure_and_stderr_and_stdout_contamination(self):
        for behavior, expected in [
            ('  raise ImportError("synthetic missing module")\n', "resolver_subprocess_failed:"),
            ('  print("noise",file=__import__("sys").stderr)\n', "resolver_stderr_contamination:"),
            ('  print("not JSON")\n', "resolver_json_invalid:"),
        ]:
            with self.subTest(expected=expected), fixture(behavior) as (_, repo, _, sha):
                self.assertIn(expected, tool.attest(repo, sha)["failed_invariant"])

    def test_alias_resolution_and_broken_cycles_without_writes(self):
        with fixture() as (root, repo, _, sha):
            original = tool.overlay
            for kind in ("external", "dangling", "cycle"):
                def install(*args):
                    result = original(*args)
                    worktree = args[1]
                    target = root / "inputs" if kind == "external" else root / "absent"
                    if kind == "cycle":
                        target = worktree / "OUTPUT"
                    (worktree / "OUTPUT").symlink_to(target, target_is_directory=True)
                    return result
                with self.subTest(kind=kind), patch.object(tool, "overlay", side_effect=install):
                    result = tool.attest(repo, sha)
                    self.assertFalse(result["passed"])
                    self.assertIn("outpath_external" if kind == "external" else "default_root_unresolved", result["failed_invariant"])

    def test_wrong_copy_import_and_wrong_caller(self):
        with fixture() as (_, repo, _, sha):
            original = tool.probe
            def wrong(*args):
                result = original(*args)
                if args[0] != repo:
                    result["package_file"] = str(repo / "__init__.py")
                return result
            with patch.object(tool, "probe", side_effect=wrong):
                self.assertFalse(tool.attest(repo, sha)["passed"])
            with self.assertRaisesRegex(ValueError, "caller_path_invalid"):
                tool.probe(repo, "src/setup/wrong.py", (), dict(os.environ))

    def test_copy_mutation_failure_preserves_installed_and_ordinary(self):
        with fixture() as (_, repo, _, sha):
            original = tool.probe
            def mutate(*args):
                result = original(*args)
                if args[0] != repo:
                    Path(result["pathing_file"]).write_text("synthetic overlay corruption")
                    raise ValueError("synthetic_probe_failure")
                return result
            with patch.object(tool, "probe", side_effect=mutate):
                result = tool.attest(repo, sha)
            self.assertFalse(result["passed"])
            self.assertEqual(result["preservation"]["installed_package"], True)
            self.assertEqual(result["preservation"]["ordinary_checkout"], True)
            self.assertEqual(result["failed_invariant"], "synthetic_probe_failure")
            self.assertIn("copied_package_changed_during_probe", result["additional_failures"])

    def test_source_head_index_worktree_drift_rejected(self):
        with fixture() as (_, repo, _, sha):
            self.assertIn("source_identity_mismatch", tool.attest(repo, "a" * 40)["failed_invariant"])
            target = repo / "src/setup/findEffectiveCharge.py"
            original = target.read_bytes()
            original_index = (repo / ".git/index").read_bytes()
            target.write_bytes(original + b"\n# synthetic drift\n")
            try:
                self.assertIn("source_worktree_changed", tool.attest(repo, sha)["failed_invariant"])
                run(repo, "git", "add", "--", str(target))
                self.assertIn("source_index_changed", tool.attest(repo, sha)["failed_invariant"])
            finally:
                # Only reset this exclusively fixture-owned index/file, never ordinary user state.
                target.write_bytes(original)
                (repo / ".git/index").write_bytes(original_index)

    def test_cli_exclusive_receipt_and_reentry(self):
        with fixture() as (_, repo, _, sha), tempfile.TemporaryDirectory(prefix="kaonlt-safe-receipts-") as output:
            out = Path(output).resolve()
            stdout = io.StringIO()
            with patch.object(tool.uuid, "uuid4", return_value=type("Id", (), {"hex": "fixed"})()), redirect_stdout(stdout):
                status = tool.main(["--repo", str(repo), "--source-commit", sha, "--output-dir", str(out)])
                self.assertEqual(status, 0, stdout.getvalue())
                first = list(out.iterdir())[0].read_bytes()
                with redirect_stderr(io.StringIO()):
                    self.assertEqual(tool.main(["--repo", str(repo), "--source-commit", sha, "--output-dir", str(out)]), 1)
                self.assertEqual(list(out.iterdir())[0].read_bytes(), first)
            result = json.loads(stdout.getvalue())
            self.assertEqual(result["sha256"], hashlib.sha256(first).hexdigest())
            self.assertFalse(json.loads(first)["farm_execution_authorized"])

    def test_cli_failure_receipt_is_blocked(self):
        with fixture('  if not args: self.OUTPATH=str(Path(self.VOLATILEPATH)/"shared")\n') as (_, repo, _, sha), tempfile.TemporaryDirectory(prefix="kaonlt-safe-receipts-") as output:
            out = Path(output).resolve()
            with redirect_stdout(io.StringIO()):
                self.assertEqual(tool.main(["--repo", str(repo), "--source-commit", sha, "--output-dir", str(out)]), 1)
            receipt = json.loads(next(out.iterdir()).read_bytes())
            self.assertEqual(receipt["status"], "BLOCKED")
            self.assertEqual(receipt["failed_invariant"], "default_root_outpath_external")

    def test_receipt_directory_hostile_paths(self):
        with fixture() as (root, repo, _, sha):
            out = root / "receipts"; out.mkdir(mode=0o700)
            link = root / "alias"; link.symlink_to(out, target_is_directory=True)
            for path in (link, repo, root / "absent", out / ".." / "receipts"):
                with self.subTest(path=path), redirect_stderr(io.StringIO()):
                    self.assertEqual(tool.main(["--repo", str(repo), "--source-commit", sha, "--output-dir", str(path)]), 1)
            self.assertFalse(list(out.iterdir()))

    def test_cli_never_publishes_inside_package_or_protected_inputs(self):
        with fixture() as (root, repo, package, sha):
            sentinel = (root / "inputs/sentinel").read_bytes()
            protected = [package, package / "PATH_TO_DIR", root / "inputs/UTILPATH",
                         root / "inputs/SIMCPATH", root / "inputs/SKIMPATH",
                         root / "inputs/REPLAYPATH", root / "inputs/VOLATILEPATH"]
            for path in protected:
                before_package = tool.tree_bytes(package)
                before_repo = tool.repository_state(repo)
                before_contents = sorted(path.iterdir())
                stderr = io.StringIO()
                with self.subTest(path=path), redirect_stderr(stderr), redirect_stdout(io.StringIO()):
                    status = tool.main(["--repo", str(repo), "--source-commit", sha, "--output-dir", str(path)])
                    self.assertEqual(status, 1)
                    self.assertIn("receipt_directory_protected_input", stderr.getvalue())
                self.assertEqual(tool.tree_bytes(package), before_package)
                self.assertEqual(tool.repository_state(repo), before_repo)
                self.assertEqual(sorted(path.iterdir()), before_contents)
                self.assertEqual((root / "inputs/sentinel").read_bytes(), sentinel)

    def test_resolved_protected_alias_and_external_background_inputs(self):
        with fixture() as (root, repo, package, sha), tempfile.TemporaryDirectory(prefix="kaonlt-protected-alias-") as output:
            outside = Path(output).resolve()
            alias = root / "input-alias"
            alias.symlink_to(outside, target_is_directory=True)
            for env in ({"LT_BG_SAMPLE_BASE": str(alias)},
                        {"LT_BG_SIGMA0_CENTER_LOW_ROOT": str(alias / "input.root")},
                        {"BACKGROUND_SIMC_PATH": str(alias)}):
                stderr = io.StringIO()
                before = tool.tree_bytes(package)
                with self.subTest(env=env), patch.dict(os.environ, env), redirect_stderr(stderr), redirect_stdout(io.StringIO()):
                    self.assertEqual(tool.main(["--repo", str(repo), "--source-commit", sha,
                                                "--output-dir", str(outside)]), 1)
                self.assertIn("receipt_directory_protected_input", stderr.getvalue())
                self.assertFalse(list(outside.iterdir()))
                self.assertEqual(tool.tree_bytes(package), before)
            # A receipt-directory alias itself is rejected before probing.
            with redirect_stderr(io.StringIO()), patch.object(tool, "attest") as attest:
                self.assertEqual(tool.main(["--repo", str(repo), "--source-commit", sha,
                                            "--output-dir", str(alias)]), 1)
                attest.assert_not_called()

    def test_publication_requires_complete_protected_discovery_and_rechecks(self):
        with fixture() as (_, repo, _, sha), tempfile.TemporaryDirectory(prefix="kaonlt-safe-receipts-") as output:
            out = Path(output).resolve()
            with patch.object(tool, "probe", side_effect=ValueError("synthetic_first_failure")), redirect_stderr(io.StringIO()) as stderr:
                self.assertEqual(tool.main(["--repo", str(repo), "--source-commit", sha,
                                            "--output-dir", str(out)]), 1)
            failure = json.loads(stderr.getvalue())
            self.assertEqual(failure["failed_invariant"], "synthetic_first_failure")
            self.assertEqual(failure["receipt_publication_failure"], "receipt_protected_paths_not_verified")
            self.assertFalse(list(out.iterdir()))

    def test_subprocess_failure_and_nonfile_source(self):
        with fixture() as (_, repo, _, sha):
            with patch.object(tool, "probe", side_effect=subprocess.TimeoutExpired("synthetic", 60)):
                result = tool.attest(repo, sha)
            self.assertFalse(result["passed"])
            self.assertTrue(result["preservation"]["ordinary_checkout"])


if __name__ == "__main__":
    unittest.main()
