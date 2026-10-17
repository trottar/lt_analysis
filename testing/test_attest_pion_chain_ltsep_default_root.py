"""Synthetic packages/repositories only; never attest the real JLab package."""
from contextlib import contextmanager, ExitStack, redirect_stdout, redirect_stderr
from copy import deepcopy
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from testing import attest_pion_chain_ltsep_default_root as tool
from testing import pion_chain_attempt_output_isolation as isolation

ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = "3ba464006090c2cea14978e31d39f715c9f95fa6"
HISTORY_HELPERS = (tool.TOOL, tool.TEST,
                   "testing/pion_chain_attempt_output_isolation.py",
                   "testing/test_pion_chain_attempt_output_isolation.py")


def run(repo, *args):
    return subprocess.run(args, cwd=repo, check=True, capture_output=True).stdout


def history_commit(repo, *paths):
    """Mutate only the newly created synthetic repository, never ROOT."""
    run(repo, "git", "add", "--", *paths)
    run(repo, "git", "-c", "user.name=fixture", "-c", "user.email=fixture@example.invalid",
        "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null",
        "commit", "-qm", "synthetic history successor")
    sha = tool.git(repo, "rev-parse", "HEAD").decode().strip()
    run(repo, "git", "update-ref", "refs/remotes/origin/test", sha)
    return sha


@contextmanager
def history_fixture(published=False):
    """Real Git history; optionally descend from the exact published commit.

    The published fixture borrows objects read-only. Its refs, index, commits
    and checkout live entirely in TemporaryDirectory; no fetch or parent
    checkout/index/ref mutation is performed.
    """
    with tempfile.TemporaryDirectory(prefix="kaonlt-receipt-history-") as name, ExitStack() as stack:
        repo = Path(name)
        run(repo, "git", "init", "-q", "-b", "test")
        if published:
            objects = Path(tool.git(ROOT, "rev-parse", "--git-path", "objects").decode().strip())
            if not objects.is_absolute():
                objects = ROOT / objects
            (repo / ".git/objects/info/alternates").write_text(str(objects.resolve()) + "\n")
            run(repo, "git", "update-ref", "refs/heads/test", PUBLISHED)
            run(repo, "git", "read-tree", PUBLISHED)
            names = run(repo, "git", "ls-files", "-z", "--", "src", "farm_env", "testing",
                        "background_samples", "run_Prod_Analysis.sh", "set_SymLinks.sh")
            subprocess.run(["git", "checkout-index", "-z", "--stdin"], cwd=repo,
                           input=names, check=True, capture_output=True)
            pin = PUBLISHED
        else:
            for source in set((*tool.FROZEN, *isolation.SOURCE_PATHS, *HISTORY_HELPERS)):
                target = repo / source
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes((ROOT / source).read_bytes())
            pin = history_commit(repo, ".")
            stack.enter_context(patch.object(tool, "BASELINE", pin))
            stack.enter_context(patch.object(isolation, "SOURCE_HEAD", pin))
        yield repo, pin


class HistoryCompatibilityTests(unittest.TestCase):
    def test_committed_isolation_helper_and_test_successors(self):
        with history_fixture() as (repo, _):
            expected = tool.source_identity(repo, tool.git(repo, "rev-parse", "HEAD").decode().strip())
            for name in HISTORY_HELPERS[2:]:
                path = repo / name
                path.write_bytes(path.read_bytes() + b"\n# synthetic isolation-only successor\n")
                sha = history_commit(repo, name)
                before = tool.repository_state(repo)
                self.assertEqual(tool.source_identity(repo, sha), expected)
                isolation.checked_source(repo)
                self.assertEqual(tool.repository_state(repo), before)

    def test_combined_successor_from_published_baseline(self):
        with history_fixture(published=True) as (repo, _):
            for name in HISTORY_HELPERS:
                (repo / name).write_bytes((ROOT / name).read_bytes())
            sha = history_commit(repo, *HISTORY_HELPERS)
            self.assertNotEqual(sha, PUBLISHED)
            run(repo, "git", "merge-base", "--is-ancestor", PUBLISHED, sha)
            changed = set(tool.git(repo, "diff", "--name-only", PUBLISHED, sha).decode().splitlines())
            self.assertEqual(changed, set(HISTORY_HELPERS))
            before = tool.repository_state(repo)
            identities = tool.source_identity(repo, sha)
            _, protected_hashes, head = isolation.checked_source(repo)
            self.assertEqual(head, sha)
            for name in isolation.SOURCE_PATHS:
                self.assertEqual(protected_hashes[name], tool.digest(tool.git(repo, "show", isolation.SOURCE_HEAD + ":" + name)))
            self.assertTrue(set(tool.FROZEN).issubset(identities))
            attempt = repo / "attempt"; attempt.mkdir()
            receipt = isolation.preflight(repo, attempt, {
                "LTANAPATH": str(attempt / "analysis"), "OUTPATH": str(attempt / "analysis/OUTPUT"),
                "SKIMPATH": str(attempt / "skim"), "VOLATILEPATH": str(repo.parent / "shared"),
                "ANATYPE": "Kaon"})
            self.assertEqual(receipt["failed_invariant"], isolation.FIRST_SHARED_WRITE)
            self.assertFalse(receipt["passed"])
            self.assertFalse(receipt["farm_execution_authorized"])
            self.assertEqual(tool.repository_state(repo), before)

    def test_unauthorized_committed_sources_remain_rejected_by_both(self):
        for name in ("src/unapproved.py", "set_SymLinks.sh", "farm_env/ltsep_paths.py",
                     "testing/unrelated.py"):
            with self.subTest(name=name), history_fixture() as (repo, _):
                path = repo / name
                path.write_bytes((path.read_bytes() if path.exists() else b"") + b"\n# forbidden\n")
                sha = history_commit(repo, name)
                with self.assertRaisesRegex(ValueError, "unapproved_source_successor"):
                    tool.source_identity(repo, sha)
                with self.assertRaisesRegex(ValueError, "protected_source_changed|unapproved_source_change"):
                    isolation.checked_source(repo)

    def test_protected_worktree_and_index_remain_rejected_by_both(self):
        for indexed in (False, True):
            with self.subTest(indexed=indexed), history_fixture() as (repo, pin):
                path = repo / "set_SymLinks.sh"
                original = path.read_bytes()
                path.write_bytes(original + b"\n# forbidden\n")
                if indexed:
                    run(repo, "git", "add", "--", "set_SymLinks.sh")
                    path.write_bytes(original)  # Dirty index alone cannot hide.
                with self.assertRaisesRegex(ValueError, "source_index_changed" if indexed else "source_worktree_changed"):
                    tool.source_identity(repo, pin)
                with self.assertRaisesRegex(ValueError, "protected_source_changed"):
                    isolation.checked_source(repo)

    def test_wrong_ancestry_remains_rejected_by_both(self):
        with history_fixture() as (repo, _):
            run(repo, "git", "checkout", "--orphan", "unrelated")
            (repo / "docs").mkdir()
            (repo / "docs/orphan.md").write_text("unrelated history\n")
            sha = history_commit(repo, ".")
            run(repo, "git", "branch", "-M", "test")
            with self.assertRaises(subprocess.CalledProcessError):
                tool.source_identity(repo, sha)
            with self.assertRaisesRegex(ValueError, "scientific_source_not_verified_ancestor"):
                isolation.checked_source(repo)


@contextmanager
def fixture(behavior="", pathing_extra="", duplicate=False, volatile=None):
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
        if volatile is not None:
            fields["VOLATILEPATH"] = str(volatile)
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
    @contextmanager
    def canonical_destination(self):
        with tempfile.TemporaryDirectory(prefix="kaonlt-canonical-storage-") as name:
            root = Path(name).resolve()
            physical = root / "physical"
            physical.mkdir()
            out = physical / "globus"
            out.mkdir()
            out.chmod(0o2775)
            alias = root / "infrastructure"
            alias.symlink_to(physical, target_is_directory=True)
            logical = alias / "globus"
            original = Path.lstat

            def facility_lstat(path, *args, **kwargs):
                info = original(path, *args, **kwargs)
                if path == alias:
                    # Analogue of facility-owned infrastructure: no actual
                    # privileged ownership/permission changes in the fixture.
                    values = list(info)
                    values[4] = 0
                    return os.stat_result(values)
                return info

            with patch.object(tool.primitives.accepted, "GLOBUS", logical), patch.object(Path, "lstat", facility_lstat):
                yield root, logical, out, alias

    def test_canonical_infrastructure_and_2775_publication(self):
        with self.canonical_destination() as (storage, logical, physical, _):
            metadata = physical.stat()
            sentinel = physical / "existing.zip"
            sentinel.write_bytes(b"preserved bundle")
            with fixture(volatile=storage) as (_, repo, _, sha):
                # One resolver attestation per CLI, including publication.
                for destination in (logical, physical):
                    stdout = io.StringIO()
                    with patch.object(tool, "attest", wraps=tool.attest) as attest, redirect_stdout(stdout):
                        self.assertEqual(tool.main(["--repo", str(repo), "--source-commit", sha,
                                                    "--output-dir", str(destination)]), 0)
                    attest.assert_called_once()
                    summary = json.loads(stdout.getvalue())
                    path = Path(summary["receipt"])
                    raw = path.read_bytes()
                    receipt = json.loads(raw)
                    self.assertEqual(summary["sha256"], hashlib.sha256(raw).hexdigest())
                    self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o600)
                    self.assertFalse(receipt["farm_execution_authorized"])
                    policy = receipt["receipt_destination"]
                    self.assertEqual(policy["role"], "canonical_globus_bundle_root")
                    self.assertEqual(policy["resolved_target"], str(physical))
                    self.assertEqual(policy["identity"], [metadata.st_dev, metadata.st_ino])
                    self.assertEqual(policy["broad_storage_exceptions"], [str(storage)])
                    self.assertIn(str(storage / "OUTPUT"), policy["protected_realpaths"])
            self.assertEqual(sentinel.read_bytes(), b"preserved bundle")
            after = physical.stat()
            self.assertEqual((after.st_dev, after.st_ino, after.st_mode, after.st_uid, after.st_gid),
                             (metadata.st_dev, metadata.st_ino, metadata.st_mode, metadata.st_uid, metadata.st_gid))

    def test_canonical_alias_negative_controls(self):
        with self.canonical_destination() as (root, logical, out, alias):
            other = root / "other"; other.mkdir(); other.chmod(0o2775)
            alternate = root / "alternate"; alternate.symlink_to(out.parent, target_is_directory=True)
            leaf = root / "leaf"; leaf.symlink_to(out, target_is_directory=True)
            dangling = root / "dangling"; dangling.symlink_to(root / "absent", target_is_directory=True)
            cycle = root / "cycle"; cycle.symlink_to(cycle, target_is_directory=True)
            for path in (other, alternate / "globus", leaf, dangling / "globus", cycle / "globus",
                         logical / ".." / "globus", root / "missing"):
                with self.subTest(path=path), redirect_stderr(io.StringIO()), patch.object(tool, "attest") as attest:
                    self.assertEqual(tool.main(["--repo", str(ROOT), "--source-commit", "a" * 40,
                                                "--output-dir", str(path)]), 1)
                    attest.assert_not_called()
            self.assertFalse(list(out.iterdir()))
            # An unprivileged alias is not an approved facility ancestor.
            original = Path.lstat
            def untrusted(path, *args, **kwargs):
                info = original(path, *args, **kwargs)
                if path == alias:
                    values = list(info); values[4] = os.getuid() + 1
                    return os.stat_result(values)
                return info
            with patch.object(Path, "lstat", untrusted), self.assertRaisesRegex(ValueError, "alias_not_root_owned"):
                tool.output_directory(logical, ROOT)

    def test_canonical_protected_roles_never_inherit_storage_exception(self):
        with self.canonical_destination() as (storage, logical, out, _), fixture(volatile=storage) as (_, repo, _, sha):
            receipt = tool.attest(repo, sha)
            self.assertTrue(receipt["passed"])
            receipt["receipt_destination"] = tool.output_directory(logical, repo)
            for field in tool.FIELDS:
                if field in {"VOLATILEPATH", "USER", "HOST", "ANATYPE"}:
                    continue
                candidate = deepcopy(receipt)
                candidate["baseline_callers"][0]["paths"][field] = str(out)
                with self.subTest(field=field):
                    protected = tool.protected_receipt_locations(repo, candidate, out)
                    with self.assertRaisesRegex(ValueError, "protected_input"):
                        tool.output_directory(logical, repo, protected)
            for role in ("installed_package", "overlay_package"):
                candidate = deepcopy(receipt)
                candidate[role] = {"path": str(out)}
                protected = tool.protected_receipt_locations(repo, candidate, out)
                with self.subTest(role=role), self.assertRaisesRegex(ValueError, "protected_input"):
                    tool.output_directory(logical, repo, protected)
            # An actual background input under the canonical root blocks it,
            # including alias resolution and ancestor/descendant overlap.
            link = storage / "input-alias"; link.symlink_to(out, target_is_directory=True)
            candidate = deepcopy(receipt)
            candidate["configured_external_inputs"] = {"LT_BG_SAMPLE_BASE": str(link / "science")}
            with self.assertRaisesRegex(ValueError, "protected_input"):
                tool.output_directory(logical, repo, tool.protected_receipt_locations(repo, candidate, out))
            sibling = storage / "unrelated"; sibling.mkdir(mode=0o700)
            with self.assertRaisesRegex(ValueError, "protected_input"):
                tool.output_directory(sibling, repo, tool.protected_receipt_locations(repo, deepcopy(receipt)))
            self.assertFalse(list(out.iterdir()))

    def test_canonical_owner_inode_and_alias_rechecks(self):
        with self.canonical_destination() as (_, logical, out, alias):
            policy = tool.output_directory(logical, ROOT)
            original = Path.stat
            def wrong_owner(path, *args, **kwargs):
                info = original(path, *args, **kwargs)
                if path == logical:
                    values = list(info); values[4] = os.getuid() + 1
                    return os.stat_result(values)
                return info
            with patch.object(Path, "stat", wrong_owner), self.assertRaisesRegex(ValueError, "not_owned"):
                tool.output_directory(logical, ROOT)
            out.rename(out.with_name("previous"))
            out.mkdir(); out.chmod(0o2775)
            with self.assertRaisesRegex(ValueError, "directory_replaced"):
                tool.output_directory(logical, ROOT, expected=policy)
            policy = tool.output_directory(logical, ROOT)
            alias.unlink(); alias.symlink_to(str(out.parent) + "/.", target_is_directory=True)
            with self.assertRaisesRegex(ValueError, "directory_replaced"):
                tool.output_directory(logical, ROOT, expected=policy)

    def test_canonical_race_and_interruption_publish_no_receipt(self):
        with self.canonical_destination() as (storage, logical, out, _), fixture(volatile=storage) as (_, repo, _, sha):
            original = tool.attest
            def replace(*args):
                receipt = original(*args)
                out.rename(out.with_name("previous"))
                out.mkdir(); out.chmod(0o2775)
                return receipt
            with patch.object(tool, "attest", side_effect=replace), redirect_stderr(io.StringIO()) as stderr:
                self.assertEqual(tool.main(["--repo", str(repo), "--source-commit", sha,
                                            "--output-dir", str(logical)]), 1)
            self.assertIn("receipt_directory_replaced", stderr.getvalue())
            self.assertFalse(list(out.iterdir()))
            self.assertFalse(list(out.with_name("previous").iterdir()))
            with patch.object(tool, "attest", side_effect=KeyboardInterrupt), self.assertRaises(KeyboardInterrupt):
                tool.main(["--repo", str(repo), "--source-commit", sha, "--output-dir", str(logical)])
            self.assertFalse(list(out.iterdir()))

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
