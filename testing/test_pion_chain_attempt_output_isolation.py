"""Bounded first-blocker reproduction; no real launcher, ROOT, network or farm.

All mutation-bearing Git/link operations below use synthetic repositories and
fake ltsep trees inside a newly created, user-owned TemporaryDirectory.
Later full-launcher write coverage is deliberately not claimed after BLOCKED.
"""
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from testing import pion_chain_attempt_output_isolation as isolation
from testing import run_e8_4_fix5_canonical_five_plot_gate as reviewed


ROOT = Path(__file__).resolve().parents[1]


def run(tree, *args, env=None):
    return subprocess.run(args, cwd=tree, env=env, check=True,
                          capture_output=True, text=True).stdout


def files(tree):
    return {str(p.relative_to(tree)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in tree.rglob("*") if p.is_file() and not p.is_symlink()}


def runtime_fixture(directory):
    """Local copy of the reviewed fake-ltsep fixture, with identical lookup.

    This intentionally imports no existing test suite or analysis module.
    Root derives OUTPATH from fixed LTANAPATH; changing caller alone does not
    change that field. SetPath selects by ANALYSISPATH and expands ${USER}.
    """
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
    fields = {k: str(root / (k.lower() + "-fixture")) for k in reviewed.PATH_FIELDS}
    fields.update(LTANAPATH=str(repo), ANALYSISPATH=str(root), USER="resolved-user",
                  HOST="synthetic-host", ANATYPE="Kaon")
    config = package / "PATH_TO_DIR/user.path"
    config.write_text("".join(k + "=" + v + "\n" for k, v in fields.items()), encoding="utf-8")
    return repo, worktree, package, config


def commit_fixture(repo, message):
    """Commit only within a newly created synthetic repository."""
    run(repo, "git", "add", "--all")
    run(repo, "git", "-c", "user.name=fixture", "-c", "user.email=fixture@example.invalid",
        "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null",
        "commit", "--quiet", "-m", message)
    return isolation.git(repo, "rev-parse", "HEAD").decode().strip()


@contextmanager
def source_history():
    """Real synthetic commits with the exact pinned protected file contents."""
    with tempfile.TemporaryDirectory(prefix="kaonlt-source-history-test-") as directory:
        repo = Path(directory)
        run(repo, "git", "init", "--quiet")
        for name in isolation.SOURCE_PATHS:
            path = repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(isolation.git(ROOT, "show", f"{isolation.SOURCE_HEAD}:{name}"))
        pin = commit_fixture(repo, "synthetic protected baseline")
        # Substitute only this synthetic history's anchor; all ancestry and
        # content checks still execute real Git, without mocked results.
        with patch.object(isolation, "SOURCE_HEAD", pin):
            yield repo, pin


@contextmanager
def sandbox():
    with tempfile.TemporaryDirectory(prefix="kaonlt-chain-isolation-test-") as directory:
        base = Path(directory).resolve()
        assert base.stat().st_uid == os.getuid()
        repo, tree, package, config = runtime_fixture(directory)
        values = reviewed.read_path_config(config.read_bytes())
        for role in ("VOLATILEPATH", "SIMCPATH", "SKIMPATH"):
            Path(values[role]).mkdir()
        volatile = Path(values["VOLATILEPATH"])
        for leaf in ("OUTPUT/Analysis/KaonLT", "log", "worksim"):
            (volatile / leaf).mkdir(parents=True)
        (Path(values["SKIMPATH"]) / "Analysis/KaonLT").mkdir(parents=True)
        background = base / "external-background"
        background.mkdir()
        for protected in (volatile, background, Path(values["SIMCPATH"]), package):
            (protected / "sentinel").write_bytes(b"immutable synthetic input\n")
        (repo / "background_samples").mkdir()
        (repo / "background_samples/background_samples.conf").write_text(
            'BACKGROUND_SIMC_PATH="' + str(background) + '"\n', encoding="utf-8")
        (repo / "input/kaon").mkdir(parents=True)
        (repo / "input/kaon/sentinel").write_bytes(b"synthetic input\n")
        (repo / "set_SymLinks.sh").write_bytes(isolation.git(ROOT, "show", "HEAD:set_SymLinks.sh"))
        (repo / "src").mkdir()
        (repo / "src/placeholder").write_text("tracked fixture\n")
        run(repo, "git", "init", "--quiet")
        run(repo, "git", "add", "--all")  # Entirely synthetic index.
        run(repo, "git", "-c", "user.name=fixture", "-c", "user.email=fixture@example.invalid",
            "-c", "commit.gpgsign=false", "commit", "--quiet", "-m", "synthetic")
        # Remove only the duplicate farm_env created by the reviewed fixture.
        assert tree.is_relative_to(base) and not tree.is_symlink()
        shutil.rmtree(tree)
        run(repo, "git", "worktree", "add", "--detach", str(tree), "HEAD")
        ambient = {"PATH": os.environ["PATH"], "HOME": str(base),
                   "PYTHONPATH": str(package.parent), "PYTHONDONTWRITEBYTECODE": "1",
                   "USER": "resolved-user"}
        with patch.dict(os.environ, ambient, clear=True):
            baseline = reviewed.path_fields(repo, repo)
            # A new caller still sees the original fixed LTANAPATH before the
            # overlay patches only that one copied configuration value.
            assert reviewed.path_fields(tree, tree)["LTANAPATH"] == str(repo)
            env, overlay = reviewed.prepare_runtime_overlay(repo, tree)
            roles = reviewed.probe_paths(tree, env, baseline)[0]["paths"]
        yield base, repo, tree, roles, env, package, overlay
        # TemporaryDirectory unlinks symlinks; it never traverses their targets.


@unittest.skipUnless(os.name == "posix", "requires local POSIX Git/symlink semantics")
class IsolationTests(unittest.TestCase):
    def test_source_order_and_no_analysis_import(self):
        texts, hashes, head = isolation.checked_source(ROOT)
        self.assertEqual(len(hashes), len(isolation.SOURCE_PATHS))
        self.assertEqual(head, isolation.git(ROOT, "rev-parse", "HEAD").decode().strip())
        self.assertIn('declare -a EPS=("low" "high")', texts["run_Prod_Analysis.sh"])
        self.assertNotIn("ROOT", isolation.__dict__)

    def test_later_memory_and_helper_test_commit_keeps_frozen_source_identity(self):
        with source_history() as (repo, pin):
            baseline_text, baseline_hashes, baseline_head = isolation.checked_source(repo)
            self.assertEqual(baseline_head, pin)
            doc = repo / "docs/memory/fixture.md"
            doc.parent.mkdir(parents=True)
            doc.write_text("synthetic documentation successor\n")
            (repo / "testing/test_pion_chain_attempt_output_isolation.py").write_text(
                "# synthetic test-only successor\n")
            later = commit_fixture(repo, "synthetic documentation and test successor")
            texts, hashes, head = isolation.checked_source(repo)
            self.assertNotEqual(later, pin)
            self.assertEqual(head, later)
            self.assertEqual(texts, baseline_text)
            self.assertEqual(hashes, baseline_hashes)
            attempt = repo / "attempt"; attempt.mkdir()
            receipt = isolation.preflight(repo, attempt, {})
            self.assertEqual(receipt["scientific_source_HEAD"], pin)
            self.assertEqual(receipt["checked_source_HEAD"], later)
            self.assertEqual(receipt["failed_invariant"], "required_path_role_missing")
            self.assertFalse(receipt["passed"])

    def test_changed_protected_source_rejected_in_worktree_index_and_commit(self):
        with source_history() as (repo, pin):
            name = "set_SymLinks.sh"
            path = repo / name
            path.write_bytes(path.read_bytes() + b"\n# unapproved change\n")
            with self.assertRaisesRegex(ValueError, "protected_source_changed:" + name):
                isolation.checked_source(repo)
            run(repo, "git", "add", "--", name)
            # Restore only the synthetic worktree bytes: a dirty index alone
            # must still be rejected, not hidden by a matching worktree.
            original = isolation.git(repo, "show", f"{pin}:{name}")
            changed = path.read_bytes()
            path.write_bytes(original)
            with self.assertRaisesRegex(ValueError, "protected_source_changed:" + name):
                isolation.checked_source(repo)
            path.write_bytes(changed)
            commit_fixture(repo, "synthetic unapproved protected change")
            with self.assertRaisesRegex(ValueError, "protected_source_changed:" + name):
                isolation.checked_source(repo)

    def test_unapproved_other_source_commit_is_rejected(self):
        with source_history() as (repo, _):
            path = repo / "src/unapproved.py"
            path.write_text("# synthetic unapproved source\n")
            commit_fixture(repo, "synthetic unapproved source addition")
            with self.assertRaisesRegex(ValueError, "unapproved_source_change:src/unapproved.py"):
                isolation.checked_source(repo)

    def test_unrelated_history_and_missing_pin_fail_ancestry(self):
        with source_history() as (repo, pin):
            run(repo, "git", "checkout", "--orphan", "unrelated")
            # An orphan commit can carry identical protected bytes but cannot
            # inherit scientific authority without the required ancestor.
            commit_fixture(repo, "synthetic unrelated history")
            with self.assertRaisesRegex(ValueError, "scientific_source_not_verified_ancestor"):
                isolation.checked_source(repo)
            with patch.object(isolation, "SOURCE_HEAD", "0" * 40):
                with self.assertRaisesRegex(ValueError, "scientific_source_not_verified_ancestor"):
                    isolation.checked_source(repo)

    def test_real_detached_clean_then_unchanged_setup_reproduces_blocker(self):
        with sandbox() as (base, repo, tree, roles, env, package, overlay):
            ordinary_before = (isolation.git(repo, "rev-parse", "HEAD"),
                               isolation.git(repo, "show-ref"), (repo / ".git/index").read_bytes(),
                               isolation.git(repo, "status", "--porcelain=v1"), files(repo))
            installed_before = files(package)
            volatile = Path(roles["VOLATILEPATH"])
            sentinels = {p: p.read_bytes() for p in base.rglob("sentinel")}
            alias_target = tree.parent / "candidate-output"
            alias_target.mkdir()
            (tree / "OUTPUT").symlink_to(alias_target, target_is_directory=True)
            before = os.readlink(tree / "OUTPUT")
            # The actual launcher mkdirs precede clean. This path is initially
            # alias-confined; setup subsequently recreates a shared OUTPUT link.
            Path(roles["OUTPATH"]).mkdir(parents=True)
            clean = run(tree, "git", "clean", "-fdx")
            self.assertIn("OUTPUT", clean)
            self.assertFalse((tree / "OUTPUT").is_symlink())
            self.assertTrue(alias_target.is_dir())
            bg_output = volatile / "OUTPUT/Analysis/KaonLT/background_samples"
            self.assertFalse(bg_output.exists())
            setup = run(tree, "bash", "./set_SymLinks.sh", "kaon", env=env)
            self.assertIn("now setup", setup)
            self.assertTrue(bg_output.is_dir())  # Forbidden real-role write reproduced only here.
            external_links = {}
            for name in ("OUTPUT", "log", "src/data", "simc_gfortran/OUTPUTS",
                         "simc_gfortran/input", "simc_gfortran/worksim",
                         "background_samples/OUTPUTS", "background_samples/simc_gfortran/worksim"):
                self.assertTrue((tree / name).is_symlink(), name)
                external_links[name] = str((tree / name).resolve(strict=True))
            for name in ("OUTPUT", "log", "src/data", "simc_gfortran/OUTPUTS",
                         "simc_gfortran/worksim", "background_samples/OUTPUTS",
                         "background_samples/simc_gfortran/worksim"):
                self.assertTrue(Path(external_links[name]).is_relative_to(volatile), name)
            self.assertTrue(Path(roles["OUTPATH"]).resolve().is_relative_to(volatile))
            # input's *link entry* was created in the external SIMC tree even
            # though its destination is inside the disposable worktree.
            self.assertTrue((Path(roles["SIMCPATH"]) / "input").is_symlink())
            self.assertEqual(files(package), installed_before)
            self.assertEqual({p: p.read_bytes() for p in sentinels}, sentinels)
            self.assertTrue(reviewed.verify_ltsep_preservation(overlay)["passed"])
            ordinary_after = (isolation.git(repo, "rev-parse", "HEAD"),
                              isolation.git(repo, "show-ref"), (repo / ".git/index").read_bytes(),
                              isolation.git(repo, "status", "--porcelain=v1"), files(repo))
            self.assertEqual(ordinary_before, ordinary_after)
            print(json.dumps({"fixture_only": True, "passed": False,
                "failed_invariant": isolation.FIRST_SHARED_WRITE,
                "farm_execution_authorized": False,
                "pre_clean_OUTPUT": "attempt/candidate-output",
                "post_clean_OUTPUT": "absent",
                "post_setup_OUTPUT": "external/VOLATILEPATH/OUTPUT",
                "untracked_alias_deleted": bool(before),
                "external_link_roles": sorted(external_links),
                "sentinels_preserved": True, "ordinary_git_preserved": True,
                "fake_installed_ltsep_preserved": True}, sort_keys=True))

    def test_read_only_receipt_stops_at_direct_background_mkdir(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            attempt = base / "attempt"; attempt.mkdir()
            skim = base / "skims/KaonLT"; skim.mkdir(parents=True)
            roles = {"LTANAPATH": str(attempt / "analysis"),
                     "OUTPATH": str(attempt / "analysis/OUTPUT/Analysis/KaonLT"),
                     "VOLATILEPATH": str(base / "shared"), "SKIMPATH": str(skim), "ANATYPE": "Kaon"}
            before = list(base.rglob("*"))
            receipt = isolation.preflight(ROOT, attempt, roles)
            self.assertEqual(receipt["failed_invariant"], isolation.FIRST_SHARED_WRITE)
            self.assertFalse(receipt["passed"])
            self.assertFalse(receipt["farm_execution_authorized"])
            self.assertEqual(before, list(base.rglob("*")))
            self.assertEqual(receipt["role_resolved_write_targets"][-1]["owner"], "external")
            self.assertEqual(receipt, json.loads(json.dumps(receipt)))
            self.assertEqual(receipt, isolation.preflight(ROOT, attempt, roles))

    def test_broken_cycle_relative_and_alias_escape(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            attempt = base / "attempt"; attempt.mkdir()
            (attempt / "escape").symlink_to(base, target_is_directory=True)
            record = isolation.target_record("alias", attempt / "escape/output", attempt, "write")
            self.assertEqual(record["owner"], "external")
            (attempt / "broken").symlink_to(base / "absent")
            (attempt / "cycle").symlink_to(attempt / "cycle")
            for leaf in ("broken", "cycle"):
                with self.assertRaisesRegex(ValueError, "broken_or_cyclic"):
                    isolation.target_record(leaf, attempt / leaf, attempt, "write")
            with self.assertRaisesRegex(ValueError, "relative_path"):
                isolation.target_record("relative", Path("relative"), attempt, "write")

    def test_collision_stale_and_no_debug_shortcut(self):
        with tempfile.TemporaryDirectory() as directory:
            attempt = Path(directory)
            for command in (("./run_Prod_Analysis.sh", "-d", "4p4", "2p74"),
                            ("./run_Prod_Analysis.sh", "-i1", "4p4", "2p74")):
                receipt = isolation.preflight(ROOT, attempt, {}, command)
                self.assertEqual(receipt["failed_invariant"], "normal_launcher_command_required_no_debug_shortcut")
            (attempt / "stale").write_text("keep")
            for _ in range(2):
                receipt = isolation.preflight(ROOT, attempt, {})
                self.assertEqual(receipt["failed_invariant"], "attempt_root_collision_or_stale")
                self.assertEqual((attempt / "stale").read_text(), "keep")

    def test_permission_and_source_errors_fail_closed(self):
        with patch.object(isolation, "checked_source", side_effect=PermissionError("fixture denial")):
            receipt = isolation.preflight(ROOT, ROOT, {})
        self.assertIn("preflight_indeterminate:fixture denial", receipt["failed_invariant"])
        with patch.object(isolation, "checked_source", side_effect=ValueError("source_head_mismatch")):
            receipt = isolation.preflight(ROOT, ROOT, {})
        self.assertIn("source_head_mismatch", receipt["failed_invariant"])


if __name__ == "__main__":
    unittest.main()
