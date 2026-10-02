"""Exercise the complete owner locally with synthetic artifacts and no farm."""
from contextlib import contextmanager, redirect_stdout
from copy import deepcopy
import ast
import hashlib
import io
import json
import os
import re
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch
import zipfile

from testing import run_e8_4_fix5_left_lowe_plot_gate as owner

ROOT = Path(__file__).resolve().parents[1]
SHA = "a" * 40


def pages():
    records = [{"page_id": page_id, "scope": scope} for page_id, scope in owner.CRITICAL_PAGES]
    for page_id in owner.NEW_PAGE_IDS:
        if page_id.endswith("parent_closure"):
            records.append({"page_id": page_id, "scope": "setting"})
        else:
            records.append({"page_id": page_id, "scope": page_id[-2:], "t_index": int(page_id[-1]) - 1,
                            "represented_phi_inventory": [{"phi_index": i, "phi_edges": [-180 + i * 40, -140 + i * 40]} for i in range(9)]})
    return {"schema_version": "full_background_subtraction_page_manifest/v1",
            "setting": {"kinematic_token": owner.KINEMATIC, "epsilon_filename_token": "lowe", "phi_setting": "Left", "particle_type": "kaon"},
            "pdf_basename": owner.PDF, "pages": records, "renderer_failures": []}


def git_fixture(repo, *args):
    if args == ("branch", "--show-current"):
        return "test" if repo == ROOT else ""
    if args[0] == "rev-parse":
        return SHA
    return ""


class OwnerTests(unittest.TestCase):
    def test_real_owner_invokes_debug_launcher_with_non_destructive_cleanup_guard(self):
        tree = ast.parse((ROOT / "testing/run_e8_4_fix5_left_lowe_plot_gate.py").read_text(encoding="utf-8"))
        run = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "run_analysis")
        calls = [node for node in ast.walk(run) if isinstance(node, ast.Call)
                 and isinstance(node.func, ast.Attribute) and node.func.attr == "Popen"]
        self.assertEqual(len(calls), 1)
        self.assertEqual(ast.literal_eval(calls[0].args[0]), ["./run_Prod_Analysis.sh", "-d", "4p4", "2p74"])
        self.assertEqual(next(keyword.value.id for keyword in calls[0].keywords if keyword.arg == "cwd"), "repo")
        launcher = (ROOT / "run_Prod_Analysis.sh").read_text(encoding="utf-8")
        self.assertEqual(len(re.findall(r'^\s*git clean -fdx\s*$', launcher, re.MULTILINE)), 1)
        self.assertRegex(launcher, r'if ! validate_external_sigma0_paths_before_cleanup; then\s+exit 1\s+fi\s+'
                         r'if \[\[ \$d_flag != "true" \]\]; then\s+git clean -fdx\s+fi\s+'
                         r'\./set_SymLinks\.sh \$ParticleType')

    def test_wrong_head_branch_push_and_dirty_source_block_before_run(self):
        for change, reason in (("head", "wrong_head"), ("branch", "wrong_branch"),
                               ("push", "source_not_observed_pushed"), ("dirty", "dirty_gate_source"),
                               ("untracked", "dirty_gate_source")):
            def git(repo, *args):
                if change == "head" and args == ("rev-parse", "HEAD"):
                    return "b" * 40
                if change == "branch" and args == ("branch", "--show-current"):
                    return "other"
                if change == "push" and args == ("rev-parse", "refs/remotes/origin/test"):
                    return "b" * 40
                if change == "dirty" and args[0] == "status":
                    return " M src/cuts/rand_sub.py"
                if change == "untracked" and args[0] == "status":
                    return "?? testing/unreviewed.py"
                return git_fixture(repo, *args)
            with self.subTest(change=change), tempfile.TemporaryDirectory() as directory, patch.object(owner, "git", side_effect=git), patch.object(owner, "run_analysis") as run:
                with self.assertRaisesRegex(ValueError, reason):
                    owner.execute_gate(ROOT, Path(directory), Path(directory) / "fresh.zip", SHA)
                run.assert_not_called()

    def test_known_farm_outputs_are_recorded_without_cleanup(self):
        status = " M src/models/xmodel_kaon_pl.f\n?? src/kaon/functions/Q4p4W2p74.model"
        calls = []
        def git(repo, *args):
            calls.append((repo, args))
            if args[0] == "diff":
                self.assertIn("HEAD", args)
                self.assertEqual(set(args[-2:]), {":(exclude)" + path for path in owner.FARM_OUTPUTS})
                self.assertEqual(args[args.index("--") + 1:-2], owner.GATE_PATHS)
            return status if args[0] == "status" else git_fixture(repo, *args)
        with patch.object(owner, "git", side_effect=git):
            result = owner.preflight(ROOT, SHA)
        self.assertEqual(result["status_short"], status.splitlines())
        self.assertTrue(all(args[0] not in {"clean", "reset", "checkout", "stash", "worktree"} for _, args in calls))

    def test_new_ids_duplicates_failures_and_old_pages_fail_closed(self):
        for mutation, reason in (
            (lambda data: data["pages"].pop(), "new_page_missing_or_duplicate"),
            (lambda data: data["pages"].append(deepcopy(data["pages"][-1])), "new_page_missing_or_duplicate"),
            (lambda data: data.update(renderer_failures=["failed"]), "renderer_failures_nonempty"),
            (lambda data: data["pages"].pop(0), "critical_page_missing_or_duplicate"),
            (lambda data: data.update(pdf_basename="wrong.pdf"), "page_manifest_pdf_invalid"),
            (lambda data: data["pages"][-2].update(represented_phi_inventory=[]), "new_page_child_inventory_invalid"),
            (lambda data: data["pages"][-2].update(scope="setting"), "new_page_scope_invalid"),
            (lambda data: data["pages"][-1].update(scope="t1"), "parent_closure_scope_invalid"),
        ):
            data = pages(); mutation(data)
            with self.subTest(reason=reason), self.assertRaisesRegex(ValueError, reason):
                owner.verify_pages(data)

    def exercise(self, directory, *, fail_run=False, stale=False, fail_collection=False):
        root = Path(directory); outdir = root / "artifacts"; outdir.mkdir()
        output = root / "fresh.zip"
        profile = owner.resolved_profile(ROOT, SHA)
        candidates = {}
        for name in owner.CANDIDATES:
            (outdir / name).write_text('{"fixture": true}', encoding="utf-8")
            candidates[name] = owner.sha256(outdir / name)
        events = []; git_calls = []
        def git(repo, *args):
            git_calls.append((repo, args))
            return git_fixture(repo, *args)
        def run(repo, log_path):
            self.assertEqual(repo, ROOT)
            events.append("run")
            log_path.write_text("Left / lowe debug analysis completed.\nFull high-epsilon processing is intentionally skipped in -d debug mode.\n")
            if fail_run:
                return 1
            for scope in ("global", "settings"):
                for declaration in profile["artifacts"][scope]:
                    name = declaration["basename_template"].format(kinematic=owner.KINEMATIC, phi="Left", epsilon="lowe")
                    if name in candidates or name == owner.SUMMARY:
                        continue
                    content = json.dumps(pages()).encode() if name == owner.PAGE_MANIFEST else b'%PDF-fixture' if name == owner.PDF else b'{}' if declaration["kind"] == "json" else b'fixture,csv\n'
                    (outdir / name).write_bytes(content)
                    if stale and name == owner.PDF:
                        os.utime(outdir / name, (1.0, 1.0))
            return 0
        def collect(**kwargs):
            events.append("collect")
            self.assertEqual((kwargs["phi"], kwargs["epsilon"]), ("Left", "lowe"))
            self.assertNotEqual(kwargs["repo_root"], ROOT)
            self.assertEqual(kwargs["repo_root"].name, "source")
            self.assertEqual(kwargs["outdir"], outdir)
            effective = json.loads(Path(kwargs["profile_path"]).read_text())
            self.assertEqual(effective, profile)
            if fail_collection:
                return {"returncode": 1}
            def command_runner(command, cwd=None):
                self.assertEqual(cwd, kwargs["repo_root"])
                stdout = SHA + "\n" if command[:3] == ["git", "rev-parse", "HEAD"] else ""
                return {"command": command, "returncode": 0, "stdout": stdout, "stderr": ""}
            return real_collect(**kwargs, command_runner=command_runner)
        real_collect = owner.collector.collect_validation_bundle
        @contextmanager
        def module(worktree):
            self.assertNotEqual(worktree, ROOT)
            yield owner.collector
        with patch.object(owner, "git", side_effect=git), patch.object(owner, "collection_module", side_effect=module), patch.object(owner, "resolved_profile", return_value=profile), patch.object(owner, "CANDIDATES", candidates), patch.object(owner, "run_analysis", side_effect=run), patch.object(owner.collector, "collect_validation_bundle", side_effect=collect):
            # The template is read from the real repository; no source operation
            # or production launcher is performed by this synthetic gate.
            try:
                result = owner.execute_gate(ROOT, outdir, output, SHA)
            finally:
                self.assertEqual(events, ["run"] if fail_run or stale else ["run", "collect"])
                additions = [args for repo, args in git_calls if args[:2] == ("worktree", "add")]
                removals = [args for repo, args in git_calls if args[:2] == ("worktree", "remove")]
                self.assertEqual(len(additions), 0 if fail_run or stale else 1)
                self.assertEqual(len(removals), len(additions))
                if additions:
                    self.assertEqual(additions[0][2], "--detach")
                    self.assertEqual(additions[0][-1], SHA)
                    self.assertEqual(removals[0], ("worktree", "remove", "--force", additions[0][-2]))
        return result, events, outdir

    def test_complete_run_verify_real_generic_collection_and_zip_integrity(self):
        with tempfile.TemporaryDirectory() as directory:
            result, events, outdir = self.exercise(directory)
            self.assertEqual(events, ["run", "collect"])
            self.assertTrue(result.is_file())
            with zipfile.ZipFile(result) as archive:
                manifest = json.loads(archive.read("manifest.json"))
                self.assertTrue(manifest["complete"])
                self.assertEqual(manifest["required_analysis_commit"], SHA)
                self.assertEqual(manifest["requested_settings"], [{"phi": "Left", "epsilon": "lowe"}])
                self.assertEqual(len(manifest["settings"]), 1)
                self.assertEqual(set(manifest["global_artifacts"]), {"candidate_f3", "candidate_f4", "run_summary"})
                self.assertEqual(set(manifest["settings"][0]["artifacts"]),
                                 {"procedure_pdf", "page_manifest", "full_analysis", "correction_ledger_json", "correction_ledger_csv"})
            summary = json.loads((outdir / owner.SUMMARY).read_text())
            self.assertEqual(summary["pre_run_provenance"]["head"], SHA)

    def test_failed_analysis_stale_pdf_collection_failure_stop(self):
        for kwargs, reason in (({"fail_run": True}, "analysis_failed"), ({"stale": True}, "artifact_stale"), ({"fail_collection": True}, "collection_failed")):
            with self.subTest(reason=reason), tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(ValueError, reason):
                    self.exercise(directory, **kwargs)

    def test_existing_zip_refused_before_preflight_or_run(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(owner, "preflight") as preflight:
            output = Path(directory) / "existing.zip"; output.write_bytes(b"existing")
            with self.assertRaisesRegex(ValueError, "output_zip_already_exists"):
                owner.execute_gate(ROOT, Path(directory), output, SHA)
            preflight.assert_not_called()

    def test_detached_source_module_is_loaded_from_its_own_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            worktree = Path(directory); testing = worktree / "testing"; testing.mkdir()
            for name in ("collect_pion_hgcer_validation_bundle.py", "pion_hgcer_validation_bundle_profile.json"):
                (testing / name).write_bytes((ROOT / "testing" / name).read_bytes())
            with owner.collection_module(worktree) as module:
                self.assertEqual(Path(module.__file__), testing / "collect_pion_hgcer_validation_bundle.py")
                self.assertEqual(module.resolve_settings("Left", "lowe"), (("Left", "lowe"),))

    def test_worktree_identity_failure_still_removes_only_owned_tree(self):
        for failure, reason in (("head", "head_mismatch"), ("branch", "not_detached"), ("status", "dirty")):
            calls = []
            def git(repo, *args):
                calls.append((repo, args))
                if args[0] == "rev-parse": return "b" * 40 if failure == "head" else SHA
                if args[0] == "branch": return "test" if failure == "branch" else ""
                if args[0] == "status": return " M src/main.py" if failure == "status" else ""
                return ""
            with self.subTest(failure=failure), patch.object(owner, "git", side_effect=git):
                with self.assertRaisesRegex(ValueError, reason), owner.clean_collection_worktree(ROOT, SHA):
                    self.fail("invalid worktree reached collection")
            self.assertEqual(calls[0][1][:3], ("worktree", "add", "--detach"))
            self.assertEqual(calls[-1], (ROOT, ("worktree", "remove", "--force", calls[0][1][-2])))

    def test_main_prints_exactly_one_zip_only_after_success_and_none_on_failure(self):
        for error in (False, True):
            stdout = io.StringIO()
            with patch.object(owner, "execute_gate", side_effect=ValueError("failed") if error else None, return_value=Path("/fresh.zip")), redirect_stdout(stdout):
                result = owner.main(["--source-commit", SHA])
            self.assertEqual(result, 1 if error else 0)
            self.assertEqual(stdout.getvalue(), "" if error else "/fresh.zip\n")


if __name__ == "__main__":
    unittest.main()
