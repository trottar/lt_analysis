"""Exercise the complete owner locally with synthetic artifacts and no farm."""
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from copy import deepcopy
import ast
import hashlib
import io
import json
import os
import re
import sys
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch
import zipfile

from testing import run_e8_4_fix5_left_lowe_plot_gate as owner

ROOT = Path(__file__).resolve().parents[1]
SHA = "a" * 40

for relative_path in ("src/cuts", "src/utility"):
    path = str(ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

import full_background_subtraction_plots as plots
import pion_hgcer_refinement_checkpoint as checkpoint


def producer_setting():
    return checkpoint.build_pion_hgcer_refinement_checkpoint(
        setting={"kinematic_token": owner.KINEMATIC, "Q2": 4.4, "W": 2.74,
                 "epsilon_setting": "low", "epsilon_filename_token": "lowe",
                 "phi_setting": "Left", "particle_type": "kaon"},
        phase_a={}, method_a={}, method_b={},
    )["setting"]


def pages():
    records = [{"page_id": page_id, "scope": scope} for page_id, scope in owner.CRITICAL_PAGES]
    for page_id in owner.NEW_PAGE_IDS:
        if page_id.endswith("parent_closure"):
            records.append({"page_id": page_id, "scope": "setting"})
        else:
            records.append({"page_id": page_id, "scope": page_id[-2:], "t_index": int(page_id[-1]) - 1,
                            "represented_phi_inventory": [{"phi_index": i, "phi_edges": [-180 + i * 40, -140 + i * 40]} for i in range(9)]})
    return plots.build_full_background_subtraction_page_manifest_artifact(
        setting=producer_setting(), pdf_basename=owner.PDF,
        pages=records, renderer_failures=[],
    )


def git_fixture(repo, *args):
    if args == ("branch", "--show-current"):
        return "test" if repo == ROOT else ""
    if args[0] == "rev-parse":
        return SHA
    return ""


class OwnerTests(unittest.TestCase):
    def test_current_checkpoint_to_manifest_to_owner_preserves_provenance(self):
        data = pages()
        before = deepcopy(data)
        self.assertEqual(data["setting"], producer_setting())
        self.assertEqual(data["setting"]["Q2"], 4.4)
        self.assertEqual(data["setting"]["W"], 2.74)
        self.assertEqual(data["setting"]["epsilon_setting"], "low")
        self.assertEqual(owner.verify_pages(data), len(data["pages"]))
        self.assertEqual(data, before)

    def test_wrong_or_missing_setting_identity_fails_closed(self):
        for key, wrong in (("kinematic_token", "Q3p0W2p32"),
                           ("epsilon_filename_token", "highe"), ("phi_setting", "Right"),
                           ("particle_type", "pion"), ("epsilon_setting", "high")):
            for missing in (False, True):
                data = pages()
                if missing:
                    del data["setting"][key]
                else:
                    data["setting"][key] = wrong
                before = deepcopy(data)
                with self.subTest(key=key, missing=missing), self.assertRaisesRegex(
                        ValueError, "^page_manifest_setting_invalid$"):
                    owner.verify_pages(data)
                self.assertEqual(data, before)
        for invalid in (None, [], "low", 4):
            data = pages(); data["setting"] = invalid
            with self.subTest(setting=invalid), self.assertRaisesRegex(
                    ValueError, "^page_manifest_setting_invalid$"):
                owner.verify_pages(data)

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
            (lambda data: data.update(schema_version="wrong"), "page_manifest_schema_invalid"),
            (lambda data: data["pages"].pop(), "new_page_missing_or_duplicate"),
            (lambda data: data["pages"].append(deepcopy(data["pages"][-1])), "new_page_missing_or_duplicate"),
            (lambda data: data.update(renderer_failures=["failed"]), "renderer_failures_nonempty"),
            (lambda data: data["pages"].pop(0), "critical_page_missing_or_duplicate"),
            (lambda data: data.update(pdf_basename="wrong.pdf"), "page_manifest_pdf_invalid"),
            (lambda data: data["pages"][-2].update(represented_phi_inventory=[]), "new_page_child_inventory_invalid"),
            (lambda data: data["pages"][-2].update(scope="setting"), "new_page_scope_invalid"),
            (lambda data: data["pages"][-2].update(t_index=99), "new_page_t_identity_invalid"),
            (lambda data: data["pages"][-1].update(scope="t1"), "parent_closure_scope_invalid"),
        ):
            data = pages(); mutation(data)
            with self.subTest(reason=reason), self.assertRaisesRegex(ValueError, reason):
                owner.verify_pages(data)

    def exercise(self, directory, *, fail_run=False, stale=False, fail_collection=False,
                 fail_preflight=False, missing_markers=False, fail_zip=False):
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
            status = json.loads(owner.gate_status_path(outdir, output).read_text())
            self.assertEqual((status["status"], status["stage"]), ("running", "analysis"))
            self.assertTrue(status["collector_source_preflight_completed"])
            self.assertTrue(status["analysis_started"])
            self.assertFalse(status["analysis_completed"])
            self.assertFalse(output.exists())
            events.append("run")
            log_path.write_text("Left / lowe debug analysis completed.\nFull high-epsilon processing is intentionally skipped in -d debug mode.\n")
            if fail_run:
                return 1
            if missing_markers:
                log_path.write_text("analysis returned zero without completion markers\n")
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
        real_checks = owner.collector.collect_source_checks
        real_load_profile = owner.collector.load_validation_profile
        def load_profile(path):
            if Path(path).as_posix().endswith(owner.PROFILE):
                self.assertNotEqual(Path(path).parents[1], ROOT)
                return real_load_profile(ROOT / owner.PROFILE)
            return real_load_profile(path)
        def source_checks(repo, command_runner=None, **kwargs):
            self.assertNotEqual(repo, ROOT)
            self.assertEqual(kwargs["required_analysis_commit"], SHA)
            self.assertEqual(kwargs["allowed_committed_files"], [])
            if command_runner is None:
                events.append("source_preflight")
                def command_runner(command, cwd=None):
                    self.assertEqual(cwd, repo)
                    failed = fail_preflight and "py_compile" in command
                    return {"command": command, "returncode": 7 if failed else 0,
                            "stdout": "", "stderr": "literal syntax failure" if failed else ""}
            else:
                events.append("source_recheck")
            return real_checks(repo, command_runner, **kwargs)
        real_collect = owner.collector.collect_validation_bundle
        @contextmanager
        def module(worktree):
            self.assertNotEqual(worktree, ROOT)
            yield owner.collector
        with patch.object(owner, "git", side_effect=git), patch.object(owner, "collection_module", side_effect=module), patch.object(owner, "resolved_profile", return_value=profile), patch.object(owner, "CANDIDATES", candidates), patch.object(owner, "run_analysis", side_effect=run), patch.object(owner.collector, "collect_validation_bundle", side_effect=collect), patch.object(owner.collector, "load_validation_profile", side_effect=load_profile), patch.object(owner.collector, "collect_source_checks", side_effect=source_checks), patch.object(owner, "verify_zip", side_effect=ValueError("literal zip failure") if fail_zip else owner.verify_zip):
            # The template is read from the real repository; no source operation
            # or production launcher is performed by this synthetic gate.
            try:
                result = owner.execute_gate(ROOT, outdir, output, SHA)
            finally:
                expected_events = ["source_preflight"]
                if not fail_preflight:
                    expected_events += ["run"]
                    if not (fail_run or stale or missing_markers):
                        expected_events += ["collect"]
                        if not fail_collection:
                            expected_events += ["source_recheck"]
                self.assertEqual(events, expected_events)
                additions = [args for repo, args in git_calls if args[:2] == ("worktree", "add")]
                removals = [args for repo, args in git_calls if args[:2] == ("worktree", "remove")]
                self.assertEqual(len(additions), 1 if fail_preflight or fail_run or stale or missing_markers else 2)
                self.assertEqual(len(removals), len(additions))
                if additions:
                    for addition, removal in zip(additions, removals):
                        self.assertEqual(addition[2], "--detach")
                        self.assertEqual(addition[-1], SHA)
                        self.assertEqual(removal, ("worktree", "remove", "--force", addition[-2]))
        return result, events, outdir

    def test_complete_run_verify_real_generic_collection_and_zip_integrity(self):
        with tempfile.TemporaryDirectory() as directory:
            result, events, outdir = self.exercise(directory)
            self.assertEqual(events, ["source_preflight", "run", "collect", "source_recheck"])
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
            status = json.loads(owner.gate_status_path(outdir, result).read_text())
            self.assertEqual((status["status"], status["stage"], status["failure_reason"]), ("success", "complete", None))
            for flag in ("analysis_started", "analysis_completed", "artifact_verification_completed",
                         "collector_source_preflight_completed", "collection_completed", "zip_verification_completed"):
                self.assertIs(status[flag], True)
            self.assertEqual(status["source_commit"], SHA)
            self.assertEqual(status["expected_zip_path"], result.as_posix())
            self.assertEqual(status["analysis_log_path"], (outdir / "fresh.log").as_posix())
            self.assertLessEqual(status["started_at_utc"], status["updated_at_utc"])
            with zipfile.ZipFile(result) as archive:
                self.assertFalse(any("gate-status" in name for name in archive.namelist()))
            summary = json.loads((outdir / owner.SUMMARY).read_text())
            self.assertEqual(summary["pre_run_provenance"]["head"], SHA)
            self.assertEqual(summary["run_log_sha256"], owner.sha256(outdir / "fresh.log"))

    def test_failed_analysis_stale_pdf_collection_failure_stop(self):
        for kwargs, reason in (({"fail_run": True}, "analysis_failed"), ({"stale": True}, "artifact_stale"), ({"fail_collection": True}, "collection_failed")):
            with self.subTest(reason=reason), tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(ValueError, reason):
                    self.exercise(directory, **kwargs)

    def test_each_failure_persists_literal_stage_and_completed_flags(self):
        for kwargs, stage, reason, completed, artifacts, collection in (
            ({"fail_preflight": True}, "collector_source_preflight", "collector_source_check_failed:py_compile:returncode=7:literal syntax failure", False, False, False),
            ({"fail_run": True}, "analysis", "analysis_failed", False, False, False),
            ({"missing_markers": True}, "completion_markers", "left_lowe_completion_missing", True, False, False),
            ({"stale": True}, "verify_artifacts", "artifact_stale:" + owner.PDF, True, False, False),
            ({"fail_collection": True}, "collection", "collection_failed", True, True, False),
            ({"fail_zip": True}, "verify_zip", "literal zip failure", True, True, True),
        ):
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as directory:
                with self.assertRaises(ValueError) as caught:
                    self.exercise(directory, **kwargs)
                self.assertEqual(str(caught.exception), reason)
                root = Path(directory)
                status = json.loads((root / "artifacts/fresh-gate-status.json").read_text())
                self.assertEqual((status["status"], status["stage"], status["failure_reason"]), ("failed", stage, reason))
                self.assertEqual(status["analysis_completed"], completed)
                self.assertEqual(status["artifact_verification_completed"], artifacts)
                self.assertEqual(status["collection_completed"], collection)
                self.assertFalse(status["zip_verification_completed"])
                self.assertEqual(status["analysis_started"], stage != "collector_source_preflight")
                self.assertEqual(status["analysis_log_path"], (root / "artifacts/fresh.log").as_posix())
                if stage == "collector_source_preflight":
                    self.assertFalse((root / "fresh.zip").exists())
                    self.assertFalse((root / "artifacts/fresh.log").exists())

    def test_detached_preflight_rejects_effective_identity_profile_and_paths(self):
        profile = owner.resolved_profile(ROOT, SHA)
        for change, reason in (("profile", "collector_preflight_profile_mismatch"),
                               ("commit", "collector_preflight_required_commit_mismatch"),
                               ("paths", "collector_preflight_unexpected_committed_files")):
            effective = deepcopy(profile)
            if change == "profile":
                effective["validation_profile"] = "wrong"
            if change == "commit":
                effective["source_identity"]["required_analysis_commit"] = "b" * 40
            detached = deepcopy(profile)
            checks = [{"name": "required_analysis_commit_ancestor", "returncode": 0, "stdout": "", "stderr": ""},
                      {"name": "committed_files_after_required_analysis_commit", "returncode": 0,
                       "stdout": "src/unexpected.py" if change == "paths" else "", "stderr": ""}]
            with self.subTest(change=change), patch.object(owner, "git", side_effect=git_fixture), patch.object(owner, "collection_module") as module:
                module.return_value.__enter__.return_value = owner.collector
                with patch.object(owner.collector, "load_validation_profile", return_value=detached), patch.object(owner.collector, "collect_source_checks", return_value=("checks", checks)), self.assertRaisesRegex(ValueError, reason):
                    owner.collector_source_preflight(ROOT, SHA, effective)

    def test_status_publication_is_atomic_and_refuses_existing_attempt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); output = root / "fresh.zip"
            with patch.object(owner.os, "link", wraps=owner.os.link) as link, patch.object(owner.os, "replace", wraps=owner.os.replace) as replace:
                status = owner.GateStatus(root, output, SHA)
                link.assert_called_once()
                self.assertEqual(replace.call_count, 0)
                status.update("profile")
                self.assertEqual(replace.call_count, 1)
            before = status.path.read_bytes()
            with patch.object(owner, "preflight") as preflight, self.assertRaises(FileExistsError):
                owner.execute_gate(ROOT, root, output, SHA)
            self.assertEqual(status.path.read_bytes(), before)
            preflight.assert_not_called()
            self.assertEqual(list(root.glob(".gate-status-*")), [])

    def test_unexpected_collector_import_failure_is_persisted_before_analysis(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(owner, "preflight", return_value={}), patch.object(owner, "resolved_profile", return_value={"artifacts": {"global": [], "settings": []}}), patch.object(owner, "CANDIDATES", {}), patch.object(owner, "collector_source_preflight", side_effect=ImportError("missing collector dependency")), patch.object(owner, "run_analysis") as run:
            root = Path(directory)
            with self.assertRaisesRegex(ImportError, "missing collector dependency"):
                owner.execute_gate(ROOT, root, root / "fresh.zip", SHA)
            run.assert_not_called()
            status = json.loads((root / "fresh-gate-status.json").read_text())
            self.assertEqual((status["stage"], status["failure_reason"]), ("collector_source_preflight", "missing collector dependency"))
            self.assertEqual(status["status"], "failed")

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
            stdout = io.StringIO(); stderr = io.StringIO()
            with patch.object(owner, "execute_gate", side_effect=ValueError("failed") if error else None, return_value=Path("/fresh.zip")), redirect_stdout(stdout), redirect_stderr(stderr):
                result = owner.main(["--source-commit", SHA])
            self.assertEqual(result, 1 if error else 0)
            self.assertEqual(stdout.getvalue(), "" if error else "/fresh.zip\n")
            if error:
                self.assertIn("Owner gate status: ", stderr.getvalue())
                self.assertIn("-gate-status.json", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
