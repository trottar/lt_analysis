#!/usr/bin/env python3
"""Own the reviewed Left/lowe run -> verify -> package gate.

After actual-diff review, user push and pushed-state synchronization, supply
that exact pushed SHA as --source-commit. The reviewed declarative profile is
instantiated with that SHA through the collector's existing profile_path API;
no tracked re-pin, collector change or additional source phase is needed.
All subprocess/log output goes to stderr. Success stdout is one ZIP path.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import contextmanager, redirect_stdout
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import zipfile

try:
    from . import collect_pion_hgcer_validation_bundle as collector
except ImportError:
    import collect_pion_hgcer_validation_bundle as collector

REPO = Path("/group/c-kaonlt/USERS/trottar/lt_analysis")
ARTIFACTS = REPO / "OUTPUT/Analysis/KaonLT"
GLOBUS = Path("/volatile/hallc/c-kaonlt/trottar/globus")
PROFILE = "testing/pion_hgcer_validation_bundle_profile_e8_4_fix5.json"
BASE_HEAD = "da38444e7aa60efd62d6638780776344daf40276"
KINEMATIC = "Q4p4W2p74"
PDF = "Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf"
PAGE_MANIFEST = PDF.replace(".pdf", "-manifest.json")
SUMMARY = "Q4p4W2p74_e8_4_fix5_left_lowe_run-summary.json"
NEW_PAGE_IDS = tuple(
    "full_background.e8_4.{}.t{}".format(group, index)
    for group in ("method_a_vs_simc", "baseline_method_a_simc", "yield_summary")
    for index in (1, 2, 3)
) + ("full_background.e8_4.parent_closure",)
OLD_PAGES = tuple(("full_background.e8_4." + group, "t" + str(index))
                  for group in ("pion_consequence", "final_mm", "signed_difference", "yield_impact")
                  for index in (1, 2, 3)) + (
    ("full_background.e8_4.authority", "setting"),
    ("full_background.e8_4.setting_summary", "setting"),
)
CRITICAL_PAGES = OLD_PAGES + tuple(
    ("full_background.e8_2.final_mm.t" + str(index), "t" + str(index)) for index in (1, 2, 3)
) + tuple(("full_background.e8_3." + group, "t" + str(index))
          for group in ("mm", "tphi") for index in (1, 2, 3))
FARM_OUTPUTS = {"src/models/xmodel_kaon_pl.f", "src/kaon/functions/Q4p4W2p74.model"}
GATE_PATHS = ("src", "testing", "tools", "run_Prod_Analysis.sh")
CANONICAL_SETTINGS = tuple({"phi": phi, "epsilon": epsilon} for phi, epsilon in (
    ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
    ("Center", "highe"), ("Right", "highe")))
CANDIDATES = {
    "Q4p4W2p74_kaon_pion-background_hgcer-method-a-acceptance-map-current-baseline-candidate.json":
        "eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d",
    "Q4p4W2p74_kaon_pion-background_hgcer-method-a-parent-preserving-correction-current-baseline-candidate.json":
        "1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902",
}


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(repo, *arguments):
    result = subprocess.run(["git", *arguments], cwd=repo, capture_output=True, text=True)
    require(result.returncode == 0, "git_failed:{}:{}".format(arguments, result.stderr))
    return result.stdout.rstrip("\r\n")


def preflight(repo, source_commit):
    require(re.fullmatch(r"[0-9a-f]{40}", source_commit) is not None, "source_commit_requires_full_sha")
    require(git(repo, "branch", "--show-current") == "test", "wrong_branch")
    require(git(repo, "rev-parse", "HEAD") == source_commit, "wrong_head")
    require(git(repo, "rev-parse", "refs/remotes/origin/test") == source_commit, "source_not_observed_pushed")
    require(source_commit != BASE_HEAD, "repair_source_not_committed")
    git(repo, "merge-base", "--is-ancestor", BASE_HEAD, source_commit)
    status = git(repo, "status", "--porcelain=v1", "--untracked-files=all")
    for line in status.splitlines():
        path = line[3:]
        # Only explicitly known farm model outputs and unrelated output/review
        # artifacts are harmless. Everything else, including untracked source,
        # is blocked; no checkout cleanup is ever performed.
        require(path in FARM_OUTPUTS or path.startswith("OUTPUT/")
                or (line[:2] == "??" and re.fullmatch(r"kaonlt_review\([^/]+\)\.diff", path)),
                "dirty_gate_source:" + line)
    git(repo, "diff", "--check", "HEAD", "--", *GATE_PATHS,
        *(":(exclude)" + path for path in sorted(FARM_OUTPUTS)))
    return {"branch": "test", "head": source_commit, "origin_test": source_commit,
            "status_short": status.splitlines(), "observed_at_utc": datetime.now(timezone.utc).isoformat()}


def resolved_profile(repo, source_commit):
    profile = collector.load_validation_profile(repo / PROFILE)
    require(profile["source_identity"]["required_analysis_commit"] == BASE_HEAD, "profile_candidate_base_mismatch")
    require(profile["settings"] == list(CANONICAL_SETTINGS), "profile_scope_mismatch")
    require(collector.resolve_settings("Left", "lowe", profile) == (("Left", "lowe"),), "profile_left_lowe_not_authorized")
    require(profile["source_identity"]["allowed_committed_files"] == [], "profile_source_allowlist_not_empty")
    profile["source_identity"]["required_analysis_commit"] = source_commit
    return profile


@contextmanager
def clean_collection_worktree(repo, source_commit):
    """Isolate collector source checks; removal prunes only this registration."""
    with tempfile.TemporaryDirectory(prefix="kaonlt-e8-4-fix5-collection-") as directory:
        parent = Path(directory).resolve()
        worktree = parent / "source"
        git(repo, "worktree", "add", "--detach", str(worktree), source_commit)
        try:
            require(git(worktree, "rev-parse", "HEAD") == source_commit, "collection_worktree_head_mismatch")
            require(git(worktree, "branch", "--show-current") == "", "collection_worktree_not_detached")
            require(git(worktree, "status", "--porcelain=v1", "--untracked-files=all") == "", "collection_worktree_dirty")
            yield worktree
        finally:
            require(worktree.resolve() == parent / "source" and parent != repo.resolve()
                    and worktree.resolve() != repo.resolve(), "collection_cleanup_path_invalid")
            # Source checks can create ignored bytecode. Force applies only to
            # the uniquely owned temporary tree, never the ordinary checkout.
            git(repo, "worktree", "remove", "--force", str(worktree))


@contextmanager
def collection_module(worktree):
    """Execute the unchanged collector from the exact detached source tree."""
    name = "_kaonlt_e8_4_fix5_detached_collector"
    spec = importlib.util.spec_from_file_location(name, worktree / "testing/collect_pion_hgcer_validation_bundle.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.modules.pop(name, None)


def run_analysis(repo, log_path):
    with log_path.open("xb") as log:
        process = subprocess.Popen(["./run_Prod_Analysis.sh", "-d", "4p4", "2p74"], cwd=repo,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(process.stdout.readline, b""):
            log.write(line); log.flush()
            sys.stderr.write(line.decode("utf-8", errors="replace")); sys.stderr.flush()
        return process.wait()


def verify_pages(payload):
    require(payload.get("schema_version") == "full_background_subtraction_page_manifest/v1", "page_manifest_schema_invalid")
    require(payload.get("setting") == {"kinematic_token": KINEMATIC, "epsilon_filename_token": "lowe", "phi_setting": "Left", "particle_type": "kaon"}, "page_manifest_setting_invalid")
    require(payload.get("pdf_basename") == PDF, "page_manifest_pdf_invalid")
    require(payload.get("renderer_failures") == [], "renderer_failures_nonempty")
    pages = payload.get("pages")
    require(isinstance(pages, list) and all(isinstance(page, dict) for page in pages), "page_inventory_invalid")
    ids = Counter(page.get("page_id") for page in pages)
    scopes = Counter((page.get("page_id"), page.get("scope")) for page in pages)
    require(all(ids[page_id] == 1 for page_id in NEW_PAGE_IDS), "new_page_missing_or_duplicate")
    require(all(scopes[key] == 1 for key in CRITICAL_PAGES), "critical_page_missing_or_duplicate")
    for page in pages:
        if page.get("page_id") == NEW_PAGE_IDS[-1]:
            require(page.get("scope") == "setting", "parent_closure_scope_invalid")
        elif page.get("page_id") in NEW_PAGE_IDS:
            index = int(page["page_id"][-1]) - 1
            require(page.get("scope") == "t" + str(index + 1), "new_page_scope_invalid")
            require(page.get("t_index") == index, "new_page_t_identity_invalid")
            require([child.get("phi_index") for child in page.get("represented_phi_inventory", [])] == list(range(9)), "new_page_child_inventory_invalid")
    return len(pages)


def verify_artifacts(outdir, profile, started_ns, before):
    records = {}
    for scope in ("global", "settings"):
        for declaration in profile["artifacts"][scope]:
            name = declaration["basename_template"].format(kinematic=KINEMATIC, phi="Left", epsilon="lowe")
            if name == SUMMARY:
                continue
            path = outdir / name
            require(path.is_file() and path.stat().st_size > 0, "artifact_missing:" + name)
            if name not in CANDIDATES:
                require(path.stat().st_mtime_ns >= started_ns, "artifact_stale:" + name)
                require((path.stat().st_mtime_ns, path.stat().st_size, sha256(path)) != before.get(name), "artifact_not_refreshed:" + name)
            if declaration["kind"] == "json":
                collector._strict_json_payload(path)
            if name in CANDIDATES:
                require(sha256(path) == CANDIDATES[name], "candidate_identity_mismatch:" + name)
            records[name] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    require((outdir / PDF).read_bytes().startswith(b"%PDF-"), "procedure_pdf_invalid")
    count = verify_pages(collector._strict_json_payload(outdir / PAGE_MANIFEST))
    return records, count


def verify_zip(path, source_commit, records):
    with zipfile.ZipFile(path) as archive:
        require(archive.testzip() is None, "zip_integrity_failed")
        manifest = json.loads(archive.read("manifest.json"))
        require(manifest.get("complete") is True and manifest.get("errors") == [], "bundle_incomplete")
        require(manifest.get("git_head") == source_commit and manifest.get("required_analysis_commit") == source_commit, "bundle_source_identity_mismatch")
        require(manifest.get("requested_settings") == [{"phi": "Left", "epsilon": "lowe"}], "bundle_scope_mismatch")
        require(manifest.get("validation_profile") == "phase_e8_4_fix5_left_lowe_shareable_pages/v1", "bundle_profile_mismatch")
        entries = list(manifest.get("global_artifacts", {}).values())
        require(len(manifest.get("settings", [])) == 1, "bundle_setting_inventory_invalid")
        entries += list(manifest["settings"][0].get("artifacts", {}).values())
        seen = set()
        for entry in entries:
            name = Path(entry["archive_path"]).name
            require(name in records and name not in seen, "bundle_artifact_inventory_invalid")
            raw = archive.read(entry["archive_path"])
            digest = hashlib.sha256(raw).hexdigest()
            require(digest == entry["sha256"] == records[name]["sha256"] and len(raw) == entry["byte_size"] == records[name]["bytes"], "bundle_artifact_hash_mismatch:" + name)
            seen.add(name)
        require(seen == set(records), "bundle_artifact_missing")


def gate_status_path(outdir, output):
    return outdir / (output.stem + "-gate-status.json")


class GateStatus:
    """Atomically publish this attempt's diagnostic, separate from the bundle."""
    def __init__(self, outdir, output, source_commit):
        self.path = gate_status_path(outdir, output)
        now = datetime.now(timezone.utc).isoformat()
        self.record = {
            "schema_version": "e8_4_fix5_owner_gate_status/v1",
            "source_commit": source_commit, "kinematic": KINEMATIC,
            "phi": "Left", "epsilon": "lowe",
            "expected_zip_path": output.as_posix(),
            "analysis_log_path": (outdir / (output.stem + ".log")).as_posix(),
            "started_at_utc": now, "updated_at_utc": now,
            "status": "running", "stage": "preflight", "failure_reason": None,
            "analysis_started": False, "analysis_completed": False,
            "artifact_verification_completed": False,
            "collector_source_preflight_completed": False,
            "collection_completed": False, "zip_verification_completed": False,
        }
        self._persist(initial=True)

    def _persist(self, *, initial=False):
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=self.path.parent,
                                             prefix=".gate-status-", delete=False) as handle:
                temporary = Path(handle.name)
                handle.write(json.dumps(self.record, indent=2, sort_keys=True) + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            if initial:
                # Exclusive atomic publication: an existing attempt is never replaced.
                os.link(temporary, self.path)
            else:
                os.replace(temporary, self.path)
        finally:
            if temporary is not None and temporary.exists():
                temporary.unlink()

    def update(self, stage=None, **fields):
        if stage is not None:
            self.record["stage"] = stage
        self.record.update(fields)
        self.record["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
        self._persist()


def collector_source_preflight(repo, source_commit, profile):
    """Reuse the unchanged collector checks before spending analysis time."""
    with clean_collection_worktree(repo, source_commit) as worktree:
        with collection_module(worktree) as detached_collector, redirect_stdout(sys.stderr):
            identity = profile["source_identity"]
            require(identity["required_analysis_commit"] == source_commit,
                    "collector_preflight_required_commit_mismatch")
            detached_profile = detached_collector.load_validation_profile(worktree / PROFILE)
            detached_profile["source_identity"]["required_analysis_commit"] = source_commit
            require(detached_profile == profile, "collector_preflight_profile_mismatch")
            _, checks = detached_collector.collect_source_checks(
                worktree, required_analysis_commit=source_commit,
                allowed_committed_files=identity["allowed_committed_files"],
            )
            for check in checks:
                require(check["returncode"] == 0,
                        "collector_source_check_failed:{}:returncode={}:{}".format(
                            check["name"], check["returncode"], check["stderr"] +
                            ("\nstdout: " + check["stdout"] if check["stdout"] else ""),
                        ))
            ancestor, _, unexpected = detached_collector._committed_identity(
                checks, identity["allowed_committed_files"],
                identity["allowed_non_analysis_path_prefixes"],
            )
            require(ancestor, "collector_preflight_required_commit_not_present")
            require(not unexpected, "collector_preflight_unexpected_committed_files:" +
                    json.dumps(unexpected, sort_keys=True))
    return checks


def execute_gate(repo, outdir, output, source_commit):
    status = GateStatus(outdir, output, source_commit)
    try:
        require(not output.exists(), "output_zip_already_exists")
        require(output.parent.is_dir(), "output_directory_missing")
        provenance = preflight(repo, source_commit)
        status.update("profile")
        profile = resolved_profile(repo, source_commit)
        status.update("candidate_identity")
        for name, expected in CANDIDATES.items():
            require((outdir / name).is_file() and sha256(outdir / name) == expected, "candidate_identity_mismatch:" + name)
        before = {path.name: (path.stat().st_mtime_ns, path.stat().st_size, sha256(path))
                  for path in outdir.iterdir() if path.is_file() and path.name in {
                      declaration["basename_template"].format(kinematic=KINEMATIC, phi="Left", epsilon="lowe")
                      for scope in ("global", "settings") for declaration in profile["artifacts"][scope]}}
        log_path = outdir / (output.stem + ".log")
        require(not log_path.exists(), "run_log_already_exists")
        status.update("collector_source_preflight")
        checks = collector_source_preflight(repo, source_commit, profile)
        status.update(collector_source_preflight_completed=True, collector_source_checks=checks)
        started_ns = time.time_ns()
        status.update("analysis", analysis_started=True)
        require(run_analysis(repo, log_path) == 0, "analysis_failed")
        status.update("completion_markers", analysis_completed=True)
        text = log_path.read_text(encoding="utf-8", errors="replace")
        require("Left / lowe debug analysis completed." in text and
                "Full high-epsilon processing is intentionally skipped in -d debug mode." in text, "left_lowe_completion_missing")
        status.update("verify_artifacts")
        records, page_count = verify_artifacts(outdir, profile, started_ns, before)
        status.update(artifact_verification_completed=True)
        require(git(repo, "rev-parse", "HEAD") == source_commit, "head_changed_during_run")
        status.update("write_run_summary")
        summary = {"schema_version": "e8_4_fix5_left_lowe_run_summary/v1", "source_commit": source_commit,
                   "profile_template_sha256": sha256(repo / PROFILE), "pre_run_provenance": provenance,
                   "post_run_status_short": git(repo, "status", "--porcelain=v1", "--untracked-files=all").splitlines(),
                   "analysis_returncode": 0, "analysis_command": ["./run_Prod_Analysis.sh", "-d", "4p4", "2p74"],
                   "run_log_path": str(log_path), "run_log_sha256": sha256(log_path),
                   "started_ns": started_ns, "finished_at_utc": datetime.now(timezone.utc).isoformat(),
                   "page_count": page_count, "renderer_failures": [], "artifacts": records}
        summary_path = outdir / SUMMARY
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        records[SUMMARY] = {"sha256": sha256(summary_path), "bytes": summary_path.stat().st_size}
        status.update("collection")
        with tempfile.TemporaryDirectory(prefix="kaonlt-e8-4-fix5-profile-") as directory:
            effective_profile = Path(directory) / "profile.json"
            effective_profile.write_text(json.dumps(profile), encoding="utf-8")
            with clean_collection_worktree(repo, source_commit) as worktree:
                with collection_module(worktree) as detached_collector, redirect_stdout(sys.stderr):
                    result = detached_collector.collect_validation_bundle(outdir=outdir, kinematic=KINEMATIC,
                        output=output, profile_path=effective_profile, repo_root=worktree,
                        phi="Left", epsilon="lowe")
                require(result["returncode"] == 0, "collection_failed")
                status.update("verify_zip", collection_completed=True)
                verify_zip(output, source_commit, records)
                status.update("collection_cleanup", zip_verification_completed=True)
        status.update("complete", status="success")
        return output
    except Exception as exc:
        try:
            status.update(status="failed", failure_reason=str(exc))
        except OSError as persistence_error:
            print("Owner status persistence failed: {}".format(persistence_error), file=sys.stderr)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True, help="exact reviewed/pushed full SHA, verified against HEAD and local origin/test")
    args = parser.parse_args(argv)
    output = GLOBUS / ("KaonLT_E8_4_Fix5_Left_lowe_Q4p4W2p74_" + datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-%f") + ".zip")
    try:
        result = execute_gate(REPO, ARTIFACTS, output, args.source_commit)
    except Exception as exc:
        print("E.8.4.Fix.5 gate failed: {}".format(exc), file=sys.stderr)
        print("Owner gate status: {}".format(gate_status_path(ARTIFACTS, output).as_posix()), file=sys.stderr)
        return 1
    print(result.as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
