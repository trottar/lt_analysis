#!/usr/bin/env python3
"""Own the isolated E.8.2 baseline Left/lowe audit; no local farm authorization.

The v4 profile authorizes canonical five; this owner requests only Left/lowe.
Review, push, synchronization and farm-readiness audit precede farm execution.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import redirect_stdout
import csv
from datetime import datetime, timezone
import hashlib
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
    from . import run_e8_4_fix5_canonical_five_plot_gate as isolation
except ImportError:
    import run_e8_4_fix5_canonical_five_plot_gate as isolation

# Reuse operational isolation only. No candidate/lineage/page gates are called.
collector = isolation.collector
require = isolation.require
sha256 = isolation.sha256
git = isolation.git
owned_worktree = isolation.owned_worktree
collection_module = isolation.collection_module
prepare_runtime_overlay = isolation.prepare_runtime_overlay
probe_paths = isolation.probe_paths
external_symlink_preflight = isolation.external_symlink_preflight
verify_preservation = isolation.verify_preservation
verify_ltsep_preservation = isolation.verify_ltsep_preservation

BASE_HEAD = "7b8eb3cb231d289d18de9c168960ddddcaa39254"
KINEMATIC = "Q4p4W2p74"
PROFILE = "testing/pion_hgcer_validation_bundle_profile_e8_2_left_lowe.json"
PROFILE_ID = "phase_e8_2_left_lowe_baseline_scientific_audit/v1"
SUMMARY = KINEMATIC + "_e8_2_left_lowe_scientific-audit-run-summary.json"
PDF = "Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf"
COMMAND = ["./run_Prod_Analysis.sh", "-d", "4p4", "2p74"]
DECLARED_SETTINGS = [dict(row) for row in isolation.CANONICAL_SETTINGS]
REQUESTED_SETTINGS = [{"phi": "Left", "epsilon": "lowe"}]
STAGES = {"random": "random_subtraction", "dummy": "dummy_subtraction",
          "proton": "slow_proton_cleaning", "pion": "baseline_pion_subtraction",
          "final_mm": "final_baseline_mm", "stage_yields": "baseline_stage_yields"}
PAGE_IDS = {f"full_background.e8_2.{stage}.t{n}": (n, semantic)
            for stage, semantic in STAGES.items() for n in (1, 2, 3)}
SETTING_KEYS = {"procedure_pdf", "page_manifest", "full_analysis",
                "correction_ledger_json", "correction_ledger_csv"}
ARTIFACTS = {"global": [{"key": "run_summary",
    "basename_template": "{kinematic}_e8_2_left_lowe_scientific-audit-run-summary.json",
    "kind": "json", "required": True}], "settings": [
    {"key": key, "basename_template": basename, "kind": kind, "required": True}
    for key, basename, kind in (
        ("procedure_pdf", "{phi}_kaon_rand_sub_{kinematic}_{epsilon}_full-background-subtraction.pdf", "file"),
        ("page_manifest", "{phi}_kaon_rand_sub_{kinematic}_{epsilon}_full-background-subtraction-manifest.json", "json"),
        ("full_analysis", "kaon_FullAnalysis_{kinematic}_{epsilon}.json", "json"),
        ("correction_ledger_json", "kaon_FullAnalysis_{kinematic}_{epsilon}_correction_ledger_no_empirical_residual.json", "json"),
        ("correction_ledger_csv", "kaon_FullAnalysis_{kinematic}_{epsilon}_correction_ledger_no_empirical_residual.csv", "file"))]}
BOUNDARIES = {"scope": "E.8.2 baseline audit only", "method_a_required": False,
              "method_b_numerically_excluded": True, "production_promotion": False,
              "e8_3_acceptance": False, "e8_4_acceptance": False,
              "canonical_five_acceptance": False, "absolute_simc_amplitude_claim": False}


def preflight(repo, source_commit):
    require(re.fullmatch(r"[0-9a-f]{40}", source_commit) is not None,
            "source_commit_requires_full_sha")
    state = isolation.snapshot(repo)
    require(state["branch"] == "test", "wrong_branch")
    require(state["head"] == source_commit, "wrong_head")
    require(state["origin_test"] == source_commit, "source_not_observed_pushed")
    require(source_commit != BASE_HEAD, "owner_source_not_committed")
    git(repo, "merge-base", "--is-ancestor", BASE_HEAD, source_commit)
    for line in state["porcelain"].splitlines():
        path = line[3:]
        require(path in isolation.FARM_OUTPUTS or path.startswith("OUTPUT/") or
                (line[:2] == "??" and re.fullmatch(r"kaonlt_review(?:\([^/]+\))?\.diff", path)),
                "dirty_gate_source:" + line)
    git(repo, "diff", "--check", "HEAD", "--", *isolation.GATE_PATHS,
        *(":(exclude)" + name for name in sorted(isolation.FARM_OUTPUTS)))
    return state


def resolved_profile(repo, source_commit):
    require(re.fullmatch(r"[0-9a-f]{40}", source_commit) is not None,
            "source_commit_requires_full_sha")
    profile = collector.load_validation_profile(Path(repo) / PROFILE)
    require(profile["settings"] == DECLARED_SETTINGS and
            collector.resolve_settings("Left", "lowe", profile) == (("Left", "lowe"),),
            "profile_inventory_invalid")
    require(profile["validation_profile"] == PROFILE_ID and
            profile["collection_mode"] == "generic_artifacts", "profile_mode_invalid")
    require(profile["source_identity"] == {"required_analysis_commit": BASE_HEAD,
            "allowed_committed_files": [], "allowed_non_analysis_path_prefixes": ["docs/memory/"]},
            "profile_source_identity_invalid")
    require(profile["artifacts"] == ARTIFACTS, "profile_artifact_inventory_invalid")
    profile["source_identity"]["required_analysis_commit"] = source_commit
    return profile


def artifact_names(profile):
    return {entry["basename_template"].format(kinematic=KINEMATIC, phi="Left", epsilon="lowe")
            for scope in ("global", "settings") for entry in profile["artifacts"][scope]}


def collector_source_preflight(repo, source_commit, profile):
    with owned_worktree(repo, source_commit, "source") as source:
        require(resolved_profile(source, source_commit) == profile, "collector_profile_mismatch")
        with collection_module(source) as module, redirect_stdout(sys.stderr):
            _, checks = module.collect_source_checks(source,
                required_analysis_commit=source_commit, allowed_committed_files=[])
            for check in checks:
                require(check["returncode"] == 0,
                        "collector_source_check_failed:" + check["name"] + ":" + check["stderr"])
            ancestor, _, unexpected = module._committed_identity(checks, [], ["docs/memory/"])
            require(ancestor and not unexpected, "collector_source_identity_invalid")
    return checks


def prepare_debug_output_link(worktree, paths, outdir):
    """Link OUTPUT only in the caller's owned disposable analysis worktree.

    Use the existing probe's authority before the debug launcher's early mkdir.
    Never replace an existing path or discover alternate runtime paths.
    """
    worktree, outdir = Path(worktree), Path(outdir)
    require(worktree.is_absolute() and worktree.is_dir() and not worktree.is_symlink(),
            "debug_output_worktree_invalid")
    require(worktree.resolve() == worktree and
            Path(paths.get("LTANAPATH", "")).resolve() == worktree,
            "debug_output_ltanapath_mismatch")
    volatile = Path(paths.get("VOLATILEPATH", ""))
    require(volatile.is_absolute(), "debug_output_volatilepath_not_absolute")
    anatype = paths.get("ANATYPE")
    require(isinstance(anatype, str) and re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", anatype),
            "debug_output_anatype_invalid")
    child = Path(paths.get("OUTPATH", ""))
    require(child.is_absolute() and child == worktree / "OUTPUT/Analysis" / (anatype + "LT"),
            "debug_output_lexical_outpath_mismatch")
    target = volatile / "OUTPUT"
    require(target.is_dir(), "debug_output_external_output_missing")
    require(outdir.is_absolute() and
            (target / "Analysis" / (anatype + "LT")).resolve() == outdir,
            "debug_output_artifact_root_mismatch")
    link = worktree / "OUTPUT"
    require(not os.path.lexists(link), "debug_output_path_already_exists")
    link.symlink_to(target, target_is_directory=True)
    require(link.is_symlink() and link.is_dir(), "debug_output_link_invalid")
    require(link.resolve() == target.resolve() and os.readlink(link) == str(target),
            "debug_output_link_target_mismatch")
    require(child.resolve() == outdir, "debug_output_child_resolution_mismatch")
    return {"link_path": str(link), "literal_target": str(target),
            "resolved_target": str(link.resolve()), "resolved_child_outpath": str(child.resolve()),
            "passed": True}


def run_analysis(worktree, log_path, env):
    with log_path.open("xb") as log:
        process = subprocess.Popen(COMMAND, cwd=worktree, env=env,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in iter(process.stdout.readline, b""):
            log.write(line)
            log.flush()
            sys.stderr.write(line.decode("utf-8", errors="replace"))
            sys.stderr.flush()
        return process.wait()


def verify_completion(log_path, returncode):
    require(returncode == 0, "analysis_failed")
    text = log_path.read_text(encoding="utf-8", errors="replace")
    require("Left / lowe debug analysis completed." in text and
            "Full high-epsilon processing is intentionally skipped in -d debug mode." in text,
            "left_lowe_completion_missing")
    require("High Epsilon Completed!" not in text, "unexpected_high_epsilon_execution")


def verify_pages(payload):
    require(isinstance(payload, dict) and
            payload.get("schema_version") == "full_background_subtraction_page_manifest/v1",
            "page_manifest_schema_invalid")
    setting = payload.get("setting")
    expected = {"kinematic_token": KINEMATIC, "epsilon_setting": "low",
                "epsilon_filename_token": "lowe", "phi_setting": "Left", "particle_type": "kaon"}
    require(isinstance(setting, dict) and all(setting.get(k) == v for k, v in expected.items()),
            "page_manifest_setting_invalid")
    require(payload.get("pdf_basename") == PDF, "page_manifest_pdf_invalid")
    require(payload.get("renderer_failures") == [], "renderer_failures_nonempty")
    pages = payload.get("pages")
    require(isinstance(pages, list) and all(isinstance(p, dict) for p in pages), "page_inventory_invalid")
    counts = Counter(p.get("page_id") for p in pages)
    require(all(counts[p] == 1 for p in PAGE_IDS), "e8_2_page_missing_or_duplicate")
    positions = []
    for position, page in enumerate(pages):
        if page.get("page_id") not in PAGE_IDS:
            continue
        n, semantic = PAGE_IDS[page["page_id"]]
        require(page.get("scope") == f"t{n}" and type(page.get("t_index")) is int and
                page["t_index"] == n - 1, "e8_2_parent_identity_invalid")
        require(page.get("semantic_stage") == semantic, "e8_2_semantic_stage_invalid")
        children = page.get("represented_phi_inventory")
        require(isinstance(children, list) and all(isinstance(c, dict) for c in children) and
                [c.get("phi_index") for c in children] == list(range(9)) and
                all(type(c.get("phi_index")) is int for c in children), "e8_2_phi_inventory_invalid")
        require(page.get("invalid_unavailable_children") == [], "e8_2_invalid_children")
        require(page.get("authoritative") is False and page.get("presentation_only") is True,
                "e8_2_presentation_flags_invalid")
        positions.append(position)
    handoff = [i for i, p in enumerate(pages) if p.get("page_id") == "full_background.e8.handoff"]
    require(len(handoff) == 1 and handoff[0] > max(positions), "e8_handoff_missing_or_misordered")
    return len(pages)


def verify_artifacts(outdir, profile, started_ns, before):
    records = {}
    for name in sorted(artifact_names(profile) - {SUMMARY}):
        path = outdir / name
        require(path.is_file() and path.stat().st_size > 0, "artifact_missing:" + name)
        digest = sha256(path)
        require(path.stat().st_mtime_ns >= started_ns and
                (path.stat().st_mtime_ns, path.stat().st_size, digest) != before.get(name),
                "artifact_stale:" + name)
        if path.suffix == ".json":
            collector._strict_json_payload(path)
        records[name] = {"sha256": digest, "bytes": path.stat().st_size}
    require((outdir / PDF).read_bytes().startswith(b"%PDF-"), "procedure_pdf_invalid")
    count = verify_pages(collector._strict_json_payload(outdir / PDF.replace(".pdf", "-manifest.json")))
    full = collector._strict_json_payload(outdir / f"kaon_FullAnalysis_{KINEMATIC}_lowe.json")
    expected = {"ParticleType": "kaon", "EPSSET": "low", "Q2": "4p4", "W": "2p74",
                "OutFilename": f"FullAnalysis_{KINEMATIC}_lowe"}
    require(all(full.get("inpDict", {}).get(k) == v for k, v in expected.items()),
            "full_analysis_identity_invalid")
    require(isinstance(full.get("histlist"), list) and
            [row.get("phi_setting") for row in full["histlist"]] == ["Left"],
            "full_analysis_setting_inventory_invalid")
    ledger = collector._strict_json_payload(outdir /
        f"kaon_FullAnalysis_{KINEMATIC}_lowe_correction_ledger_no_empirical_residual.json")
    require(all(ledger.get(k) == v for k, v in {"active_profile": "no_empirical_residual",
            "particle_type": "kaon", "epsset": "low", "q2": "4p4", "w": "2p74"}.items()) and
            [row.get("phi_setting") for row in ledger.get("settings", [])] == ["Left"],
            "ledger_identity_invalid")
    with (outdir / f"kaon_FullAnalysis_{KINEMATIC}_lowe_correction_ledger_no_empirical_residual.csv").open(
            encoding="utf-8", newline="") as handle:
        totals = [r.get("phi_setting") for r in csv.DictReader(handle) if r.get("row_kind") == "setting_total"]
    require(totals == ["Left"], "ledger_csv_setting_inventory_invalid")
    return records, count


def verify_zip(path, source_commit, records, profile):
    require(profile["settings"] == DECLARED_SETTINGS, "bundle_profile_declaration_invalid")
    with zipfile.ZipFile(path) as archive:
        require(archive.testzip() is None, "zip_integrity_failed")
        manifest = json.loads(archive.read("manifest.json"))
        require(manifest.get("complete") is True and manifest.get("errors") == [], "bundle_incomplete")
        require(manifest.get("git_head") == source_commit and
                manifest.get("required_analysis_commit") == source_commit, "bundle_source_identity_mismatch")
        require(manifest.get("validation_profile") == PROFILE_ID and
                manifest.get("requested_kinematic") == KINEMATIC, "bundle_profile_mismatch")
        require(manifest.get("requested_settings") == REQUESTED_SETTINGS, "bundle_requested_inventory_invalid")
        settings = manifest.get("settings", [])
        require(len(settings) == 1 and
                all(settings[0].get(k) == v for k, v in {**REQUESTED_SETTINGS[0], "kinematic": KINEMATIC}.items()),
                "bundle_setting_inventory_invalid")
        require(set(manifest.get("global_artifacts", {})) == {"run_summary"} and
                set(settings[0].get("artifacts", {})) == SETTING_KEYS, "bundle_artifact_inventory_invalid")
        expected_paths = {"manifest.json", "source_state.txt", "source_checks.txt", "global/", "Left_lowe/"}
        entries = list(manifest["global_artifacts"].values()) + list(settings[0]["artifacts"].values())
        seen = set()
        for entry in entries:
            name = Path(entry["archive_path"]).name
            prefix = "global/" if name == SUMMARY else "Left_lowe/"
            require(name in records and entry["archive_path"] == prefix + name and name not in seen,
                    "bundle_artifact_inventory_invalid")
            raw = archive.read(entry["archive_path"])
            require(hashlib.sha256(raw).hexdigest() == entry["sha256"] == records[name]["sha256"] and
                    len(raw) == entry["byte_size"] == records[name]["bytes"], "bundle_artifact_hash_mismatch:" + name)
            seen.add(name)
            expected_paths.add(entry["archive_path"])
        require(seen == set(records) and set(archive.namelist()) == expected_paths and
                len(archive.namelist()) == len(expected_paths), "bundle_archive_inventory_invalid")


class GateStatus(isolation.accepted.GateStatus):
    def __init__(self, outdir, output, source_commit):
        self.path = isolation.accepted.gate_status_path(outdir, output)
        now = datetime.now(timezone.utc).isoformat()
        self.record = {"schema_version": "e8_2_left_lowe_scientific_audit_gate_status/v1",
            "source_commit": source_commit, "kinematic": KINEMATIC, "phi": "Left", "epsilon": "lowe",
            "status": "running", "stage": "preflight", "failure_reason": None,
            "started_at_utc": now, "updated_at_utc": now, "expected_zip_path": str(output),
            "analysis_started": False, "analysis_completed": False,
            "artifact_verification_completed": False, "collector_source_preflight_completed": False,
            "collection_completed": False, "zip_verification_completed": False,
            "ordinary_checkout_preservation": {"passed": False},
            "installed_ltsep_preservation": {"passed": False}}
        self._persist(initial=True)


def execute_gate(repo, outdir, output, source_commit):
    repo, outdir, output = (Path(p).resolve() for p in (repo, outdir, output))
    status = GateStatus(outdir, output, source_commit)
    provenance = overlay = None
    try:
        require(not os.path.lexists(output), "output_zip_already_exists")
        require(output.parent.is_dir(), "output_directory_missing")
        provenance = preflight(repo, source_commit)
        profile = resolved_profile(repo, source_commit)
        before = {name: (p.stat().st_mtime_ns, p.stat().st_size, sha256(p))
                  for name in artifact_names(profile) if (p := outdir / name).is_file()}
        log_path = outdir / (output.stem + ".log")
        require(not os.path.lexists(log_path), "run_log_already_exists")
        status.update("collector_source_preflight")
        checks = collector_source_preflight(repo, source_commit, profile)
        status.update(collector_source_preflight_completed=True, collector_source_checks=checks)
        try:
            with owned_worktree(repo, source_commit, "analysis") as worktree:
                status.update("path_isolation")
                env, overlay = prepare_runtime_overlay(repo, worktree)
                status.update(ltsep_runtime_overlay=overlay)
                probes = probe_paths(worktree, env, overlay["baseline_paths"])
                require(outdir == (Path(probes[0]["paths"]["VOLATILEPATH"]) /
                                  "OUTPUT/Analysis/KaonLT").resolve(), "analysis_artifact_root_mismatch")
                output_link = prepare_debug_output_link(worktree, probes[0]["paths"], outdir)
                status.update(debug_output_link=output_link)
                external = external_symlink_preflight(worktree, probes[0]["paths"])
                identity = {"path": str(worktree), "head": source_commit, "detached": True, "owner_created": True}
                started_ns = time.time_ns()
                status.update("analysis", analysis_started=True, analysis_worktree=identity,
                              isolation_probes=probes, external_symlink_preflight=external)
                returncode = run_analysis(worktree, log_path, env)
                status.update("completion_markers", analysis_returncode=returncode)
                verify_completion(log_path, returncode)
                status.update("verify_artifacts", analysis_completed=True)
                records, count = verify_artifacts(outdir, profile, started_ns, before)
                status.update(artifact_verification_completed=True)
        finally:
            # Even a failed child/path/artifact gate must recheck the ordinary
            # checkout and installed ltsep after bounded worktree cleanup.
            try:
                preservation = verify_preservation(repo, provenance)
                status.update(ordinary_checkout_preservation=preservation)
            finally:
                if overlay is not None:
                    ltsep = verify_ltsep_preservation(overlay)
                    status.update(installed_ltsep_preservation=ltsep)
        status.update(worktree_cleanup_completed=True)
        summary = {"schema_version": "e8_2_left_lowe_scientific_audit_run_summary/v1",
            "source_commit": source_commit, "analysis_command": COMMAND, "analysis_returncode": returncode,
            "log_sha256": sha256(log_path), "page_count": count, "e8_2_required_pages": 18,
            "e8_2_page_verification_passed": True, "active_profile": "no_empirical_residual",
            "profile_declared_settings": DECLARED_SETTINGS, "requested_settings": REQUESTED_SETTINGS,
            "analysis_worktree": identity, "isolation_probes": probes, "external_symlink_preflight": external,
            "debug_output_link": output_link,
            "ordinary_checkout_preservation": preservation, "installed_ltsep_preservation": ltsep,
            "worktree_cleanup_completed": True, "collector_source_checks": checks,
            "artifacts": records, "boundaries": BOUNDARIES}
        summary_path = outdir / SUMMARY
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        records[SUMMARY] = {"sha256": sha256(summary_path), "bytes": summary_path.stat().st_size}
        status.update("collection")
        try:
            with tempfile.TemporaryDirectory(prefix="kaonlt-e8-2-profile-") as directory:
                effective = Path(directory) / "profile.json"
                effective.write_text(json.dumps(profile), encoding="utf-8")
                with owned_worktree(repo, source_commit, "source") as source:
                    with collection_module(source) as module, redirect_stdout(sys.stderr):
                        result = module.collect_validation_bundle(outdir=outdir, kinematic=KINEMATIC,
                            output=output, profile_path=effective, repo_root=source, phi="Left", epsilon="lowe")
                    require(result["returncode"] == 0, "collection_failed")
                    status.update("verify_zip", collection_completed=True)
                    verify_zip(output, source_commit, records, profile)
                    status.update(zip_verification_completed=True)
        finally:
            try:
                status.update(final_ordinary_checkout_preservation=verify_preservation(repo, provenance))
            finally:
                status.update(final_installed_ltsep_preservation=verify_ltsep_preservation(overlay))
        isolation.deliver_evidence(status, log_path, summary_path, output)
        return output
    except Exception as exc:
        status.update(status="failed", failure_reason=str(exc))
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--outdir", type=Path, required=True, help="configured analysis artifact root")
    parser.add_argument("--output", type=Path, required=True, help="fresh ZIP at configured transfer root")
    args = parser.parse_args(argv)
    try:
        result = execute_gate(args.repo, args.outdir, args.output, args.source_commit)
    except Exception as exc:
        print("E.8.2 gate failed: " + str(exc), file=sys.stderr)
        print("Owner gate status: " + str(isolation.accepted.gate_status_path(args.outdir, args.output)), file=sys.stderr)
        return 1
    print(result.as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
