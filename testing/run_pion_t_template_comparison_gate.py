#!/usr/bin/env python3
"""Own one isolated Center/Left lowe template comparison, then verify and collect.

Farm execution requires independent diff, publication, synchronization and
readiness review. This module imports operational helpers only, never physics.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import redirect_stdout
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import time
import zipfile

try:
    from . import run_e8_2_left_lowe_scientific_audit_gate as debug
except ImportError:
    import run_e8_2_left_lowe_scientific_audit_gate as debug

isolation = debug.isolation
collector, require, sha256, git = debug.collector, debug.require, debug.sha256, debug.git
BASE_HEAD = "554ad84630ad77920ce17954d99f73df8ed2e6ac"
PROFILE = "testing/pion_t_template_comparison_gate_profile.json"
PROFILE_ID = "pion_t_template_comparison_center_left_lowe/v1"
ENABLE_ENV = "LT_ANALYSIS_PION_T_TEMPLATE_COMPARISON"
COMMAND = ["./run_Prod_Analysis.sh", "-d", "4p4", "2p74"]
SETTINGS = [dict(row) for row in isolation.CANONICAL_SETTINGS]
COMPONENTS = ("pi_n", "pi_delta", "pi_sidis")
SOURCE_FILES = (
    "src/cuts/pion_t_template_comparison.py", "src/utility/pion_component_shapes.py",
    "src/cuts/rand_sub.py", "src/cuts/pion_component_fits.py", "src/cuts/pion_t_bin_parents.py",
    "src/cuts/apply_cuts.py", "src/utility/background_config.py", "src/utility/canonical_binning.py",
    "src/cuts/hgcer_hole.py", "src/cuts/data_coordinates.py",
    "src/utility/root_histogram_ownership.py", "run_Prod_Analysis.sh")
ARTIFACTS = {"global": [
    {"key": setting.lower() + "_comparison_" + extension,
     "basename_template": f"{setting}_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.{extension}",
     "kind": "json" if extension == "json" else "file", "required": True}
    for setting in ("Center", "Left") for extension in ("json", "pdf")], "settings": []}
BOUNDARIES = {"diagnostic_only": True, "scientifically_accepted": False,
    "production_applied": False, "ordinary_production_outputs_may_be_written": True,
    "whole_output_root_immutable": False, "pdf_legibility_verified": False,
    "method_a_detached": True, "method_b_numerically_excluded": True}


def strict_json(raw):
    def invalid(value):
        raise ValueError("nonfinite_json:" + value)
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate_json_key:" + key)
            result[key] = value
        return result
    result = json.loads(raw, parse_constant=invalid, object_pairs_hook=pairs)
    def finite(value):
        if isinstance(value, float):
            require(math.isfinite(value), "nonfinite_json_number")
        elif isinstance(value, dict):
            for item in value.values():
                finite(item)
        elif isinstance(value, list):
            for item in value:
                finite(item)
    finite(result)
    return result


def resolved_profile(repo, source_commit):
    require(re.fullmatch(r"[0-9a-f]{40}", source_commit) is not None, "source_commit_requires_full_sha")
    profile = collector.load_validation_profile(Path(repo) / PROFILE)
    expected = {"schema_version": "pion_hgcer_validation_bundle_profile/v4",
        "validation_profile": PROFILE_ID, "collection_mode": "generic_artifacts",
        "settings": SETTINGS, "artifacts": ARTIFACTS,
        "source_identity": {"required_analysis_commit": BASE_HEAD, "allowed_committed_files": [],
                            "allowed_non_analysis_path_prefixes": ["docs/memory/"]}}
    require(profile == expected, "profile_contract_mismatch")
    profile = deepcopy(profile)
    profile["source_identity"]["required_analysis_commit"] = source_commit
    return profile


def names():
    return {row["basename_template"] for row in ARTIFACTS["global"]}


def preflight(repo, source_commit):
    require(re.fullmatch(r"[0-9a-f]{40}", source_commit) is not None, "source_commit_requires_full_sha")
    before = isolation.snapshot(repo)
    require(before["branch"] == "test", "wrong_branch")
    require(before["head"] == source_commit, "wrong_head")
    require(before["origin_test"] == source_commit, "source_not_observed_pushed")
    require(source_commit != BASE_HEAD, "owner_source_not_committed")
    git(repo, "merge-base", "--is-ancestor", BASE_HEAD, source_commit)
    for line in before["porcelain"].splitlines():
        path = line[3:]
        require(path in isolation.FARM_OUTPUTS or path.startswith("OUTPUT/") or
                (line[:2] == "??" and re.fullmatch(r"kaonlt_review(?:\([^/]+\))?\.diff", path)),
                "dirty_gate_source:" + line)
    # This successor may add only the three owned files and repository memory.
    changed = git(repo, "diff", "--name-only", BASE_HEAD, source_commit).splitlines()
    allowed = {PROFILE, "testing/run_pion_t_template_comparison_gate.py",
               "testing/test_run_pion_t_template_comparison_gate.py"}
    require(all(p in allowed or p.startswith("docs/memory/") for p in changed),
            "unexpected_committed_source")
    return before


def committed_identities(repo, source_commit):
    records = {}
    for path in SOURCE_FILES:
        result = subprocess.run(["git", "show", source_commit + ":" + path], cwd=repo,
                                capture_output=True)
        require(result.returncode == 0 and result.stdout, "committed_source_missing:" + path)
        records[path] = {"sha256": hashlib.sha256(result.stdout).hexdigest(), "bytes": len(result.stdout)}
    return records


def source_preflight(repo, source_commit, profile):
    with isolation.owned_worktree(repo, source_commit, "source") as source:
        require(resolved_profile(source, source_commit) == profile, "detached_profile_mismatch")
        with isolation.collection_module(source) as module, redirect_stdout(sys.stderr):
            _, checks = module.collect_source_checks(source,
                required_analysis_commit=source_commit, allowed_committed_files=[])
            require(all(c["returncode"] == 0 for c in checks), "collector_source_check_failed")
            ancestor, _, unexpected = module._committed_identity(checks, [], ["docs/memory/"])
            require(ancestor and not unexpected, "collector_source_identity_invalid")
        identities = committed_identities(source, source_commit)
    return checks, identities


def child_environment(env):
    result = dict(env)
    for key in ("LT_ANALYSIS_DEBUG_LEFT_LOW", "LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE",
                "LT_ANALYSIS_ALLOW_UNPAIRED_CANONICAL_BINNING"):
        result.pop(key, None)
    require(not any(k.startswith("LT_ANALYSIS_") and k != ENABLE_ENV and v
                    for k, v in result.items()), "untrusted_analysis_environment")
    require(result.get(ENABLE_ENV, "1") == "1", "conflicting_comparison_environment")
    result[ENABLE_ENV] = "1"
    return result


def verify_completion(log_path, returncode):
    require(returncode == 0, "analysis_failed")
    text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    require("All available low-epsilon settings completed in debug mode." in text and
            "Full high-epsilon processing is intentionally skipped in -d debug mode." in text,
            "all_low_completion_missing")
    require("High Epsilon Completed!" not in text, "unexpected_high_epsilon_execution")


def pdf_pages(path):
    # Existing Poppler backend; do not install packages or use regex page counts.
    result = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True)
    require(result.returncode == 0, "pdf_structure_unreadable:" + result.stderr)
    pages = re.findall(r"(?m)^Pages:\s+(\d+)\s*$", result.stdout)
    require(len(pages) == 1, "pdf_page_count_unavailable")
    return int(pages[0])


def fingerprint(hist):
    require(isinstance(hist, dict) and hist.get("status") == "available" and
            isinstance(hist.get("stored_sumw2"), list) and hist["stored_sumw2"] and
            re.fullmatch(r"[0-9a-f]{64}", hist.get("content_error_fingerprint", "")) is not None,
            "frozen_histogram_identity_invalid")
    return hist["content_error_fingerprint"]


def verify_payload(payload, setting, source_commit, identities, pdf_record):
    require(payload.get("schema_version") == "pion_t_template_comparison/v1" and
            payload.get("diagnostic_only") is True and payload.get("production_applied") is False and
            payload.get("scientifically_accepted") is False, "comparison_schema_or_authority_invalid")
    require(payload.get("source_HEAD") == source_commit and
            payload.get("source_identity_rejection_reasons") == [], "comparison_source_identity_invalid")
    require(payload.get("phi_setting") == setting and payload.get("setting") in (
        {"Q2": "4p4", "W": "2p74", "EPSSET": "low", "ParticleType": "kaon"},
        {"Q2": "4p4", "W": "2p74", "EPSSET": "lowe", "ParticleType": "kaon"}) and
        payload.get("enable_environment_variable") == ENABLE_ENV, "comparison_setting_invalid")
    source = payload.get("source_files") or {}
    require(set(source) == set(SOURCE_FILES), "comparison_source_inventory_invalid")
    for path, identity in source.items():
        require(identity.get("status") == "available" and
                all(identity.get(k) == v for k, v in identities[path].items()),
                "committed_source_hash_mismatch:" + path)
    pdf = payload.get("pdf") or {}
    require(pdf.get("status") == "available" and pdf.get("sha256") == pdf_record["sha256"] and
            pdf.get("bytes") == pdf_record["bytes"] and
            Path(pdf.get("path", "")).name == pdf_record["basename"], "comparison_pdf_identity_invalid")
    edges = payload.get("canonical_t_edges")
    require(isinstance(edges, list) and len(edges) == 4 and
            all(type(x) in (int, float) and math.isfinite(x) for x in edges) and
            all(a < b for a, b in zip(edges, edges[1:])), "canonical_t_edges_invalid")
    parents = payload.get("parents")
    require(isinstance(parents, list) and len(parents) == 3 and
            sorted(p.get("t_bin_index", -1) for p in parents) == [0, 1, 2], "comparison_t_inventory_invalid")
    for parent in parents:
        ti = parent["t_bin_index"]
        require(type(ti) is int and parent.get("t_edges") == edges[ti:ti+2], "parent_t_edges_invalid")
        require(parent.get("authoritative_parent_unchanged") is True and
                parent.get("production_applied") is False and parent.get("scientifically_accepted") is False and
                parent.get("alignment_mode") == "fixed_parent_alignment_template_shape_test" and
                re.fullmatch(r"[0-9a-f]{64}", parent.get("baseline_parent_fingerprint", "")) is not None and
                all(parent.get(k) for k in ("pion_parent_id", "canonical_interval_pair_id",
                                           "canonical_interval_pair_hash", "coordinate_fingerprint")),
                "parent_authority_or_identity_invalid")
        baseline, candidate = parent.get("baseline_A"), parent.get("candidate_B")
        require(isinstance(baseline, dict) and baseline, "baseline_A_missing")
        for side, key in (("pion", "H_pion_control_input"), ("kaon", "H_kaon_nosub_input")):
            frozen = fingerprint((parent.get("frozen_inputs") or {}).get(side))
            require(frozen == fingerprint(baseline.get(key)), "baseline_input_identity_mismatch")
            if candidate is not None:
                require(frozen == fingerprint(candidate.get(key)), "candidate_input_identity_mismatch")
        reasons = parent.get("rejection_reasons")
        require(isinstance(reasons, list) and all(isinstance(r, str) and r for r in reasons),
                "candidate_rejection_reasons_invalid")
        require(set(parent.get("templates") or {}) == set(COMPONENTS), "component_inventory_invalid")
        for component, evidence in parent["templates"].items():
            require(isinstance(evidence, dict) and isinstance(evidence.get("provenance"), dict) and
                    evidence["provenance"], "component_provenance_missing:" + component)
            provenance = evidence["provenance"]
            if provenance.get("status") == "available":
                require(re.fullmatch(r"[0-9a-f]{64}", provenance.get("sha256", "")) is not None and
                        provenance.get("path") and provenance.get("tree_name") and
                        provenance.get("weight_expression") == "getattr(evt, 'iter_weight', 1.0)" and
                        provenance.get("coordinate") == "missmass; -t; phipq*180/pi",
                        "component_provenance_invalid:" + component)
            else:
                require(provenance.get("status") == "unavailable" and provenance.get("reason") and
                        parent.get("status") == "unavailable" and reasons,
                        "component_provenance_invalid:" + component)
            if parent.get("status") == "available":
                require(evidence.get("status") == "genuinely_local" and
                        evidence.get("normalization_count") == 1 and
                        math.isclose(evidence.get("normalized_integral", 0), 1., abs_tol=1e-9) and
                        all(evidence.get(k) is True for k in ("canonical_selection_closure",
                            "weighted_child_to_parent_closure", "production_global_selection_closure")) and
                        evidence.get("t_edges") == edges[ti:ti+2], "component_normalization_invalid")
        if parent.get("status") == "available":
            require(isinstance(candidate, dict) and candidate and not reasons, "available_candidate_invalid")
            require(candidate.get("pion_component_alignment") == baseline.get("pion_component_alignment") and
                    candidate.get("resolved_subtraction_config") == baseline.get("resolved_subtraction_config"),
                    "candidate_alignment_or_configuration_changed")
        else:
            require(parent.get("status") == "unavailable" and candidate is None and reasons,
                    "unexplained_unavailable_candidate")
    expected = "available" if all(p["status"] == "available" for p in parents) else "partial_unavailable"
    require(payload.get("status") == expected, "comparison_status_invalid")
    return {"status": expected, "rejection_reasons_by_t": {
        str(p["t_bin_index"]): p["rejection_reasons"] for p in parents}, "canonical_t_edges": edges}


def verify_artifacts(outdir, source_commit, identities, started_ns, before):
    records, outcomes = {}, {}
    for name in sorted(names()):
        path = outdir / name
        require(path.is_file() and not path.is_symlink() and path.stat().st_size > 0,
                "comparison_artifact_missing_or_empty:" + name)
        info = path.stat()
        require(info.st_mtime_ns >= started_ns and before[name] is None,
                "comparison_artifact_stale:" + name)
        records[name] = {"basename": name, "bytes": info.st_size, "sha256": sha256(path),
                         "mtime_ns": info.st_mtime_ns}
        if path.suffix == ".pdf":
            require(path.read_bytes().startswith(b"%PDF-"), "comparison_pdf_signature_invalid")
            require(pdf_pages(path) == 6, "comparison_pdf_page_count_invalid")
    for setting in ("Center", "Left"):
        stem = f"{setting}_kaon_Q4p4W2p74_lowe_pion-t-template-comparison"
        payload = strict_json((outdir / (stem + ".json")).read_text(encoding="utf-8"))
        outcomes[setting] = verify_payload(payload, setting, source_commit, identities, records[stem + ".pdf"])
    require(outcomes["Center"]["canonical_t_edges"] == outcomes["Left"]["canonical_t_edges"],
            "mixed_canonical_t_edges")
    return records, outcomes


def verify_zip(path, source_commit, records):
    expected = {"manifest.json", "source_state.txt", "source_checks.txt", "global/"}
    expected.update(s["phi"] + "_" + s["epsilon"] + "/" for s in SETTINGS)
    expected.update("global/" + name for name in names())
    with zipfile.ZipFile(path) as archive:
        require(archive.testzip() is None, "zip_crc_failed")
        require(Counter(archive.namelist()) == Counter(expected), "zip_inventory_invalid")
        manifest = strict_json(archive.read("manifest.json"))
        require(manifest.get("complete") is True and manifest.get("errors") == [] and
                manifest.get("git_head") == manifest.get("required_analysis_commit") == source_commit and
                manifest.get("validation_profile") == PROFILE_ID and
                manifest.get("requested_kinematic") == "Q4p4W2p74" and
                manifest.get("requested_settings") == SETTINGS, "zip_manifest_invalid")
        require(len(manifest.get("settings", [])) == 5 and
                all(s.get("artifacts") == {} for s in manifest["settings"]), "zip_setting_artifacts_invalid")
        require(set(manifest.get("global_artifacts", {})) == {r["key"] for r in ARTIFACTS["global"]},
                "zip_global_inventory_invalid")
        for declaration in ARTIFACTS["global"]:
            name = declaration["basename_template"]
            raw = archive.read("global/" + name)
            record = manifest["global_artifacts"][declaration["key"]]
            require(len(raw) == records[name]["bytes"] and
                    hashlib.sha256(raw).hexdigest() == records[name]["sha256"] and
                    record.get("sha256") == records[name]["sha256"] and
                    record.get("byte_size") == records[name]["bytes"], "zip_artifact_identity_mismatch:" + name)


def destinations(outdir, output):
    return {"log": outdir / (output.stem + ".log"),
            "summary": outdir / (output.stem + "-run-summary.json"),
            "receipt": isolation.accepted.gate_status_path(outdir, output),
            "zip": output,
            "delivered_log": output.with_name(output.stem + ".log"),
            "delivered_summary": output.with_name(output.stem + "-run-summary.json"),
            "delivered_receipt": output.with_name(output.stem + "-gate-status.json")}


class GateStatus(isolation.accepted.GateStatus):
    def __init__(self, outdir, output, source_commit):
        super().__init__(outdir, output, source_commit)
        self.record.pop("phi", None)
        self.record.pop("epsilon", None)
        self.update(schema_version="pion_t_template_comparison_gate_status/v1",
                    requested_settings=[{"phi": s, "epsilon": "lowe"} for s in ("Center", "Left")],
                    boundaries=BOUNDARIES)


def execute_gate(repo, outdir, output, source_commit):
    repo, outdir, output = (Path(p).absolute() for p in (repo, outdir, output))
    # Use the established farm-owner boundary; no new machine-specific values.
    require(repo.resolve() == isolation.accepted.REPO.resolve() and
            outdir.resolve() == isolation.accepted.ARTIFACTS.resolve() and
            output.parent.resolve() == isolation.accepted.GLOBUS.resolve() and
            output.suffix == ".zip" and outdir.is_dir() and output.parent.is_dir(), "owner_path_boundary_invalid")
    require(not output.is_symlink(), "output_symlink_invalid")
    paths = destinations(outdir, output)
    require(len(set(map(str, paths.values()))) == len(paths), "companion_paths_overlap")
    for path in paths.values():
        require(not os.path.lexists(path), "output_collision:" + str(path))
    status = GateStatus(outdir, output, source_commit)
    before = overlay = installed = None
    succeeded = False
    try:
        before = preflight(repo, source_commit)
        profile = resolved_profile(repo, source_commit)
        child_environment(os.environ)  # Reject injected opt-ins before any path work.
        prior = {name: None if not os.path.lexists(outdir / name) else
                 {"path": str(outdir / name)} for name in names()}
        status.update(preexisting_comparison_targets=prior, destinations={k: str(v) for k,v in paths.items()})
        require(all(value is None for value in prior.values()), "comparison_artifact_preexists")
        checks, identities = source_preflight(repo, source_commit, profile)
        status.update("runtime_isolation", collector_source_preflight_completed=True,
                      collector_source_checks=checks, committed_source_identities=identities)
        # Snapshot the original package before the frozen overlay helper: its
        # own return record is unavailable if copy/import/path verification fails.
        identity = isolation.import_identity(repo)
        package = Path(identity["package_file"]).resolve().parent
        protected = [Path(identity["package_file"]), Path(identity["pathing_file"]),
                     *sorted((package / "PATH_TO_DIR").glob("*.path"))]
        installed = {"original_source_sha256": {str(p): sha256(p) for p in protected}}
        with isolation.owned_worktree(repo, source_commit, "analysis") as worktree:
            env, overlay = isolation.prepare_runtime_overlay(repo, worktree)
            env = child_environment(env)
            probes = isolation.probe_paths(worktree, env, overlay["baseline_paths"])
            external = isolation.external_symlink_preflight(worktree, probes[0]["paths"])
            if external["background_simc_path"]:
                require(Path(external["background_simc_path"]).is_dir(), "background_simc_directory_missing")
            link = debug.prepare_debug_output_link(worktree, probes[0]["paths"], outdir.resolve())
            status.update("analysis", analysis_started=True, analysis_worktree=str(worktree),
                          isolation_probes=probes, debug_output_link=link, external_symlink_preflight=external)
            started = time.time_ns()
            returncode = debug.run_analysis(worktree, paths["log"], env)
            status.update("completion_markers", analysis_returncode=returncode,
                          log_identity={"sha256": sha256(paths["log"]), "bytes": paths["log"].stat().st_size})
            verify_completion(paths["log"], returncode)
            status.update("verify_artifacts", analysis_completed=True)
            records, outcomes = verify_artifacts(outdir, source_commit, identities, started, prior)
            status.update(artifact_verification_completed=True, artifacts=records, diagnostic_outcomes=outcomes)
        status.update(worktree_cleanup_completed=True)
        status.update(ordinary_checkout_preservation=isolation.verify_preservation(repo, before),
                      installed_ltsep_preservation=isolation.verify_ltsep_preservation(installed))
        status.update("collection")
        with tempfile.TemporaryDirectory(prefix="kaonlt-pion-t-profile-") as directory:
            effective = Path(directory) / "profile.json"
            effective.write_text(json.dumps(profile), encoding="utf-8")
            with isolation.owned_worktree(repo, source_commit, "source") as source:
                require(resolved_profile(source, source_commit) == profile, "detached_profile_mismatch")
                with isolation.collection_module(source) as module, redirect_stdout(sys.stderr):
                    result = module.collect_validation_bundle(outdir=outdir, kinematic="Q4p4W2p74",
                        output=output, profile_path=effective, repo_root=source)
                require(result["returncode"] == 0, "collection_failed")
                status.update("verify_zip", collection_completed=True)
                verify_zip(output, source_commit, records)
                status.update(zip_verification_completed=True, zip_identity={"sha256": sha256(output),
                              "bytes": output.stat().st_size})
        status.update(final_ordinary_checkout_preservation=isolation.verify_preservation(repo, before),
                      final_installed_ltsep_preservation=isolation.verify_ltsep_preservation(installed))
        summary = {"schema_version": "pion_t_template_comparison_run_summary/v1",
            "source_commit": source_commit, "analysis_command": COMMAND, "analysis_returncode": returncode,
            "diagnostic_outcomes": outcomes, "artifacts": records, "boundaries": BOUNDARIES,
            "log_identity": status.record["log_identity"], "zip_identity": status.record["zip_identity"]}
        with paths["summary"].open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        isolation.deliver_evidence(status, paths["log"], paths["summary"], output)
        succeeded = True
        return output
    except Exception as exc:
        status.update(status="failed", failure_reason=str(exc))
        raise
    finally:
        # Every failure after the initial snapshot gets preservation checks,
        # including source, overlay, child, cleanup, artifact and ZIP failures.
        try:
            if before is not None and not succeeded:
                status.update(final_ordinary_checkout_preservation=isolation.verify_preservation(repo, before))
        except Exception as exc:
            status.update(status="failed", preservation_failure_reason=str(exc))
            raise
        finally:
            if installed is not None and not succeeded:
                try:
                    status.update(final_installed_ltsep_preservation=isolation.verify_ltsep_preservation(installed))
                except Exception as exc:
                    status.update(status="failed", ltsep_preservation_failure_reason=str(exc))
                    raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        result = execute_gate(args.repo, args.outdir, args.output, args.source_commit)
    except Exception as exc:
        print("Pion template gate failed: " + str(exc), file=sys.stderr)
        print("Owner gate status: " + str(isolation.accepted.gate_status_path(args.outdir, args.output)), file=sys.stderr)
        return 1
    print(result.as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
