#!/usr/bin/env python3
"""Independent baseline-five evidence owner; execution requires reviewed pushed source.

No candidate staging, lineage equivalence, production mutation, or SIMC rescaling.
Numerical integrations below are read-only checks of existing producer outputs.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import redirect_stdout
from copy import deepcopy
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time
import zipfile

import numpy as np

try:
    from . import run_e8_4_fix5_canonical_five_plot_gate as isolation
    from . import run_e8_2_left_lowe_scientific_audit_gate as baseline
except ImportError:
    import run_e8_4_fix5_canonical_five_plot_gate as isolation
    import run_e8_2_left_lowe_scientific_audit_gate as baseline

# Only operational primitives are reused, never candidate/staging/equivalence gates.
require = isolation.require
git = isolation.git
sha256 = isolation.sha256
owned_worktree = isolation.owned_worktree
collection_module = isolation.collection_module
prepare_runtime_overlay = isolation.prepare_runtime_overlay
probe_paths = isolation.probe_paths
external_symlink_preflight = isolation.external_symlink_preflight
verify_ltsep_preservation = isolation.verify_ltsep_preservation
run_analysis = isolation.run_analysis
verify_completion = isolation.verify_completion
collector = isolation.collector

BASE_HEAD = "af45126899c0d89707b29b95f7138b3388d432cd"
KINEMATIC = "Q4p4W2p74"
COMMAND = ["./run_Prod_Analysis.sh", "4p4", "2p74"]
SETTINGS = [dict(row) for row in isolation.CANONICAL_SETTINGS]
PROFILE = "testing/pion_hgcer_validation_bundle_profile_e8_2_canonical_five_yields.json"
PROFILE_ID = "phase_e8_2_canonical_five_baseline_yields/v1"
ATTEMPT_MARKER = "OWNER_ATTEMPT"
SUMMARY = KINEMATIC + "_e8_2_baseline-five_" + ATTEMPT_MARKER + "_yield-receipt.json"
SOURCE_IDENTITY = {"required_analysis_commit": BASE_HEAD, "allowed_committed_files": [],
                   "allowed_non_analysis_path_prefixes": ["docs/memory/"]}
BOUNDARIES = {"method_a_required": False, "method_a_acceptance": False,
              "method_b_numerically_excluded": True, "production_promotion": False,
              "runtime_acceptance": False, "absolute_simc_amplitude_claim": False,
              "absolute_simc_blocker": "SIMC_normfac_luminosity_and_charge_units_not_source_proven"}


def declaration(key, name, kind):
    return {"key": key, "basename_template": name, "kind": kind, "required": True}


ARTIFACTS = {"global": [declaration("run_summary", SUMMARY, "json")], "settings": [
    declaration("procedure_pdf", "{phi}_kaon_rand_sub_{kinematic}_{epsilon}_full-background-subtraction.pdf", "file"),
    declaration("page_manifest", "{phi}_kaon_rand_sub_{kinematic}_{epsilon}_full-background-subtraction-manifest.json", "json"),
    *[declaration("yield_" + kind, "{kinematic}_e8_2_baseline_" + ATTEMPT_MARKER + "_{phi}_{epsilon}_yield_" + kind + ".dat", "file")
      for kind in ("data", "simc")]]}
for _eps in ("lowe", "highe"):
    for _key, _name, _kind in (
        ("full_analysis", f"kaon_FullAnalysis_{KINEMATIC}_{_eps}.json", "json"),
        ("ledger_json", f"kaon_FullAnalysis_{KINEMATIC}_{_eps}_correction_ledger_no_empirical_residual.json", "json"),
        ("ledger_csv", f"kaon_FullAnalysis_{KINEMATIC}_{_eps}_correction_ledger_no_empirical_residual.csv", "file"),
        ("mm_support", f"kaon_xsect_support_{KINEMATIC}_{_eps}.npz", "file")):
        ARTIFACTS["global"].append(declaration(_key + "_" + _eps, _name, _kind))


def strict_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, "duplicate_json_key:" + key)
            result[key] = value
        return result
    def constant(value):
        raise ValueError("nonfinite_json:" + value)
    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=pairs, parse_constant=constant)


def names(profile):
    entries = [("global", entry, {}) for entry in profile["artifacts"]["global"]]
    entries += [(row["phi"] + "_" + row["epsilon"], entry, row)
                for row in SETTINGS for entry in profile["artifacts"]["settings"]]
    return {prefix + "/" + entry["basename_template"].format(kinematic=KINEMATIC, **row): entry
            for prefix, entry, row in entries}


def attempt_identity(output):
    """Bind evidence to the exact canonical requested ZIP, without lossy slugs.

    A full digest distinguishes equal stems at different transfer roots. The
    bounded ASCII stem keeps all resulting evidence basenames below NAME_MAX.
    This function does not create, remove or repair any external path.
    """
    output = Path(output)
    require(output.is_absolute() and output == output.resolve() and output.suffix == ".zip" and
            re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", output.stem) is not None,
            "requested_attempt_identity_invalid")
    digest = hashlib.sha256(output.as_posix().encode("utf-8")).hexdigest()
    return {"attempt_id": output.stem + "-" + digest, "requested_output": output.as_posix()}


def artifact_name(profile, key, *, phi=None, epsilon=None):
    scope = "global" if phi is None else "settings"
    entries = [e for e in profile["artifacts"][scope] if e["key"] == key]
    require(len(entries) == 1, "attempt_artifact_declaration_invalid:" + key)
    return entries[0]["basename_template"].format(kinematic=KINEMATIC, phi=phi, epsilon=epsilon)


def resolved_profile(repo, commit, output):
    require(re.fullmatch(r"[0-9a-f]{40}", commit) is not None, "source_commit_requires_full_sha")
    profile = collector.load_validation_profile(Path(repo) / PROFILE)
    require(profile["settings"] == SETTINGS and profile["artifacts"] == ARTIFACTS,
            "profile_inventory_invalid")
    require(profile["validation_profile"] == PROFILE_ID and
            profile["collection_mode"] == "generic_artifacts" and
            profile["source_identity"] == SOURCE_IDENTITY, "profile_identity_invalid")
    profile = deepcopy(profile)
    profile["source_identity"]["required_analysis_commit"] = commit
    identity = attempt_identity(output)
    for scope in ("global", "settings"):
        for entry in profile["artifacts"][scope]:
            entry["basename_template"] = entry["basename_template"].replace(ATTEMPT_MARKER, identity["attempt_id"])
    return profile


def attempt_preflight(outdir, output, profile):
    """Reject collisions before a source worktree or analysis is launched."""
    identity = attempt_identity(output)
    owned = [outdir / artifact_name(profile, "run_summary"),
             outdir / (identity["attempt_id"] + ".log"),
             outdir / (identity["attempt_id"] + "-gate-status.json"), output]
    owned += [outdir / artifact_name(profile, "yield_" + kind, **row)
              for row in SETTINGS for kind in ("data", "simc")]
    owned += [output.with_name(output.stem + suffix) for suffix in
              (".log", "-run-summary.json", "-gate-status.json")]
    require(len(set(owned)) == len(owned), "attempt_output_alias")
    for path in owned:
        require(not os.path.lexists(path), "attempt_evidence_already_exists:" + path.name)
    return identity


def snapshot(repo):
    state = isolation.snapshot(repo)
    state["refs"] = git(repo, "show-ref")
    state["index"] = git(repo, "ls-files", "--stage", "-z")
    files = git(repo, "ls-files", "--cached", "--others", "--exclude-standard", "-z").split("\0")
    state["content"] = {}
    for name in set(files) - {""}:
        path = Path(repo) / name
        state["content"][name] = ("link:" + os.readlink(path) if path.is_symlink() else
                                  sha256(path) if path.is_file() else "missing")
    return state


def preflight(repo, commit):
    require(re.fullmatch(r"[0-9a-f]{40}", commit) is not None, "source_commit_requires_full_sha")
    state = snapshot(repo)
    require(state["branch"] == "test", "wrong_branch")
    require(state["head"] == commit, "wrong_head")
    require(state["origin_test"] == commit, "source_not_observed_pushed")
    require(commit != BASE_HEAD, "owner_source_not_committed")
    git(repo, "merge-base", "--is-ancestor", BASE_HEAD, commit)
    for line in state["porcelain"].splitlines():
        name = line[3:]
        require(name in isolation.FARM_OUTPUTS or name.startswith("OUTPUT/") or
                (line[:2] == "??" and re.fullmatch(r"kaonlt_review(?:\([^/]+\))?\.diff", name)),
                "dirty_gate_source:" + line)
    git(repo, "diff", "--check", "HEAD", "--", ".",
        *(":(exclude)" + name for name in sorted(isolation.FARM_OUTPUTS)))
    return state


def preservation(repo, before):
    require(snapshot(repo) == before, "ordinary_checkout_identity_or_bytes_changed")
    return {"passed": True, "refs_index_status_and_content_checked": True}


def source_preflight(repo, commit, profile, output):
    with owned_worktree(repo, commit, "source") as source:
        require(resolved_profile(source, commit, output) == profile, "collector_profile_mismatch")
        with collection_module(source) as module, redirect_stdout(sys.stderr):
            _, checks = module.collect_source_checks(source, required_analysis_commit=commit,
                                                     allowed_committed_files=[])
            for check in checks:
                require(check["returncode"] == 0, "collector_source_check_failed:" + check["name"])
            ancestor, _, unexpected = module._committed_identity(checks, [], ["docs/memory/"])
            require(ancestor and not unexpected, "collector_source_identity_invalid")
    return checks


def fresh(path, started, before):
    require(path.is_file() and not path.is_symlink() and path.stat().st_size > 0,
            "artifact_missing_or_symlink:" + path.name)
    record = {"sha256": sha256(path), "bytes": path.stat().st_size}
    require(path.stat().st_mtime_ns >= started and
            (path.stat().st_mtime_ns, record["bytes"], record["sha256"]) != before.get(path.name),
            "artifact_stale:" + path.name)
    return record


def parent_states(hist, t_edges):
    parents = hist.get("pion_t_amplitude_table")
    require(hist.get("pion_t_parent_collection_frozen") is True and
            isinstance(parents, list) and len(parents) == 3, "frozen_parent_inventory_missing")
    states = []
    for t, parent in enumerate(parents):
        require(parent.get("t_bin") == t + 1 and parent.get("t_edges") == t_edges[t:t + 2],
                "frozen_parent_geometry_invalid")
        status = parent.get("diagnostic_application_status", {})
        fallback = status.get("fallback_mode")
        # A detached diagnostic proposal can be unavailable even when the
        # frozen production evaluation has a defined skip/zero policy.
        evaluation = status.get("production_evaluation")
        if evaluation == "accepted":
            states.append("valid")
        elif evaluation == "rejected" and fallback == "skip_bin":
            states.append("skip_bin")
        elif evaluation == "rejected" and fallback == "zero":
            states.append("zero")
        elif evaluation == "rejected" and fallback == "single_scale":
            states.append("valid")
        else:
            raise ValueError("frozen_parent_policy_not_verified:t" + str(t + 1))
    return states


def verify_pages(payload, phi, eps, hist, states):
    require({"phi": phi, "epsilon": eps} in SETTINGS, "setting_not_canonical")
    expected = {"kinematic_token": KINEMATIC, "epsilon_setting": eps[:-1],
                "epsilon_filename_token": eps, "phi_setting": phi, "particle_type": "kaon"}
    require(payload.get("schema_version") == "full_background_subtraction_page_manifest/v1" and
            all(payload.get("setting", {}).get(k) == v for k, v in expected.items()), "page_identity_invalid")
    require(payload.get("pdf_basename") == isolation.pdf_name(phi, eps) and
            payload.get("renderer_failures") == [], "page_pdf_or_renderer_invalid")
    pages = payload.get("pages", [])
    counts = Counter(p.get("page_id") for p in pages)
    require(all(counts[key] == 1 for key in baseline.PAGE_IDS), "e8_2_page_missing_or_duplicate")
    positions = []
    for pos, page in enumerate(pages):
        if page.get("page_id") not in baseline.PAGE_IDS:
            continue
        n, semantic = baseline.PAGE_IDS[page["page_id"]]
        require(page.get("scope") == f"t{n}" and type(page.get("t_index")) is int and
                page["t_index"] == n - 1 and page.get("semantic_stage") == semantic and
                page.get("setting") == phi and page.get("epsilon") == eps[:-1] and
                page.get("t_edges") == hist["t_bins"][n - 1:n + 1], "e8_2_parent_identity_invalid")
        require(page.get("represented_phi_inventory") == [
            {"phi_index": k, "phi_edges": hist["phi_bins"][k:k + 2]} for k in range(9)],
            "e8_2_phi_geometry_invalid")
        invalid = page.get("invalid_unavailable_children")
        if states[n - 1] == "skip_bin":
            require(isinstance(invalid, list) and [p.get("phi_index") for p in invalid] == list(range(9)) and
                    all(isinstance(p.get("reason"), str) and p["reason"] for p in invalid),
                    "explicit_frozen_skip_inventory_invalid")
        else:
            require(invalid == [], "e8_2_invalid_populated_children")
        require(page.get("authoritative") is False and page.get("presentation_only") is True,
                "e8_2_presentation_flags_invalid")
        positions.append(pos)
    handoff = [i for i, p in enumerate(pages) if p.get("page_id") == "full_background.e8.handoff"]
    require(len(handoff) == 1 and handoff[0] > max(positions), "e8_handoff_missing_or_misordered")
    return len(pages)


def table_name(worktree, inp, phi, kind):
    pol = "pl" if float(inp["POL"]) > 0 else "mn"
    stem = f"yield_{kind}.{pol}_Q44W274_{float(inp['EPSVAL']) * 100:.0f}_"
    if phi == "Center":
        return stem + "+0000.dat"
    runs = [int(x) for x in inp["runNum" + phi].split()]
    for i, run in enumerate(runs):
        log = Path(worktree) / "log" / f"{phi}_kaon_{run}_{inp['OutFilename'].replace('FullAnalysis_', '')}.log"
        if log.is_file():
            angle = float(f"{abs(float(inp['pThetaValCenter'][i]) - float(inp['pThetaVal' + phi][i])):.3f}")
            return stem + ("+" if phi == "Left" else "-") + str(int(angle * 1000)) + ".dat"
    raise ValueError("yield_filename_run_log_missing:" + phi)


def read_table(path, expected):
    rows = {}
    order = []
    for line in path.read_text(encoding="utf-8").splitlines():
        tokens = line.split()
        require(len(tokens) == 4 and all(re.fullmatch(r"-?\d+\.\d{4}", x) for x in tokens[:2]) and
                all(re.fullmatch(r"\d+", x) for x in tokens[2:]), "yield_table_format_invalid")
        # physics_lists.py writes one-based (phi_bin, t_bin); NPZ uses
        # zero-based [t_index, phi_index] children.
        phi_bin, t_bin = int(tokens[2]), int(tokens[3])
        require(1 <= phi_bin <= 9 and 1 <= t_bin <= 3, "yield_producer_coordinate_out_of_range")
        coordinate = (t_bin - 1, phi_bin - 1)
        require(coordinate not in rows, "duplicate_yield_coordinate")
        values = [float(x) for x in tokens[:2]]
        require(all(np.isfinite(values)) and values[1] >= 0, "yield_table_undefined_error")
        rows[coordinate] = {"yield": values[0], "error": values[1], "yield_text": tokens[0], "error_text": tokens[1],
                            "producer_phi_bin": phi_bin, "producer_t_bin": t_bin}
        order.append(coordinate)
    require(order == sorted(expected), "yield_coordinate_inventory_or_order_invalid")
    return rows


def integrate(values):
    # Match the existing sequential normal-bin Python float accumulation.
    total = 0.0
    for value in values:
        total += float(value)
    return total


def check_simc(simc, simc_errors, simc_row):
    raw_simc = integrate(simc)
    ys = max(0.0, raw_simc)
    es = 0.0 if raw_simc < 0 else integrate(float(x) ** 2 for x in simc_errors) ** 0.5
    require(format(ys, ".4f") == simc_row["yield_text"] and
            format(es, ".4f") == simc_row["error_text"], "simc_yield_support_mismatch")
    return {"simc_signed_normal_bin_integral": raw_simc,
            "simc_nonnegative_yield": ys, "text_quantization": "producer .4f, exact formatted comparison",
            "simc": simc_row, "simc_numerical_gate": "PASS"}


def check_cell(data, simc, simc_errors, data_row, simc_row):
    yd = integrate(data)
    require(format(yd, ".4f") == data_row["yield_text"], "data_yield_support_mismatch")
    return {"data_normal_bin_integral": yd, "data": data_row, "numerical_gate": "PASS",
            **check_simc(simc, simc_errors, simc_row)}


def audit_outputs(worktree, outdir, profile, started, before):
    records, cells, metadata, used_tables = {}, [], {}, set()
    copied_names = {Path(p).name for p, e in names(profile).items() if e["key"].startswith("yield_")}
    for archive, entry in names(profile).items():
        name = Path(archive).name
        if name == artifact_name(profile, "run_summary") or name in copied_names:
            continue
        records[name] = fresh(outdir / name, started, before)
        if entry["kind"] == "json":
            strict_json(outdir / name)
    for eps, phis in (("lowe", ["Left", "Center"]), ("highe", ["Left", "Center", "Right"])):
        full = strict_json(outdir / f"kaon_FullAnalysis_{KINEMATIC}_{eps}.json")
        inp = full.get("inpDict", {})
        require(all(inp.get(k) == v for k, v in {"ParticleType": "kaon", "EPSSET": eps[:-1],
                "Q2": "4p4", "W": "2p74", "OutFilename": f"FullAnalysis_{KINEMATIC}_{eps}"}.items()),
                "full_analysis_identity_invalid")
        require(int(inp.get("iter_num", 0) or 0) == 0, "unsupported_or_ambiguous_iteration")
        require(inp.get("bg_active_profile") == "no_empirical_residual" and
                inp.get("bg_stat_scale1") == 0 and inp.get("bg_stat_scale2") == 0,
                "full_analysis_runtime_profile_invalid")
        histlist = full.get("histlist", [])
        actual_phis = [h.get("phi_setting") for h in histlist]
        require(Counter(actual_phis) == Counter(phis), "full_analysis_setting_inventory_invalid")
        ledgerbase = f"kaon_FullAnalysis_{KINEMATIC}_{eps}_correction_ledger_no_empirical_residual"
        ledger = strict_json(outdir / (ledgerbase + ".json"))
        require(all(ledger.get(k) == v for k, v in {"active_profile": "no_empirical_residual",
                "particle_type": "kaon", "epsset": eps[:-1], "q2": "4p4", "w": "2p74",
                "outfilename": inp["OutFilename"]}.items()) and
                Counter(h.get("phi_setting") for h in ledger.get("settings", [])) == Counter(phis), "ledger_identity_invalid")
        with (outdir / (ledgerbase + ".csv")).open(encoding="utf-8", newline="") as handle:
            require(Counter(r.get("phi_setting") for r in csv.DictReader(handle) if r.get("row_kind") == "setting_total") == Counter(phis),
                    "ledger_csv_setting_inventory_invalid")
        with np.load(outdir / f"kaon_xsect_support_{KINEMATIC}_{eps}.npz", allow_pickle=False) as support:
            require(support["settings"].tolist() == actual_phis, "support_setting_inventory_invalid")
            t, phi_edges, mm = [support[key].tolist() for key in ("t_bins", "phi_bins", "mm_edges")]
            require(len(t) == 4 and len(phi_edges) == 10 and len(mm) >= 2 and
                    all(np.all(np.isfinite(a)) and np.all(np.diff(a) > 0) for a in (t, phi_edges, mm)), "support_axes_invalid")
            metadata[eps] = {"iteration": 0, "mm_edges": mm, "t_edges": t, "phi_edges": phi_edges,
                "persisted_input": {key: value for key, value in inp.items()
                    if key in ("Q2", "W", "ParticleType", "EPSSET", "EPSVAL", "POL", "OutFilename", "iter_num")
                    or (not key.startswith("_") and any(token in key.lower() for token in ("model", "parameter", "iter_")))},
                "model_identity_claim": "only source-persisted input metadata; no verified model freeze",
                "settings": {}}
            window = [float(inp["mm_min"]), float(inp["mm_max"])]
            require(len(mm) == 101 and [mm[0], mm[-1]] == window and ledger.get("mm_cut_window") == window,
                    "support_cut_selected_mm_window_mismatch")
            metadata[eps]["cut_selected_mm_window"] = window
            for hist in histlist:
                phi = hist["phi_setting"]
                require(hist.get("t_bins") == t and hist.get("phi_bins") == phi_edges, "support_source_axes_mismatch")
                states = parent_states(hist, t)
                pdf = outdir / isolation.pdf_name(phi, eps)
                require(pdf.read_bytes().startswith(b"%PDF-"), "procedure_pdf_invalid")
                page_count = verify_pages(strict_json(pdf.with_name(pdf.stem + "-manifest.json")), phi, eps, hist, states)
                e84 = [p for p in strict_json(pdf.with_name(pdf.stem + "-manifest.json"))["pages"]
                       if str(p.get("page_id", "")).startswith("full_background.e8_4")]
                metadata[eps]["settings"][phi] = {"parent_states": states, "frozen_parents": hist["pion_t_amplitude_table"],
                    "page_count": page_count, "e8_4_persisted_page_records": e84,
                    "e8_4_reason_claim": "persisted page records; no Method-A acceptance",
                    "e8_4_unavailable_reasons": [p.get("reason") for p in e84 if p.get("page_id") == "full_background.e8_4.unavailable"]}
                require(e84 and all(isinstance(p.get("reason"), str) and p["reason"] for p in e84
                    if p.get("page_id") == "full_background.e8_4.unavailable"), "e8_4_diagnostic_reason_not_persisted")
                arrays = {}
                for kind in ("data", "simc"):
                    for field in ("values", "errors"):
                        key = f"{kind}_mm_{phi.lower()}_{field}"
                        arr = support[key]
                        require(arr.shape == (3, 9, len(mm) - 1) and np.all(np.isfinite(arr)) and
                                (field != "errors" or np.all(arr >= 0)), "support_matrix_invalid:" + key)
                        arrays[kind + "_" + field] = arr
                expected_data = {(j, k) for j in range(3) for k in range(9) if states[j] != "skip_bin"}
                tables = {}
                for kind in ("data", "simc"):
                    source = Path(worktree) / "src/kaon/yields" / table_name(worktree, inp, phi, kind)
                    require(source not in used_tables, "duplicate_yield_source_identity")
                    used_tables.add(source)
                    identity = fresh(source, started, {})
                    tables[kind] = read_table(source, expected_data if kind == "data" else {(j, k) for j in range(3) for k in range(9)})
                    copied = artifact_name(profile, "yield_" + kind, phi=phi, epsilon=eps)
                    destination = outdir / copied
                    require(not os.path.lexists(destination), "yield_evidence_already_exists:" + copied)
                    with source.open("rb") as src, destination.open("xb") as dst:
                        shutil.copyfileobj(src, dst)
                    copied_identity = fresh(destination, started, before)
                    require(copied_identity == identity, "yield_copy_hash_mismatch")
                    records[copied] = {**identity, "source_relative_basename": str(source.relative_to(worktree)),
                        "copied_basename": copied, "source_mtime_ns": source.stat().st_mtime_ns,
                        "copied_mtime_ns": destination.stat().st_mtime_ns}
                for j in range(3):
                    for k in range(9):
                        cell = {"setting": phi, "epsilon": eps, "t_index": j, "phi_index": k,
                                "parent_state": states[j], **{key: arr[j, k].tolist() for key, arr in arrays.items()}}
                        if states[j] == "skip_bin":
                            cell.update(data=None, numerical_gate="excluded_frozen_skip_bin",
                                **check_simc(arrays["simc_values"][j, k], arrays["simc_errors"][j, k], tables["simc"][(j, k)]))
                        else:
                            cell.update(check_cell(arrays["data_values"][j, k], arrays["simc_values"][j, k],
                                                   arrays["simc_errors"][j, k], tables["data"][(j, k)], tables["simc"][(j, k)]))
                        cells.append(cell)
    return records, cells, metadata


def verify_zip(path, commit, records, profile):
    with zipfile.ZipFile(path) as archive:
        require(archive.testzip() is None, "zip_integrity_failed")
        manifest = json.loads(archive.read("manifest.json"))
        require(manifest.get("complete") is True and manifest.get("errors") == [], "bundle_incomplete")
        require(manifest.get("git_head") == commit and manifest.get("required_analysis_commit") == commit and
                manifest.get("validation_profile") == PROFILE_ID and manifest.get("requested_settings") == SETTINGS and
                manifest.get("requested_kinematic") == KINEMATIC, "bundle_source_or_inventory_invalid")
        rows = manifest.get("settings", [])
        require([{k: r.get(k) for k in ("phi", "epsilon")} for r in rows] == SETTINGS, "bundle_settings_invalid")
        require(set(manifest.get("global_artifacts", {})) == {e["key"] for e in profile["artifacts"]["global"]},
                "bundle_global_keys_invalid")
        for row in rows:
            require(row.get("kinematic") == KINEMATIC and
                    set(row.get("artifacts", {})) == {e["key"] for e in profile["artifacts"]["settings"]},
                    "bundle_setting_keys_invalid")
        for row, scope, prefix in [(manifest, "global", "global")] + [
                (row, "settings", row["phi"] + "_" + row["epsilon"]) for row in rows]:
            entries_by_key = row["global_artifacts"] if scope == "global" else row["artifacts"]
            for declaration in profile["artifacts"][scope]:
                basename = declaration["basename_template"].format(kinematic=KINEMATIC,
                    **({} if scope == "global" else {k: row[k] for k in ("phi", "epsilon")}))
                require(entries_by_key[declaration["key"]].get("archive_path") == prefix + "/" + basename,
                        "bundle_artifact_scope_invalid")
        declared = names(profile)
        entries = list(manifest.get("global_artifacts", {}).values()) + [e for row in rows for e in row.get("artifacts", {}).values()]
        require(len(entries) == len(declared), "bundle_artifact_count_invalid")
        seen = set()
        for entry in entries:
            name = entry["archive_path"]
            require(name in declared and name not in seen, "bundle_artifact_inventory_invalid")
            seen.add(name)
            raw = archive.read(name)
            record = records[Path(name).name]
            require(hashlib.sha256(raw).hexdigest() == entry["sha256"] == record["sha256"] and
                    len(raw) == entry["byte_size"] == record["bytes"], "bundle_artifact_hash_mismatch")
        expected = set(declared) | {"manifest.json", "source_state.txt", "source_checks.txt", "global/"} | {
            row["phi"] + "_" + row["epsilon"] + "/" for row in SETTINGS}
        require(set(archive.namelist()) == expected and len(archive.namelist()) == len(expected), "bundle_archive_inventory_invalid")


class GateStatus(isolation.accepted.GateStatus):
    def __init__(self, outdir, output, commit):
        identity = attempt_identity(output)
        self.path = outdir / (identity["attempt_id"] + "-gate-status.json")
        now = datetime.now(timezone.utc).isoformat()
        self.record = {"schema_version": "e8_2_baseline_five_yield_gate_status/v1",
            **identity, "source_commit": commit, "kinematic": KINEMATIC, "settings": SETTINGS,
            "expected_zip_path": str(output), "started_at_utc": now, "updated_at_utc": now,
            "status": "running", "stage": "preflight", "failure_reason": None,
            "analysis_started": False, "analysis_completed": False,
            "artifact_verification_completed": False, "collector_source_preflight_completed": False,
            "collection_completed": False, "zip_verification_completed": False}
        self._persist(initial=True)


def execute_gate(repo, outdir, output, commit):
    repo, outdir, output = [Path(p).resolve() for p in (repo, outdir, output)]
    # Validate the exact attempt profile and all owned destinations before
    # publishing status, so a collision cannot overwrite even a failed receipt.
    profile = resolved_profile(repo, commit, output)
    attempt = attempt_preflight(outdir, output, profile)
    status = GateStatus(outdir, output, commit)
    original = overlay = None
    try:
        require(not os.path.lexists(output) and output.parent.is_dir(), "output_zip_exists_or_parent_missing")
        original = preflight(repo, commit)
        before = {p.name: (p.stat().st_mtime_ns, p.stat().st_size, sha256(p))
                  for name in names(profile) if (p := outdir / Path(name).name).is_file()}
        log = outdir / (attempt["attempt_id"] + ".log")
        require(not os.path.lexists(log), "run_log_already_exists")
        checks = source_preflight(repo, commit, profile, output)
        status.update(collector_source_preflight_completed=True, collector_source_checks=checks)
        try:
            with owned_worktree(repo, commit, "analysis") as worktree:
                status.update("path_isolation")
                env, overlay = prepare_runtime_overlay(repo, worktree)
                status.update(ltsep_runtime_overlay=overlay)
                probes = probe_paths(worktree, env, overlay["baseline_paths"])
                require(outdir == (Path(probes[0]["paths"]["VOLATILEPATH"]) / "OUTPUT/Analysis/KaonLT").resolve(),
                        "analysis_artifact_root_mismatch")
                external = external_symlink_preflight(worktree, probes[0]["paths"])
                identity = {"path": str(worktree), "head": commit, "detached": True, "owner_created": True}
                started = time.time_ns()
                status.update("analysis", analysis_started=True, analysis_worktree=identity,
                              isolation_probes=probes, external_symlink_preflight=external)
                rc = run_analysis(worktree, log, env)
                status.update("completion_markers", analysis_returncode=rc)
                verify_completion(log, rc)
                status.update("yield_artifact_audit", analysis_completed=True)
                records, cells, metadata = audit_outputs(worktree, outdir, profile, started, before)
                status.update(artifact_verification_completed=True)
        finally:
            try:
                status.update(ordinary_checkout_preservation=preservation(repo, original))
            finally:
                if overlay is not None:
                    status.update(installed_ltsep_preservation=verify_ltsep_preservation(overlay))
        status.update(worktree_cleanup_completed=True)
        summary = {"schema_version": "e8_2_baseline_five_yield_receipt/v1", "source_commit": commit,
                   **attempt, "analysis_command": COMMAND, "analysis_returncode": rc, "log_sha256": sha256(log),
                   "settings": SETTINGS, "active_profile": "no_empirical_residual", "boundaries": BOUNDARIES,
                   "evidence_label": "SOURCE VERIFIED", "artifacts": records, "cells": cells, "epsilon_metadata": metadata,
                   "analysis_worktree": identity, "isolation_probes": probes, "external_symlink_preflight": external,
                   "collector_source_checks": checks, "ordinary_checkout_preservation": status.record["ordinary_checkout_preservation"],
                   "installed_ltsep_preservation": status.record["installed_ltsep_preservation"], "worktree_cleanup_completed": True}
        summary_name = artifact_name(profile, "run_summary")
        summary_path = outdir / summary_name
        with summary_path.open("x", encoding="utf-8") as handle:
            json.dump(summary, handle, allow_nan=False, separators=(",", ":"))
            handle.write("\n")
        records[summary_name] = {"sha256": sha256(summary_path), "bytes": summary_path.stat().st_size}
        try:
            with tempfile.TemporaryDirectory(prefix="kaonlt-baseline-five-profile-") as directory:
                effective = Path(directory) / "profile.json"
                effective.write_text(json.dumps(profile), encoding="utf-8")
                with owned_worktree(repo, commit, "source") as source:
                    with collection_module(source) as module, redirect_stdout(sys.stderr):
                        result = module.collect_validation_bundle(outdir=outdir, kinematic=KINEMATIC, output=output,
                                                                 profile_path=effective, repo_root=source)
                    require(result["returncode"] == 0, "collection_failed")
                    status.update("verify_zip", collection_completed=True)
                    verify_zip(output, commit, records, profile)
                    status.update(zip_verification_completed=True)
        finally:
            try:
                status.update(final_ordinary_checkout_preservation=preservation(repo, original))
            finally:
                status.update(final_installed_ltsep_preservation=verify_ltsep_preservation(overlay))
        isolation.deliver_evidence(status, log, summary_path, output)
        return output
    except Exception as exc:
        status.update(status="failed", failure_reason=str(exc))
        # Preserve the primary checkout even if source/path preflight failed
        # before the analysis-stage finally block was reached.
        try:
            if original is not None:
                status.update(final_ordinary_checkout_preservation=preservation(repo, original))
        finally:
            if overlay is not None:
                status.update(final_installed_ltsep_preservation=verify_ltsep_preservation(overlay))
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        print(execute_gate(args.repo, args.outdir, args.output, args.source_commit))
    except Exception as exc:
        print("Baseline-five gate failed: " + str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
