"""Default-off read-only cache observer. No ROOT import or event-tree iteration."""
from __future__ import annotations

import bisect
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess

ENABLE_ENV = "LT_ANALYSIS_PION_LOCAL_T_EVENT_VALIDATION"
SCHEMA = "pion_local_t_event_validation_packet/v1"
SOURCES = ("prompt", "rand", "dummy", "dummy_rand")
# Hash/size projections from the accepted A/B JSONs at 6d13abb4681b938c56def2a1f455244b0880096c.
# These pins allow capture to reject different physical inputs without a C file on the farm.
ACCEPTED_AB_JSON_SHA256 = {
    "Center": "7e34ec134ef28378ffebf1686c47915a266478531d9910bbdd211f5392cb1f24",
    "Left": "f5c8cf2483ce77d74dc2e2a497d8b7c5f1b0aa3c96b4d76c3964d813e675b812",
}
ACCEPTED_INPUT_IDENTITIES = {
    "Center": {
        "data": {"sha256": "18272b333282e0e0d811e80d07c1af86a6588601dd5c64951304f65296b4dff6", "bytes": 845758483},
        "dummy": {"sha256": "0459e71c9e0f9354df997df49b2fe35354f715449750acb9cce84d4aede7c3b1", "bytes": 10470121},
    },
    "Left": {
        "data": {"sha256": "63f3a320c3d28da036e572d6506e29107e580fba14e4fbcd237a5577569ce71f", "bytes": 1194686097},
        "dummy": {"sha256": "1376574c749b87b8bc7c6d9d2c61e1c0c77c4e58a633f26285c0b25e5f4d18f7", "bytes": 8599072},
    },
}
FOLD_POLICY = {"seed": "kaonlt-pion-event-validation-2026-10-21-v1", "count": 5,
               "identity": "source-file SHA256 + source label + tree-local entry index",
               "status": "partition only; full-data models contaminate every evaluation fold"}


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def encoded(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def file_identity(path):
    path = Path(path).resolve(strict=True)
    h = hashlib.sha256()
    before = path.stat()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    after = path.stat()
    require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns), "source_changed_during_hash")
    return {"sha256": h.hexdigest(), "bytes": after.st_size}


def backing_file_path(source, role):
    """Resolve a tree's current TFile or an existing TFile's owning file.

    GetName on a TTree names a tree, never evidence of its backing filename.
    No file opening, tree traversal, directory change or source mutation occurs.
    """
    reason = "missing_backing_file_provenance:" + role
    getter = getattr(source, "GetCurrentFile", None)
    if not callable(getter):
        getter = getattr(source, "GetFile", None)
    require(callable(getter), reason)
    try:
        backing = getter()
        require(bool(backing), reason)
        require(backing.IsOpen() and not backing.IsZombie(), reason)
        name = backing.GetName()
        require(isinstance(name, str) and bool(name.strip()), reason)
        path = Path(name).resolve(strict=True)
        require(path.is_file(), reason)
    except (AttributeError, OSError, ReferenceError, TypeError, ValueError) as exc:
        raise ValueError(reason) from exc
    return path


def verify_accepted_input_identity(identities, setting):
    require(setting in ACCEPTED_INPUT_IDENTITIES, "accepted_AB_setting_missing")
    for role in ("data", "dummy"):
        require(identities.get(role) == ACCEPTED_INPUT_IDENTITIES[setting][role],
                "accepted_AB_input_identity_mismatch:" + setting + ":" + role)


def resolved_production_weight_policy(inp, parent):
    """Observe the same resolvers used by the existing yield application.

    Frozen parents need not carry a parent_fit_configuration dictionary. When
    retained application diagnostics exist, their recorded policy must agree.
    """
    from background_config import (
        resolve_particle_subtraction_weight_denominator_floor,
        resolve_particle_subtraction_weight_clip_bounds,
    )
    floor = resolve_particle_subtraction_weight_denominator_floor(inp)
    bounds = list(resolve_particle_subtraction_weight_clip_bounds(inp))
    require(math.isfinite(floor) and floor >= 0 and len(bounds) == 2
            and all(v is None or math.isfinite(v) for v in bounds)
            and (None in bounds or bounds[0] <= bounds[1]), "production_weight_policy_invalid")
    expected = {"denom_floor": floor, "clip_min": bounds[0], "clip_max": bounds[1]}
    for key in ("final_diagnostic_application_result", "diagnostic_application_result",
                "proposed_diagnostic_application_result"):
        payload = parent.get(key)
        diagnostics = payload.get("diagnostics") if isinstance(payload, dict) else None
        if isinstance(diagnostics, dict) and any(k in diagnostics for k in expected):
            require(all(k in diagnostics and diagnostics[k] == v for k,v in expected.items()),
                    "production_weight_configuration_disagreement:" + key)
    return {"weight_denominator_floor": floor, "weight_clip_bounds": bounds}


def fold(row):
    key = [FOLD_POLICY["seed"], row["source_sha256"], row["source_label"], row["entry_index"]]
    return int(digest(key), 16) % FOLD_POLICY["count"]


def bin_index(value, edges, *, final_endpoint=False):
    require(math.isfinite(value), "nonfinite_coordinate")
    i = bisect.bisect_right(edges, value) - 1
    if final_endpoint and value == edges[-1]:
        i -= 1
    return i if 0 <= i < len(edges)-1 else None


def coefficients(data, dummy, windows):
    require(all(math.isfinite(x) and x > 0 for x in (data, dummy, windows)), "source_normalization_invalid")
    return dict(zip(SOURCES, (data, -data/windows, -dummy, dummy/windows)))


def validate_histogram(h):
    require(h.get("status") == "available", "parent_histogram_unavailable")
    edges,values,errors,variance = (h[k] for k in ("edges","values","errors","stored_sumw2"))
    require(len(edges) >= 2 and len(values) == len(errors) == len(variance) == len(edges)+1
            and all(math.isfinite(x) for x in edges+values+errors+variance)
            and all(a < b for a,b in zip(edges,edges[1:]))
            and all(x >= 0 for x in errors+variance), "parent_histogram_layout_or_finiteness")
    require(all(math.isclose(e*e,v,rel_tol=1e-10,abs_tol=1e-18) for e,v in zip(errors,variance)), "parent_histogram_Sumw2")
    packed = [float(len(edges)-1),edges[0],edges[-1],*edges]
    packed += [v for pair in zip(values,errors) for v in pair]
    require(hashlib.sha256(b"".join(struct.pack("!d",v) for v in packed)).hexdigest() == h["content_error_fingerprint"], "parent_histogram_fingerprint")


def validate_packet(packet):
    require(packet.get("schema_version") == SCHEMA and packet.get("diagnostic_only") is True
            and packet.get("production_applied") is False, "packet_schema_or_authority")
    require(packet.get("setting") in ("Center", "Left") and packet.get("kinematic") == "Q4p4W2p74"
            and packet.get("epsilon") == "lowe", "packet_setting")
    require(packet.get("fold_policy") == FOLD_POLICY, "fold_policy_changed")
    require(re.fullmatch(r"[0-9a-f]{40}", packet.get("source_HEAD", "")) is not None, "packet_source_HEAD")
    for role in ("data", "dummy"):
        identity = packet["input_identity"][role]
        require(re.fullmatch(r"[0-9a-f]{64}", identity.get("sha256", "")) is not None
                and type(identity["bytes"]) is int and identity["bytes"] > 0, "input_identity_missing")
    require(packet["input_identity"]["data"]["sha256"] != packet["input_identity"]["dummy"]["sha256"], "source_file_collision")
    for key in ("t_edges", "phi_edges"):
        edges = packet[key]
        require(len(edges) >= 2 and all(math.isfinite(x) for x in edges)
                and all(a < b for a,b in zip(edges, edges[1:])), "packet_edges")
    require(len(packet["t_edges"]) == 4 and len(packet["parents"]) == 3, "packet_parent_inventory")
    require([p["t_bin_index"] for p in packet["parents"]] == [0,1,2], "packet_parent_inventory")
    for parent in packet["parents"]:
        for side in ("pion_input", "kaon_input"):
            validate_histogram(parent[side])
    for population in ("pion", "kaon"):
        rows = packet["events"][population]
        seen = set()
        counts = {s: 0 for s in SOURCES}
        for r in rows:
            require(r.get("population") == population, "cross_population_leakage")
            s = r["source_label"]
            require(s in SOURCES and type(r["entry_index"]) is int and r["entry_index"] >= 0, "source_row_identity")
            key = (s, r["entry_index"])
            require(key not in seen, "duplicate_source_entry")
            seen.add(key); counts[s] += 1
            role = "dummy" if s.startswith("dummy") else "data"
            require(r["source_sha256"] == packet["input_identity"][role]["sha256"], "source_file_identity")
            require(r["coefficient"] == packet["coefficients"][s], "source_coefficient")
            require(r["coordinate_fingerprint"] == packet["coordinate_fingerprint"], "coordinate_mismatch")
            require(all(type(r[k]) is bool for k in ("allcuts", "nommcuts", "proton_rf_accept"))
                    and (r["allcuts"] or r["nommcuts"]), "cut_flags")
            require(all(math.isfinite(r[k]) for k in ("MM", "t", "phi", "coefficient", "post_proton_factor")), "nonfinite_event")
            ti = bin_index(r["t"], packet["t_edges"], final_endpoint=population == "pion")
            pi = bin_index(r["phi"], packet["phi_edges"], final_endpoint=population == "pion")
            require(ti is not None and r["t_index"] == ti and r["phi_index"] == pi, "event_bin_identity")
            require(r["fold"] == fold(r), "fold_identity")
            if population == "pion":
                require(r["post_proton_factor"] == 1. and r["proton_rf_accept"], "proton_on_pion_control")
        require(counts == packet["source_counts"][population], "missing_source_rows")
        require(digest(rows) == packet["event_hashes"][population], "event_hash_mismatch")
        require(all(any(r["t_index"] == ti and r["allcuts"] and r["phi_index"] is not None
                        and abs(r["coefficient"]*r["post_proton_factor"]) > 0 for r in rows)
                    for ti in range(3)), "missing_allcuts_support")
    return packet


def build_packet(hist, inp, t_edges, phi_edges, *, source_head, input_identity, parent_records, source_files):
    """Copy only scalar cache fields; never add fields to a production object."""
    cache = hist["_authoritative_pion_control_source_cache"]
    coeff = coefficients(float(hist["normfac_data"]), float(hist["normfac_dummy"]), float(hist["nWindows"]))
    events = {"pion": [], "kaon": []}
    def append(population, source, entry, mm, t, phi, ti, pi, allcuts, nommcuts, factor=1., rf=True):
        source = "dummy" if source == "dummy_prompt" else source
        role = "dummy" if source.startswith("dummy") else "data"
        row = {"population": population, "source_label": source, "entry_index": int(entry), "source_sha256": input_identity[role]["sha256"],
               "coefficient": coeff[source], "post_proton_factor": float(factor), "proton_rf_accept": bool(rf),
               "MM": float(mm), "t": float(t), "phi": float(phi), "t_index": int(ti),
               "phi_index": None if pi is None or int(pi) < 0 else int(pi),
               "allcuts": bool(allcuts), "nommcuts": bool(nommcuts),
               "coordinate_fingerprint": cache["coordinate_fingerprint"]}
        row["fold"] = fold(row)
        events[population].append(row)
    for r in cache["records"]:
        require(r["coefficient"] == coeff[r["source_label"]] and r["coordinate_fingerprint"] == cache["coordinate_fingerprint"], "authoritative_row_provenance")
        append("pion", r["source_label"], r["entry_index"], r["adj_MM"], r["adj_t"], r["phi_degrees"],
               r["t_index"], r["phi_index"], r["allcuts"], r["nommcuts"])
    for s in SOURCES:
        require(sum(r["source_label"] == s for r in events["pion"]) == cache["source_accounting"][s]["records_inside_canonical_t"], "pion_source_accounting_mismatch")
    for source, section in hist["_yield_data_event_cache"].items():
        count = len(section["adj_MM"])
        require(all(len(section[k]) == count for k in ("source_label", "source_entry_index", "adj_t", "phi_shift_deg",
                    "t_index", "phi_index", "allcuts", "nommcuts", "proton_cleaning_factor", "proton_rf_accept")), "kaon_cache_column_length")
        for i in range(len(section["adj_MM"])):
            append("kaon", str(section["source_label"][i]), section["source_entry_index"][i],
                   section["adj_MM"][i], section["adj_t"][i], section["phi_shift_deg"][i],
                   section["t_index"][i], section["phi_index"][i], section["allcuts"][i], section["nommcuts"][i],
                   section["proton_cleaning_factor"][i], section["proton_rf_accept"][i])
    # The child cache is a projection of exactly the authoritative pion rows.
    child = hist["_yield_sub_event_cache"]
    for s in SOURCES:
        expected = [r for r in events["pion"] if r["source_label"] == s and r["phi_index"] is not None]
        require(len(expected) == len(child[s]["entry_index"]), "child_cache_count_mismatch")
        for r, i in zip(expected, range(len(expected))):
            for key, field in (("entry_index", "entry_index"), ("MM", "adj_MM"), ("t", "adj_t"),
                               ("phi", "phi_degrees"), ("allcuts", "allcuts"), ("nommcuts", "nommcuts")):
                require(r[key] == child[s][field][i], "child_cache_row_mismatch")
    for rows in events.values():
        rows.sort(key=lambda r: (SOURCES.index(r["source_label"]), r["entry_index"]))
    packet = {"schema_version": SCHEMA, "diagnostic_only": True, "production_applied": False,
              "source_HEAD": source_head, "source_files": source_files,
              "setting": hist["phi_setting"], "kinematic": "Q4p4W2p74", "epsilon": "lowe",
              "input_identity": input_identity, "coefficients": coeff, "fold_policy": FOLD_POLICY,
              "t_edges": list(map(float,t_edges)), "phi_edges": list(map(float,phi_edges)),
              "coordinate_fingerprint": cache["coordinate_fingerprint"],
              "coordinate": dict(cache["kaon_data_coordinate"]),
              "pid_fingerprint": cache["physical_pion_control_mask_fingerprint"],
              "parents": parent_records, "events": events,
              "source_counts": {p: {s: sum(r["source_label"] == s for r in rows) for s in SOURCES} for p,rows in events.items()},
              "event_hashes": {p: digest(rows) for p,rows in events.items()},
              "policies": {"duplicates": "reject within population/source/tree-entry; PID populations separate",
                           "source_collision": "data/dummy files must have distinct SHA256; no cross-tree physical-event identity inferred",
                           "canonical": "pion includes final endpoint; yield cache excludes final endpoint; phi outside acceptance retained as null child",
                           "MM": "ROOT [low, high) normal bins; flow weights zero, flows explicitly retained",
                           "serialization": "sorted JSON, UTF-8, finite Python binary64 round-trip decimals"}}
    require(input_identity["data"]["sha256"] != input_identity["dummy"]["sha256"], "source_file_collision")
    validate_packet(packet)
    return packet


def observe(hist, inp, t_edges, phi_edges, outpath):
    """Runtime adapter called only after the single opt-in guard."""
    if os.environ.get(ENABLE_ENV) != "1":
        return
    if (str(inp.get("Q2")), str(inp.get("W")), str(inp.get("ParticleType")), str(inp.get("EPSSET"))) not in (
            ("4p4", "2p74", "kaon", "low"), ("4p4", "2p74", "kaon", "lowe")) or hist.get("phi_setting") not in ("Center", "Left"):
        return
    from pion_background_chain_audit import digest as parent_digest, histogram
    root = Path(__file__).resolve().parents[1]
    source_paths = {role: backing_file_path(hist[key], role)
                    for role,key in (("data","InFile_DATA"),("dummy","InFile_DUMMY"))}
    identities = {role: file_identity(path) for role,path in source_paths.items()}
    # Reject different data/dummy inputs before parent serialization or packet output.
    verify_accepted_input_identity(identities, hist["phi_setting"])
    parents = []
    for p in hist["_pion_t_bin_parent_results"]:
        parents.append({k: p[k] for k in ("t_bin_index", "t_edges", "pion_parent_id", "canonical_interval_pair_id", "canonical_interval_pair_hash", "coordinate_fingerprint")})
        parents[-1].update(baseline_parent_fingerprint=parent_digest(p),
                           pion_input=histogram(p["H_pion_control"]), kaon_input=histogram(p["H_proton_cleaned_final_rf"]),
                           production_weight_configuration=resolved_production_weight_policy(inp, p))
    from pion_t_template_comparison import REQUIRED_SOURCE_FILES
    files = (*REQUIRED_SOURCE_FILES, "src/binning/calculate_yield.py", "testing/pion_local_t_event_validation_capture.py",
             "src/binning/binning_helpers.py", "src/cuts/pion_component_subtraction.py", "src/cuts/proton_contamination_weights.py")
    packet = build_packet(hist, inp, t_edges, phi_edges,
        source_head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        input_identity=identities, parent_records=parents, source_files={p:file_identity(root/p) for p in files})
    packet["source_input_paths"] = {role: str(path) for role,path in source_paths.items()}
    destination = Path(outpath) / (hist["phi_setting"] + "_kaon_Q4p4W2p74_lowe_pion-local-t-events.json")
    require(not os.path.lexists(destination), "packet_destination_exists")
    with destination.open("xb") as stream:
        stream.write(encoded(packet))
