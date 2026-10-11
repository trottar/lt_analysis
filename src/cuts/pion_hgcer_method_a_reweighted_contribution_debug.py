"""Detached descriptive aggregates from the exact October 4 historical lineage.

No production caller, authority override, event persistence or new correction.
Shared F.4 reproduction owns the transient row order and factors. All squared
moments and Neff are descriptive, not yield errors or covariance estimates.
"""

from __future__ import annotations

import bisect
import copy
import hashlib
import json
import math
from collections.abc import Mapping

import numpy as np

import pion_hgcer_method_a_acceptance_map as f3
import pion_hgcer_method_a_parent_preserving_correction as f4
import pion_hgcer_method_a_tphi_propagation as f5
import pion_hgcer_method_a_current_lineage_diagnostic as legacy

SCHEMA = "pion_hgcer_method_a_reweighted_contribution_debug/v1"
KINEMATIC = "Q4p4W2p74"
SETTING = "Left-lowe"
SOURCES = ("prompt", "rand", "dummy", "dummy_rand")
SIGNS = {"prompt": 1, "rand": -1, "dummy": -1, "dummy_rand": 1}
F1_HASHES = {
    "Left-lowe": "10086d7d4c42389c9fdd16471c49d59cd189a30980767bcdb87b85914169ef95",
    "Left-highe": "203f4c76f1a251e3e8f231fa3a5e50c9a4fa337420efffffd803205c6d7ea218",
    "Center-lowe": "1593e22b55382b4a9e831d3a1114584e2e3057fcbc4aeeea4aa74c948edf39f8",
    "Center-highe": "5de64b850735ebe70040a703bac999bdd1ac84dc821aba8e1a21ec86ae4b9db3",
    "Right-highe": "77d98006ff81e772c466bf8d10ec83088509a0b8a4d430bded1172bc118bc15f",
}
F3_HASH = "eeb480d9b7ddffab1e97c7f9f27c3099d0780c2bbb4d42a0ee68c077a5752b3d"
F3_MAP = "3a9787fc58d26cc0816012bd1b637ad0c8f201b448d54cb1625a841131154728"
F3_ARTIFACT = "8d2d12068dfa98922d01dfafedbaa4b994bce0bf1e39ce23b1333872313ef121"
F3_ALGORITHM = "ba29630b2f40a87cbadbe751504ce48f23e2b17a08378a2c8a219f93131cb912"
F4_HASH = "1d545924eba89c7f9ffa28028e307aca9b434a89beec06863cf2893887b6b902"
F4_CORRECTION = "bce0cc12ef283b0c361818b436c4e190fc4ea5165436f6c36912cfbb6de53368"
F4_ARTIFACT = "4c935271a0b2723b58b02cc34d82a1e7d5757f7cd36b7707f400c2b28893668a"
F4_HEAD = "b349967c0d4210a78b144ce6134d3c1f15970245"
REFERENCE_HASH = "ecd547de15aa19a50b73bfd8d9dd1c159027c6438ca38b79bc45ff69a8d0e8c7"
REFERENCE_FP = "2c68804d60aaea423cf3dffcd55c3e324c812b97298149df8599ba193ee23616"
REFERENCE_WRAPPER = "7b6c1b244eb15d5d262391ca4fb288d896b9776989c9f658f2877b05dd02205e"
RECORDED_PARENTS = (
    (11019, 0.1081517597541286, 4.114018603391793),
    (19809, 0.3374726818838243, 0.8866714623885792),
    (24552, 0.10645779207740595, 2.0555330576308606),
)
MM_EDGES = tuple(np.linspace(0.8, 1.4, 41).tolist())
SUM_FIELDS = ("B", "B_rw", "A", "A_rw", "positive", "negative", "positive_rw",
              "negative_rw", "sum_b2", "sum_b_rw2", "sum_w0", "sum_w_rw",
              "sum_w02", "sum_w_rw2")
FORBIDDEN = legacy._FORBIDDEN | {
    "event_id", "event_ids", "selected_rows", "top_events", "row", "rows",
    "phi_degrees", "analysis_mm", "weight_tuples", "signed_source_coefficient",
    "signed_baseline_event_contribution", "baseline_pion_weight_w0", "input_paths",
}


class ContributionDebugError(ValueError):
    """A historical authority, alignment or aggregate invariant failed."""


def require(condition, reason):
    if not condition:
        raise ContributionDebugError(reason)


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def fingerprint(value):
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def guard(value):
    if isinstance(value, Mapping):
        for key, child in value.items():
            require(isinstance(key, str) and key.lower() not in FORBIDDEN,
                    "event_level_persistence_forbidden")
            guard(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            guard(child)
    elif isinstance(value, np.ndarray):
        raise ContributionDebugError("transient_array_persistence_forbidden")
    elif isinstance(value, float):
        require(math.isfinite(value), "persisted_nonfinite")
    elif value is not None:
        require(isinstance(value, (str, bool, int)), "persistence_type_invalid")


def seal(value):
    guard(value)
    value = copy.deepcopy(value)
    value["fingerprint"] = fingerprint(value)
    return json.loads(canonical_json(value))


def historical_f4_authority():
    """Independent fixed reference map; never mutates private F.6.3 pins."""
    return {KINEMATIC: {
        "source_file_sha256": F4_HASH, "correction_fingerprint": F4_CORRECTION,
        "artifact_fingerprint": F4_ARTIFACT, "farm_source_head": F4_HEAD,
        "f3_source_file_sha256": F3_HASH, "f3_map_fingerprint": F3_MAP,
        "f3_algorithm_fingerprint": F3_ALGORITHM, "f3_artifact_fingerprint": F3_ARTIFACT,
        "f1_source_file_sha256": dict(F1_HASHES),
    }}


def check_hashes(hashes):
    expected = {**F1_HASHES, "f3": F3_HASH, "f4": F4_HASH}
    require(set(hashes) == set(expected), "historical_input_role_inventory_mismatch")
    for role, digest in expected.items():
        require(hashes[role] == digest, "historical_raw_sha256_mismatch:" + role)


def close(a, b, label):
    require(math.isfinite(a) and math.isfinite(b) and f4._close(a, b, 1e-12), label)


def metric(numerator, denominator, *, signed=False, scale=1.0):
    valid = abs(denominator) > 1e-12 * max(1.0, scale) if signed else denominator > 0
    return {"numerator": numerator, "denominator": denominator,
            "available": valid, "value": numerator / denominator if valid else None,
            "reason": None if valid else "zero_or_near_zero_signed_support" if signed else "zero_support"}


def effective(values):
    require(all(math.isfinite(v) and v >= 0 for v in values), "neff_measure_invalid")
    total = math.fsum(values)
    numerator = total * total
    denominator = math.fsum(v * v for v in values)
    require(math.isfinite(numerator) and math.isfinite(denominator), "squared_moment_nonfinite")
    return {"count": len(values), **metric(numerator, denominator)}


def summary(w, c, b):
    """All records, including zero support. Negative support is signed."""
    rw = [x * y for x, y in zip(w, c)]
    brw = [x * y for x, y in zip(b, c)]
    require(len(w) == len(c) == len(b), "summary_alignment_invalid")
    require(all(math.isfinite(v) for seq in (w, c, b, rw, brw) for v in seq), "contribution_nonfinite")
    require(all(v >= 0 for v in w) and all(v > 0 for v in c), "weight_or_correction_invalid")
    def support(v):
        return math.fsum(v), math.fsum(abs(x) for x in v), math.fsum(x for x in v if x > 0), math.fsum(x for x in v if x < 0), math.fsum(x*x for x in v)
    B, A, pos, neg, sq = support(b)
    R, AR, posr, negr, sqr = support(brw)
    result = {"count": len(w), "zero_weight_count": sum(x == 0 for x in w),
              "B": B, "B_rw": R, "A": A, "A_rw": AR,
              "positive": pos, "negative": neg, "positive_rw": posr, "negative_rw": negr,
              "sum_b2": sq, "sum_b_rw2": sqr,
              "sum_w0": math.fsum(w), "sum_w_rw": math.fsum(rw),
              "sum_w02": math.fsum(x*x for x in w), "sum_w_rw2": math.fsum(x*x for x in rw),
              "signed_delta": R-B, "absolute_delta": AR-A,
              "cancellation": metric(abs(B), A), "cancellation_rw": metric(abs(R), AR),
              "signed_relative_delta": metric(R-B, B, signed=True, scale=A),
              "b_squared_ratio": metric(sqr, sq),
              "w_squared_ratio": metric(math.fsum(x*x for x in rw), math.fsum(x*x for x in w)),
              "neff": {name: effective(v) for name, v in (("w0", w), ("w_rw", rw), ("abs_b", [abs(x) for x in b]), ("abs_b_rw", [abs(x) for x in brw]))}}
    guard(result)
    return result


def distribution(values):
    require(all(math.isfinite(x) and x >= 0 for x in values), "distribution_invalid")
    keys = ("min", "p01", "p50", "p90", "p99", "max")
    quantiles = map(float, np.percentile(values, [0, 1, 50, 90, 99, 100])) if values else [None]*6
    return {"count": len(values), "zero_count": sum(x == 0 for x in values), **dict(zip(keys, quantiles))}


def bin_index(value, edges):
    require(math.isfinite(value), "display_coordinate_nonfinite")
    index = bisect.bisect_right(edges, value)
    return len(edges)-1 if value == edges[-1] else index


def histogram_pair(before, after):
    # A common data-derived display axis only, unrelated to physics selections.
    maximum = max([0.0, *before, *after])
    edges = np.linspace(0, maximum if maximum > 0 else 1, 41).tolist()
    require(all(a < b for a, b in zip(edges, edges[1:])), "display_weight_edges_invalid")
    def counts(values):
        result = [0]*42
        for value in values:
            result[bin_index(value, edges)] += 1
        return result
    return {"edges": edges, "baseline_counts": counts(before), "reweighted_counts": counts(after)}


def tails(rows, w, c, b):
    brw = [x*y for x, y in zip(b, c)]
    totals = summary(w, c, b)
    result = {}
    for name, ranks in (("C", c), ("abs_b_rw", [abs(x) for x in brw]), ("abs_b", [abs(x) for x in b])):
        # Identity is transient only; ties break by the source-owned sorted row order.
        order = sorted(range(len(rows)), key=lambda i: (-ranks[i], rows[i]["source_label"], rows[i]["entry_index"]))
        selected = []
        for q in (.01, .05):
            ids = order[:math.ceil(q*len(rows))]
            subset = summary([w[i] for i in ids], [c[i] for i in ids], [b[i] for i in ids])
            selected.append({"fraction": q, "selected_count": len(ids), "support": subset,
                             "baseline_abs_share": metric(subset["A"], totals["A"]),
                             "reweighted_abs_share": metric(subset["A_rw"], totals["A_rw"]),
                             "baseline_signed_share": metric(subset["B"], totals["B"], signed=True, scale=totals["A"]),
                             "reweighted_signed_share": metric(subset["B_rw"], totals["B_rw"], signed=True, scale=totals["A_rw"]),
                             "sources": {s: summary([w[i] for i in ids if rows[i]["source_label"] == s], [c[i] for i in ids if rows[i]["source_label"] == s], [b[i] for i in ids if rows[i]["source_label"] == s]) for s in SOURCES}})
        result[name] = selected
    return result


def measure_parent(rows, raw_by_id, review, parent, phi_edges):
    require(len(phi_edges) == 10 and tuple(phi_edges) == f5._PHI_EDGES, "phi_geometry_invalid")
    identities = [(r["source_label"], r["entry_index"]) for r in rows]
    require(len(identities) == len(set(identities)), "application_identity_duplicate")
    c = np.asarray(review["correction_factors"], dtype=float)
    r = np.asarray(review["raw_shape_factors"], dtype=float)
    require(c.shape == r.shape == (len(rows),), "review_factor_alignment_invalid")
    w, b, mm, phi = [], [], [], []
    for row in rows:
        identity = (row["source_label"], row["entry_index"])
        require(identity in raw_by_id, "original_mm_identity_missing")
        raw = raw_by_id[identity]
        require(raw["t_index"] == row["t_index"] and raw["phi_index"] == row["phi_index"], "original_mm_assignment_mismatch")
        source = row["source_label"]
        require(source in SIGNS, "source_class_invalid")
        coeff = float(row["signed_source_coefficient"])
        require(math.isfinite(coeff) and coeff*SIGNS[source] > 0, "source_sign_invalid")
        weight, baseline = float(row["baseline_pion_weight_w0"]), float(row["signed_baseline_event_contribution"])
        close(baseline, coeff*weight, "signed_baseline_identity_mismatch")
        w.append(weight); b.append(baseline)
        mm.append(bin_index(float(raw["analysis_MM"]), MM_EDGES))
        pi = bin_index(float(row["phi_degrees"]), phi_edges)-1
        require(pi == row["phi_index"] and 0 <= pi < 9, "phi_assignment_invalid")
        phi.append(pi)
    c = c.tolist()
    totals = summary(w, c, b)
    require(totals["count"] == parent["application_event_count"], "f4_parent_count_mismatch")
    require(totals["B"] == parent["baseline_parent_sum"] and totals["B_rw"] == parent["adjusted_parent_sum"], "f4_parent_exact_sum_mismatch")
    require(math.fsum(x*y for x, y in zip(b, r)) == parent["raw_shape_parent_sum"], "f4_raw_shape_sum_mismatch")
    require(parent["closure_passed"] is True and abs(totals["signed_delta"]) <= parent["closure_relative_scale"], "f4_parent_closure_failed")
    def subset(indices):
        return summary([w[i] for i in indices], [c[i] for i in indices], [b[i] for i in indices])
    groups = [[[] for _ in range(42)] for _ in range(9)]
    for i, (p, m) in enumerate(zip(phi, mm)):
        groups[p][m].append(i)
    joint = []
    for p in range(9):
        cells = []
        for m in range(42):
            ids = groups[p][m]
            cells.append({**subset(ids), "sources": {s: subset([i for i in ids if rows[i]["source_label"] == s]) for s in SOURCES}})
        joint.append(cells)
    source = {s: {"totals": subset([i for i, row in enumerate(rows) if row["source_label"] == s]),
                  "tails": tails([row for row in rows if row["source_label"] == s], [w[i] for i, row in enumerate(rows) if row["source_label"] == s], [c[i] for i, row in enumerate(rows) if row["source_label"] == s], [b[i] for i, row in enumerate(rows) if row["source_label"] == s])} for s in SOURCES}
    phi_marginal = [subset([i for i, pi in enumerate(phi) if pi == p]) for p in range(9)]
    mm_marginal = [subset([i for i, mi in enumerate(mm) if mi == m]) for m in range(42)]
    # Review arrays carry no IDs: the unchanged F.4 sorted row order is the
    # pairing authority. Independently check its persisted source/phi sums.
    f4_sources = {s["source_label"]: s for s in parent["source_diagnostics"]}
    require(len(f4_sources) == len(parent["source_diagnostics"]) and set(f4_sources) == {s for s in SOURCES if source[s]["totals"]["count"]}, "f4_source_inventory_mismatch")
    for name, recorded in f4_sources.items():
        actual = source[name]["totals"]
        require(actual["count"] == recorded["event_count"], "f4_source_count_mismatch")
        for field, old in (("B", "baseline_signed_sum"), ("B_rw", "adjusted_signed_sum"), ("signed_delta", "signed_delta")):
            close(actual[field], recorded[old], "f4_source_sum_mismatch:" + name)
    f4_phi = f5._parent_phi_diagnostics(parent)
    require(set(f4_phi) == {i for i, s in enumerate(phi_marginal) if s["count"]}, "f4_phi_inventory_mismatch")
    for index, recorded in f4_phi.items():
        actual = phi_marginal[index]
        require(actual["count"] == recorded["event_count"], "f4_phi_count_mismatch")
        for field, old in (("B", "baseline_signed_sum"), ("B_rw", "adjusted_signed_sum"), ("signed_delta", "signed_delta")):
            close(actual[field], recorded[old], "f4_phi_sum_mismatch")
    # Independent unchanged October 4 producer semantics, from the SAME inputs.
    previous = legacy._sensitivity_parent(rows, raw_by_id, review, parent, phi_edges)
    compare_legacy(phi_marginal, previous["phi"]["bins"], "historical_phi")
    compare_legacy(mm_marginal, previous["coordinate_bins"]["analysis_MM"]["bins"], "historical_mm")
    brw = [x*y for x, y in zip(b, c)]
    wrw = [x*y for x, y in zip(w, c)]
    return {"canonical_t_index": parent["canonical_t_index"], "totals": totals,
            "f4_parent": {k: parent[k] for k in ("application_event_count", "baseline_parent_sum", "adjusted_parent_sum", "parent_normalization", "closure_relative_scale", "closure_residual")},
            "distributions": {name: distribution(v) for name, v in (("w0", w), ("w_rw", wrw), ("abs_b", [abs(x) for x in b]), ("abs_b_rw", [abs(x) for x in brw]), ("C", c))},
            "histograms": {"weights": histogram_pair(w, wrw), "contributions": histogram_pair([abs(x) for x in b], [abs(x) for x in brw])},
            "tails": tails(rows, w, c, b), "sources": source,
            "phi": phi_marginal, "mm": mm_marginal, "joint": joint,
            "historical_1d_check": "unchanged October 4 producer on exact reproduced lineage"}


def compare_legacy(measured, recorded, label):
    require(len(measured) == len(recorded), label + "_inventory_mismatch")
    for a, b in zip(measured, recorded):
        require(a["count"] == b["count"], label + "_count_mismatch")
        for new, old in (("B", "baseline_signed_sum"), ("B_rw", "adjusted_signed_sum"), ("A", "absolute_baseline_support"), ("signed_delta", "signed_delta")):
            close(a[new], b[old], label + "_sum_mismatch:" + new)


def validate_reference(reference, sha):
    require(sha == REFERENCE_HASH, "historical_reference_raw_sha256_mismatch")
    core = dict(reference); wrapper = core.pop("artifact_fingerprint", None)
    require(wrapper == REFERENCE_WRAPPER and fingerprint(core) == wrapper, "historical_reference_wrapper_mismatch")
    diagnostic = dict(reference["diagnostic"]); fp = diagnostic.pop("fingerprint", None)
    require(fp == REFERENCE_FP and fingerprint(diagnostic) == fp, "historical_reference_fingerprint_mismatch")
    check_hashes({**diagnostic["input_file_hashes"]["f1"], "f3": diagnostic["input_file_hashes"]["f3"], "f4": diagnostic["input_file_hashes"]["f4"]})
    return reference["diagnostic"]


def build_diagnostic(f1_artifacts, f3_artifact, f4_artifact, *, hashes,
                     provenance, reference=None, reference_sha=None):
    """No caller-supplied authority. CLI verifies raw bytes before this gate."""
    check_hashes(hashes)
    before = fingerprint([f1_artifacts, f3_artifact, f4_artifact])
    try:
        _, persisted, authority = f5._validate_f4_artifact(f4_artifact, hashes["f4"], historical_f4_authority())
        parsed = f3._validate_f1_artifacts(f1_artifacts)
        f5._common_geometry(f1_artifacts)
        map_value, _ = f4._validate_f3_artifact(f3_artifact, parsed, F1_HASHES, 1e-12)
        require(map_value["fingerprint"] == F3_MAP and map_value["algorithm_fingerprint"] == F3_ALGORITHM and f3_artifact["artifact_fingerprint"] == F3_ARTIFACT, "historical_f3_fingerprint_mismatch")
        # Recover ONLY from already hash-pinned, intrinsically validated F.4.
        embedded = persisted["f3_runtime_authority"]
        require(embedded["kinematic_token"] == KINEMATIC and embedded["accepted_authority_match"] is True, "historical_f3_embedded_authority_invalid")
        original_f3 = embedded["accepted"]
        checked = f4._f3_runtime_authority(KINEMATIC, F3_HASH, map_value, F3_ARTIFACT, {KINEMATIC: original_f3})
        require(checked == embedded, "historical_f3_embedded_authority_mismatch")
        computed, reviews = f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(
            f1_artifacts, f3_artifact, f1_input_file_hashes=F1_HASHES,
            f3_input_file_sha256=F3_HASH, accepted_f3_runtime_authority_by_kinematic={KINEMATIC: copy.deepcopy(original_f3)})
        require(computed == persisted and computed["fingerprint"] == F4_CORRECTION, "historical_f4_exact_reproduction_mismatch")
        rows_by_parent = f4._raw_application_rows(f1_artifacts, parsed, 1e-12)
        parents = {(p["setting_id"], p["canonical_t_index"]): p for p in persisted["parents"]}
        review_map = {(p["setting_id"], p["canonical_t_index"]): p for p in reviews}
        require(len(reviews) == len(review_map) == 15 and set(review_map) == set(parents) == set(rows_by_parent), "review_parent_inventory_mismatch")
        selected = next(a for a in f1_artifacts if a["setting"]["phi_setting"] == "Left" and a["setting"]["epsilon_filename_token"] == "lowe")
        require(all(a["setting"]["kinematic_token"] == KINEMATIC for a in f1_artifacts), "historical_input_kinematic_mismatch")
        raw = selected["contract"]
        raw_by_id = {(r["source_label"], r["entry_index"]): r for r in raw["application_records"]}
        require(len(raw_by_id) == len(raw["application_records"]), "raw_identity_duplicate")
        measured = []
        accounted = set()
        for t in range(3):
            key = (SETTING, t)
            rows = rows_by_parent[key]
            for row in rows:
                require(row["t_index"] == t, "parent_row_alignment_mismatch")
                identity = (row["source_label"], row["entry_index"])
                require(identity not in accounted, "parent_row_identity_duplicate")
                accounted.add(identity)
            p = measure_parent(rows, raw_by_id, review_map[key], parents[key], raw["phi_edges"])
            count, B, normalization = RECORDED_PARENTS[t]
            require(p["totals"]["count"] == count and p["totals"]["B"] == B and parents[key]["parent_normalization"] == normalization, "october4_recorded_parent_mismatch:t" + str(t))
            measured.append(p)
        require(accounted == set(raw_by_id), "original_rows_not_all_accounted")
        if reference is not None:
            old = validate_reference(reference, reference_sha)
            old_parents = old["measurement_b"]["parents"]
            require([p["canonical_t_index"] for p in old_parents] == [0, 1, 2], "historical_reference_parent_inventory_mismatch")
            for p, previous in zip(measured, old_parents):
                compare_legacy(p["phi"], previous["phi"]["bins"], "reference_phi")
                compare_legacy(p["mm"], previous["coordinate_bins"]["analysis_MM"]["bins"], "reference_mm")
        require(before == fingerprint([f1_artifacts, f3_artifact, f4_artifact]), "production_objects_changed")
    except (KeyError, TypeError, StopIteration, ValueError, OverflowError) as exc:
        raise ContributionDebugError("historical_lineage_failed:" + str(exc)) from exc
    result = seal({
        "schema_version": SCHEMA, "status": "non_authoritative", "non_authoritative": True,
        "scope": {"kinematic": KINEMATIC, "setting_id": SETTING, "parents": [0, 1, 2], "primary_parents": [0, 2], "reference_parent": 1},
        "role": "detached descriptive historical background contributions; no final yields",
        "production_objects_mutated": False, "method_b_numerical_dependency": False,
        "child_renormalization_performed": False, "correction_modified": False,
        "input_file_hashes": dict(hashes), "input_authority": authority,
        "f1_contract_fingerprints": {p["setting_id"]: p["fingerprints"]["fingerprint"] for p in parsed},
        "f3_map_fingerprint": F3_MAP, "f3_algorithm_fingerprint": F3_ALGORITHM,
        "f3_artifact_fingerprint": F3_ARTIFACT,
        "f4_exact_reproduction": {"payload_identical": True, "fingerprint": F4_CORRECTION},
        "historical_evidence": {"sha256": REFERENCE_HASH, "fingerprint": REFERENCE_FP,
                                "wrapper_fingerprint": REFERENCE_WRAPPER, "reference_json_checked": reference is not None},
        "provenance": copy.deepcopy(provenance), "parents": measured,
        "display_geometry": {"phi_edges": list(raw["phi_edges"]), "mm_edges": list(MM_EDGES), "shape": [3, 9, 42],
                             "boundary": "left closed/right open; equality at last edge in last normal bin; explicit flows"},
        "derivation": {"baseline": "b=s*w0", "reweighted": "b_rw=s*w0*C; w_rw=w0*C", "moments": "math.fsum; all finite application records including zeros",
                       "percentiles": "unweighted NumPy linear p01,p50,p90,p99; no product inference",
                       "tail_policy": "ceil(q*N) for q=0.01,0.05; descending rank; ties by source_label,entry_index transient only; zeros retained",
                       "neff": "(sum nonnegative measure)^2 / sum measure^2; descriptive only, not Poisson errors, covariance or yield uncertainty",
                       "signed_ratio_policy": "unavailable if abs(denominator)<=1e-12*max(1,absolute_support)",
                       "reference_policy": "exact raw lineage pins and recorded counts/B/normalization always required; unchanged historical 1D producer cross-check always; supplied reference additionally checked"},
        "checks": {"historical_authority": True, "f4_exact_reproduction": True, "original_rows_accounted": True,
                   "historical_1d_closure": True, "parent_closure": True, "objects_preserved": True},
    })
    validate_diagnostic(result)
    return result


def reconstruct(parts, total, label):
    require(sum(p["count"] for p in parts) == total["count"], label + "_count_mismatch")
    require(sum(p["zero_weight_count"] for p in parts) == total["zero_weight_count"], label + "_zero_count_mismatch")
    for field in SUM_FIELDS + ("signed_delta", "absolute_delta"):
        close(math.fsum(p[field] for p in parts), total[field], label + "_sum_mismatch:" + field)


def validate_summary(s):
    require(type(s["count"]) is int and 0 <= s["zero_weight_count"] <= s["count"], "summary_count_invalid")
    for suffix in ("", "_rw"):
        close(s["B"+suffix], s["positive"+suffix]+s["negative"+suffix], "signed_support_mismatch")
        close(s["A"+suffix], s["positive"+suffix]-s["negative"+suffix], "absolute_support_mismatch")
        require(s["positive"+suffix] >= 0 and s["negative"+suffix] <= 0, "support_sign_invalid")
    close(s["signed_delta"], s["B_rw"]-s["B"], "delta_invalid")
    close(s["absolute_delta"], s["A_rw"]-s["A"], "absolute_delta_invalid")
    for field in ("sum_b2", "sum_b_rw2", "sum_w0", "sum_w_rw", "sum_w02", "sum_w_rw2"):
        require(s[field] >= 0, "nonnegative_moment_invalid")
    for name, total, square in (("w0", "sum_w0", "sum_w02"), ("w_rw", "sum_w_rw", "sum_w_rw2"), ("abs_b", "A", "sum_b2"), ("abs_b_rw", "A_rw", "sum_b_rw2")):
        expected = {"count": s["count"], **metric(s[total]**2, s[square])}
        require(s["neff"][name] == expected, "neff_moment_mismatch")
        if expected["available"]:
            require(expected["value"] <= s["count"] + 1e-9, "neff_exceeds_count")
    for key, expected in (("cancellation", metric(abs(s["B"]), s["A"])), ("cancellation_rw", metric(abs(s["B_rw"]), s["A_rw"])), ("signed_relative_delta", metric(s["signed_delta"], s["B"], signed=True, scale=s["A"])), ("b_squared_ratio", metric(s["sum_b_rw2"], s["sum_b2"])), ("w_squared_ratio", metric(s["sum_w_rw2"], s["sum_w02"]))):
        require(s[key] == expected, "summary_ratio_mismatch:" + key)


def validate_tails(tables, total):
    require(set(tables) == {"C", "abs_b_rw", "abs_b"}, "tail_rank_inventory_invalid")
    for entries in tables.values():
        require(len(entries) == 2, "tail_fraction_inventory_invalid")
        for q, tail in zip((.01, .05), entries):
            require(tail["fraction"] == q and tail["selected_count"] == math.ceil(q*total["count"]) == tail["support"]["count"], "tail_count_invalid")
            validate_summary(tail["support"])
            require(set(tail["sources"]) == set(SOURCES), "tail_source_inventory_invalid")
            for s in tail["sources"].values():
                validate_summary(s)
            reconstruct(list(tail["sources"].values()), tail["support"], "tail_sources")
            for key, numerator, denominator, signed in (("baseline_abs_share", "A", "A", False), ("reweighted_abs_share", "A_rw", "A_rw", False), ("baseline_signed_share", "B", "B", True), ("reweighted_signed_share", "B_rw", "B_rw", True)):
                expected = metric(tail["support"][numerator], total[denominator], signed=signed, scale=total["A_rw"] if denominator == "B_rw" else total["A"])
                require(tail[key] == expected, "tail_share_mismatch")


def validate_diagnostic(value):
    """Independent persisted-only checker; a valid fingerprint alone is insufficient."""
    guard(value)
    core = dict(value); observed = core.pop("fingerprint", None)
    require(observed == fingerprint(core), "diagnostic_fingerprint_mismatch")
    require(value["schema_version"] == SCHEMA and value["status"] == "non_authoritative" and value["non_authoritative"] is True, "diagnostic_schema_or_role_invalid")
    require(value["scope"] == {"kinematic": KINEMATIC, "setting_id": SETTING, "parents": [0, 1, 2], "primary_parents": [0, 2], "reference_parent": 1}, "persisted_scope_invalid")
    for key in ("production_objects_mutated", "method_b_numerical_dependency", "child_renormalization_performed", "correction_modified"):
        require(value[key] is False, "frozen_role_invalid:" + key)
    check_hashes(value["input_file_hashes"])
    a = value["input_authority"]
    expected = historical_f4_authority()[KINEMATIC]
    require(a.get("kinematic_token") == KINEMATIC and a.get("accepted_authority_match") is True, "persisted_authority_invalid")
    require(a.get("accepted") == {k: expected[k] for k in ("source_file_sha256", "correction_fingerprint", "artifact_fingerprint", "farm_source_head")}, "persisted_authority_pins_invalid")
    require(a.get("observed") == {k: expected[k] for k in ("source_file_sha256", "correction_fingerprint", "artifact_fingerprint")}, "persisted_observed_pins_invalid")
    require(value["f3_map_fingerprint"] == F3_MAP and value["f3_algorithm_fingerprint"] == F3_ALGORITHM and value["f3_artifact_fingerprint"] == F3_ARTIFACT, "persisted_f3_pins_invalid")
    require(set(value["f1_contract_fingerprints"]) == set(F1_HASHES), "persisted_f1_fingerprint_inventory_invalid")
    require(all(type(v) is str and len(v) == 64 and all(c in "0123456789abcdef" for c in v) for v in value["f1_contract_fingerprints"].values()), "persisted_f1_fingerprint_invalid")
    evidence = value["historical_evidence"]
    require(evidence["sha256"] == REFERENCE_HASH and evidence["fingerprint"] == REFERENCE_FP and evidence["wrapper_fingerprint"] == REFERENCE_WRAPPER and type(evidence["reference_json_checked"]) is bool, "persisted_reference_pins_invalid")
    require(set(value["checks"]) == {"historical_authority", "f4_exact_reproduction", "original_rows_accounted", "historical_1d_closure", "parent_closure", "objects_preserved"} and all(v is True for v in value["checks"].values()), "persisted_check_receipt_invalid")
    require(value["f4_exact_reproduction"] == {"payload_identical": True, "fingerprint": F4_CORRECTION}, "persisted_f4_reproduction_invalid")
    require(value["display_geometry"]["shape"] == [3, 9, 42] and value["display_geometry"]["mm_edges"] == list(MM_EDGES) and value["display_geometry"]["phi_edges"] == list(f5._PHI_EDGES), "persisted_geometry_invalid")
    require([p["canonical_t_index"] for p in value["parents"]] == [0, 1, 2], "persisted_parent_inventory_invalid")
    for t, p in enumerate(value["parents"]):
        total = p["totals"]
        validate_summary(total)
        count, B, normalization = RECORDED_PARENTS[t]
        require(total["count"] == count and total["B"] == B and p["f4_parent"]["parent_normalization"] == normalization, "persisted_historical_parent_mismatch")
        require(abs(total["signed_delta"]) <= p["f4_parent"]["closure_relative_scale"], "persisted_parent_closure_failed")
        require(total["B_rw"] == p["f4_parent"]["adjusted_parent_sum"] and total["B"] == p["f4_parent"]["baseline_parent_sum"], "persisted_f4_parent_sum_mismatch")
        require(len(p["joint"]) == len(p["phi"]) == 9 and len(p["mm"]) == 42 and all(len(row) == 42 for row in p["joint"]), "joint_shape_invalid")
        cells = [s for row in p["joint"] for s in row]
        for s in cells + p["phi"] + p["mm"]:
            validate_summary(s)
        reconstruct(cells, total, "joint_parent")
        reconstruct(p["phi"], total, "phi_parent")
        reconstruct(p["mm"], total, "mm_parent")
        for i in range(9):
            reconstruct(p["joint"][i], p["phi"][i], "joint_phi")
        for i in range(42):
            reconstruct([row[i] for row in p["joint"]], p["mm"][i], "joint_mm")
        require(set(p["sources"]) == set(SOURCES), "source_inventory_invalid")
        for cell in cells:
            require(set(cell["sources"]) == set(SOURCES), "joint_source_inventory_invalid")
            for s in cell["sources"].values():
                validate_summary(s)
            reconstruct(list(cell["sources"].values()), cell, "joint_sources")
        for name, source in p["sources"].items():
            validate_summary(source["totals"])
            reconstruct([cell["sources"][name] for cell in cells], source["totals"], "joint_source_parent")
            validate_tails(source["tails"], source["totals"])
        reconstruct([s["totals"] for s in p["sources"].values()], total, "source_parent")
        validate_tails(p["tails"], total)
        for distribution_value in p["distributions"].values():
            require(distribution_value["count"] == count and 0 <= distribution_value["zero_count"] <= count, "distribution_count_invalid")
            quantiles = [distribution_value[k] for k in ("min", "p01", "p50", "p90", "p99", "max")]
            require(all(x >= 0 for x in quantiles) and all(a <= b for a, b in zip(quantiles, quantiles[1:])), "distribution_quantiles_invalid")
        for histogram in p["histograms"].values():
            require(len(histogram["edges"]) == 41 and all(a < b for a, b in zip(histogram["edges"], histogram["edges"][1:])), "histogram_edges_invalid")
            for name in ("baseline_counts", "reweighted_counts"):
                require(len(histogram[name]) == 42 and all(type(x) is int and x >= 0 for x in histogram[name]) and sum(histogram[name]) == count, "histogram_count_invalid")
    return True
