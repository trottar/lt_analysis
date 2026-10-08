"""Read-only, standard-library F.1 forensic parser; never a science acceptance gate.

The canonical five-setting validators require scientific dependencies. This
single-setting parser checks the serialized v2 fingerprints and event semantics,
but cannot validate upstream fits, Phase-A ROOT objects or F.4 reconstruction.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "cuts"))
import pion_hgcer_method_a_acceptance_contract as f1  # standard library only

SCHEMA = "e8_4_f1_baseline_reproducibility_diagnostic/v1"
REVIEWED_F1_SHA = "eb6f659da0511f0f6ea867420fda0ee6feb77e6698509ccd31f80a92cc541c07"
REVIEWED_F4_SHA = "79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7"
SETTING = {"kinematic_token": "Q4p4W2p74", "Q2": 4.4, "W": 2.74,
           "epsilon_setting": "low", "epsilon_filename_token": "lowe",
           "phi_setting": "Left", "particle_type": "kaon"}
FEATURES = ("SHMS_delta", "SHMS_xptar", "SHMS_yptar", "P_hgcer_xAtCer", "P_hgcer_yAtCer")
PROVENANCE = ("phase_a_contract_fingerprint", "phase_a_pion_event_population_fingerprint",
              "method_a_fingerprint", "method_a_event_population_fingerprint",
              "part1_config_fingerprint", "coordinate_fingerprint")
CHILD = ("source_label", "entry_index", "t_index", "phi_index", "phi_low", "phi_high", "phi_status")
VARIABLES = ("t_index", "phi_status", "phi_index", "phi_low", "phi_high", "phi_degrees",
             "analysis_MM", "signed_source_coefficient", "baseline_pion_weight_w0",
             "signed_baseline_event_contribution")
SAMPLE_LIMIT = 10
# Same formula/default as F.4 _close / signed_baseline_identity_relative_tolerance.
# This checks persisted b=c*w0 integrity only; comparisons below remain exact.
IDENTITY_TOLERANCE = 1.0e-12


def require(condition, message):
    if not condition:
        raise ValueError(message)


def number(value, label):
    require(type(value) in (int, float), label + ": finite number required")
    require(math.isfinite(value), label + ": finite number required")
    return float(value)


def integer(value, label):
    require(type(value) is int and value >= 0, label + ": nonnegative integer required")
    return value


def pairs(items):
    result = {}
    for key, value in items:
        require(key not in result, "duplicate JSON key: " + key)
        result[key] = value
    return result


def reject_constant(value):
    raise ValueError("nonfinite JSON: " + value)


def finite_tree(value):
    if isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)
    elif type(value) in (int, float):
        number(value, "JSON number")


def read_pinned(path, expected):
    require(isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected), "invalid SHA-256")
    raw = Path(path).read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    require(digest == expected, "raw SHA-256 mismatch")
    obj = json.loads(raw.decode("utf-8"), object_pairs_hook=pairs, parse_constant=reject_constant)
    require(isinstance(obj, dict), "JSON object required")
    finite_tree(obj)
    return obj, digest


def edges(value, label):
    require(isinstance(value, list) and len(value) >= 2, label + ": edges required")
    result = [number(x, label) for x in value]
    require(all(a < b for a, b in zip(result, result[1:])), label + ": unordered edges")
    return result


def identity(row):
    require(isinstance(row, dict), "event object required")
    source = row.get("source_label")
    require(isinstance(source, str) and bool(source.strip()), "source_label missing")
    return source, integer(row.get("entry_index"), "entry_index")


def setting_identity(value):
    # rand_sub persists main.py's filename tokens; scientific fixtures also use
    # numeric Q2/W. Accept only these two explicit representations of this point.
    setting = f1._artifact_setting(value)
    for key, token in (("Q2", "4p4"), ("W", "2p74")):
        actual = setting[key]
        require((type(actual) is str and actual == token) or
                (type(actual) in (int, float) and actual == SETTING[key]), "setting " + key + " invalid")
        setting[key] = SETTING[key]
    return setting


def assignment(row, bounds, prefix, coordinate=None):
    index = integer(row.get(prefix + "_index"), prefix + "_index")
    require(index < len(bounds) - 1, prefix + ": invalid index")
    low, high = number(row.get(prefix + "_low"), prefix), number(row.get(prefix + "_high"), prefix)
    require(low == bounds[index] and high == bounds[index + 1], prefix + ": geometry mismatch")
    if coordinate is not None:
        x = number(row.get(coordinate), coordinate)
        require(low <= x < high or (index == len(bounds) - 2 and x == high), prefix + ": assignment mismatch")
    return index


def check_fingerprints(contract):
    """Exact serialized v2 projection from F.1 producer / F.2 fingerprint audit."""
    for key in PROVENANCE:
        require(isinstance(contract.get(key), str) and bool(contract[key].strip()), "missing " + key)
    require(isinstance(contract.get("host_state"), str) and bool(contract["host_state"]), "host_state missing")
    require(contract.get("source_target_state") == "post_proton_noRF", "source_target_state invalid")
    metadata = contract.get("feature_metadata")
    expected = {"primary_acceptance_features": list(FEATURES),
                "training_population": "prompt_noRF_nommcuts_P_hgcer_npeSum_gt_0",
                "training_low_definition": "0_lt_P_hgcer_npeSum_le_2",
                "training_control_definition": "P_hgcer_npeSum_gt_2",
                "application_population": "authoritative_physical_pion_control_P_hgcer_npeSum_gt_2",
                "parent_coordinate": "canonical_t", "downstream_yield_coordinates": ["canonical_t", "canonical_phi"],
                "phi_is_training_feature": False, "method_b_numerical_dependency": False,
                "probability_map_constructed": False, "weight_adjustment_constructed": False,
                "future_normalization_policy": "future_parent_t_only_no_tphi_child_renormalization",
                "absolute_leakage_probability_claimed": False}
    require(f1._hash(metadata) == f1._hash(expected), "feature_metadata invalid")
    training, application = contract["method_a_training_records"], contract["application_records"]
    child = [{key: row[key] for key in CHILD} for row in application]
    values = {"method_a_training_population_fingerprint": f1._hash(training),
              "application_population_fingerprint": f1._hash(application),
              "acceptance_feature_metadata_fingerprint": f1._hash(metadata),
              "application_child_assignment_projection_fingerprint": f1._hash(child)}
    for key, value in values.items():
        require(contract.get(key) == value, key + " mismatch")
    summary = contract.get("method_a_training_summary")
    require(isinstance(summary, dict) and isinstance(summary.get("by_t_delta"), list), "training summary missing")
    projection = {"schema_version": contract["schema_version"], "fingerprint_schema_version": contract["fingerprint_schema_version"],
                  **{key: contract[key] for key in PROVENANCE}, "host_state": contract["host_state"],
                  "source_target_state": contract["source_target_state"],
                  **{key: contract[key] for key in ("t_edges", "delta_edges", "phi_edges")},
                  **values, "method_a_closure": summary["by_t_delta"], "feature_metadata": metadata}
    require(contract.get("fingerprint_inputs") == projection, "fingerprint_inputs mismatch")
    require(contract.get("fingerprint") == f1._hash(projection), "contract fingerprint mismatch")


def validate_f1(artifact):
    """Forensic single-setting checks, not the canonical five-setting validator."""
    require(isinstance(artifact, dict), "F.1 object required")
    finite_tree(artifact)
    require(artifact.get("schema_version") == f1.METHOD_A_ACCEPTANCE_EVENT_CONTRACT_ARTIFACT_SCHEMA_VERSION, "F.1 wrapper schema invalid")
    require(setting_identity(artifact.get("setting")) == SETTING, "only complete Q4p4W2p74 Left-lowe setting is supported")
    contract = artifact.get("contract")
    require(isinstance(contract, dict), "F.1 contract missing")
    for body in (artifact, contract):
        for key in ("non_authoritative", "production_objects_mutated", "refinement_applied", "production_application_performed", "event_application_performed"):
            require(body.get(key) is (key == "non_authoritative"), key + " invalid")
    require(contract.get("schema_version") == f1.METHOD_A_ACCEPTANCE_EVENT_CONTRACT_SCHEMA_VERSION, "F.1 contract schema invalid")
    require(contract.get("fingerprint_schema_version") == f1.METHOD_A_ACCEPTANCE_EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION, "F.1 fingerprint schema invalid")
    for key, value in {"status": "available", "available": True, "reason": None, "diagnostic_stage": "complete", "method_b_numerical_dependency": False, "future_weight_adjustment_constructed": False}.items():
        require(key in contract and type(contract[key]) is type(value) and contract[key] == value, key + " invalid")
    t_edges, phi_edges = edges(contract.get("t_edges"), "t"), edges(contract.get("phi_edges"), "phi")
    require(len(t_edges) == 4, "three canonical t parents required")
    edges(contract.get("delta_edges"), "delta")
    controls = set()
    for population in ("method_a_training_records", "application_records"):
        rows = contract.get(population)
        require(isinstance(rows, list) and bool(rows), population + " missing")
        seen = set()
        for row in rows:
            key = identity(row)
            require(key not in seen, "duplicate event identity")
            seen.add(key)
            assignment(row, t_edges, "t")
            for feature in FEATURES:
                number(row.get(feature), feature)
            npe = number(row.get("P_hgcer_npeSum"), "npe")
            if population == "method_a_training_records":
                require(key[0] == "prompt" and row.get("nommcuts") is True and npe > 0, "training selection invalid")
                require(row.get("response_class") == ("low" if npe <= 2 else "control"), "training response class invalid")
                if npe > 2:
                    controls.add(key)
            else:
                require(npe > 2, "application not physical pion control")
    rows = contract["application_records"]
    for row in rows:
        require(row["source_label"] != "prompt" or identity(row) in controls, "prompt application missing training control")
        assignment(row, t_edges, "t", "analysis_t")
        number(row.get("analysis_MM"), "analysis_MM")
        phi = number(row.get("phi_degrees"), "phi_degrees")
        if row.get("phi_status") == "inside_phi":
            assignment(row, phi_edges, "phi", "phi_degrees")
        else:
            require(row.get("phi_status") == "outside_phi" and all(key in row and row[key] is None for key in ("phi_index", "phi_low", "phi_high")), "outside phi geometry invalid")
            require(phi < phi_edges[0] or phi > phi_edges[-1], "outside phi assignment invalid")
        c, w, b = [number(row.get(key), key) for key in ("signed_source_coefficient", "baseline_pion_weight_w0", "signed_baseline_event_contribution")]
        require(w >= 0, "negative baseline weight")
        product = c * w
        require(math.isfinite(product) and abs(b - product) <= IDENTITY_TOLERANCE * max(1.0, abs(b), abs(product)), "signed baseline identity mismatch")
    require({row["t_index"] for row in rows if row["phi_status"] == "inside_phi"} == {0, 1, 2}, "eligible canonical t inventory incomplete")
    check_fingerprints(contract)
    return sorted(rows, key=identity)


def distribution(values, include_values=False):
    counts = Counter(values)
    result = {"min": min(values) if values else None, "max": max(values) if values else None,
              "zero_count": counts.get(0, 0), "distinct_count": len(counts)}
    if include_values:
        # Fixed-size value sample; never export a full event/factor array.
        result["distinct_value_sample"] = [{"value": v, "count": counts[v]} for v in sorted(counts)[:SAMPLE_LIMIT]]
        result["sample_truncated"] = len(counts) > SAMPLE_LIMIT
    return result


def summary(rows):
    baseline = [r["signed_baseline_event_contribution"] for r in rows]
    return {"count": len(rows), "signed_baseline_sum": math.fsum(baseline),
            "absolute_baseline_sum": math.fsum(abs(x) for x in baseline),
            "coefficient_distribution": distribution([r["signed_source_coefficient"] for r in rows], True),
            "baseline_weight_distribution": distribution([r["baseline_pion_weight_w0"] for r in rows]),
            "analysis_MM_range": distribution([r["analysis_MM"] for r in rows])}


def inventory(rows):
    eligible = [r for r in rows if r["phi_status"] == "inside_phi"]
    sources = sorted({r["source_label"] for r in rows})
    def grouped(population):
        return [{"canonical_t_index": t, **summary([r for r in population if r["t_index"] == t]),
                 "by_source": [{"source_label": s, **summary([r for r in population if r["t_index"] == t and r["source_label"] == s])} for s in sources]}
                for t in range(3)]
    return {"all_application": {**summary(rows), "by_t": grouped(rows)},
            "f4_consumed_inside_phi": {**summary(eligible), "by_t": grouped(eligible)},
            "canonical_phi_groups": [{"canonical_t_index": t, "phi_status": status, "phi_index": p,
                                      **summary([r for r in rows if (r["t_index"], r["phi_status"], r["phi_index"]) == (t, status, p)])}
                                     for t, status, p in sorted({(r["t_index"], r["phi_status"], r["phi_index"]) for r in rows}, key=lambda x: (x[0], x[1], -1 if x[2] is None else x[2]))]}


def validate_f4(artifact, t_edges):
    """Pinned wrapper/aggregate parser only; no factors or F.5 acceptance."""
    require(isinstance(artifact, dict), "F.4 object required")
    finite_tree(artifact)
    require(artifact.get("schema_version") == "pion_hgcer_method_a_parent_preserving_correction_artifact/v1", "reviewed F.4 wrapper schema invalid")
    body = artifact.get("correction")
    require(isinstance(body, dict) and body.get("schema_version") == "pion_hgcer_method_a_parent_preserving_correction/v1", "reviewed F.4 correction schema invalid")
    require(body.get("fingerprint_schema_version") == "pion_hgcer_method_a_parent_preserving_correction_fingerprint/v1", "F.4 fingerprint schema invalid")
    for obj in (artifact, body):
        for key, expected in {"non_authoritative": True, "production_objects_mutated": False, "production_application_performed": False, "method_b_numerical_dependency": False, "event_correction_persisted": False, "child_renormalization_performed": False, "correction_applied_to_production": False}.items():
            require(obj.get(key) is expected, "F.4 " + key + " invalid")
    require(body.get("status") == "available" and body.get("available") is True and body.get("accepted_basis") == "hgcer3", "F.4 unavailable/basis invalid")
    require(isinstance(body.get("fingerprint_inputs"), dict) and body.get("fingerprint") == f1._hash(body["fingerprint_inputs"]), "F.4 correction fingerprint invalid")
    # Fingerprint inputs bind the complete parent aggregates used here.
    require(body["fingerprint_inputs"].get("parents") == body.get("parents"), "F.4 parent fingerprint mismatch")
    provenance = artifact.get("provenance")
    require(isinstance(provenance, dict) and "input_paths" in provenance, "F.4 provenance missing")
    require(artifact.get("artifact_fingerprint") == f1._hash({"schema_version": artifact["schema_version"], "correction_fingerprint": body["fingerprint"], "input_paths": provenance["input_paths"]}), "F.4 wrapper fingerprint invalid")
    parents = body.get("parents")
    require(isinstance(parents, list) and len(parents) == 15, "F.4 canonical-five parent inventory invalid")
    expected = {(s, t) for s in ("Left-lowe", "Left-highe", "Center-lowe", "Center-highe", "Right-highe") for t in range(3)}
    seen, selected = set(), {}
    for row in parents:
        require(isinstance(row, dict), "F.4 parent object invalid")
        key = (row.get("setting_id"), integer(row.get("canonical_t_index"), "F.4 t index"))
        require(key in expected and key not in seen, "F.4 parent identity invalid")
        seen.add(key)
        setting = setting_identity(row.get("setting"))
        phi, epsilon = key[0].split("-")
        require(setting == {**SETTING, "phi_setting": phi, "epsilon_filename_token": epsilon,
                            "epsilon_setting": "low" if epsilon == "lowe" else "high"}, "F.4 parent setting mismatch")
        require(number(row.get("canonical_t_low"), "F.4 t low") == t_edges[key[1]] and
                number(row.get("canonical_t_high"), "F.4 t high") == t_edges[key[1] + 1], "F.4 parent t geometry mismatch")
        integer(row.get("application_event_count"), "F.4 application count")
        number(row.get("baseline_parent_sum"), "F.4 signed sum")
        require(number(row.get("absolute_baseline_parent_sum"), "F.4 absolute sum") >= 0, "F.4 absolute sum invalid")
        if key[0] == "Left-lowe":
            selected[key[1]] = row
    return selected


def delta(reference, current):
    difference = current - reference
    return {"reference": reference, "current": current, "exact_equal": current == reference,
            "signed_difference": difference, "absolute_difference": abs(difference),
            "relative_difference": difference / abs(reference) if reference else None}


def reference_comparison(reference, current):
    a, b = {identity(r): r for r in reference}, {identity(r): r for r in current}
    common = sorted(a.keys() & b.keys())
    def changes(keys):
        result = {}
        for variable in VARIABLES:
            changed = [k for k in keys if a[k][variable] != b[k][variable]]
            numeric = [b[k][variable] - a[k][variable] for k in changed if type(a[k][variable]) in (int, float) and type(b[k][variable]) in (int, float)]
            result[variable] = {"changed_count": len(changed), "exact_equal": not changed,
                                "first_changed_ids": [list(k) for k in changed[:SAMPLE_LIMIT]],
                                "numeric_delta_count": len(numeric), "absolute_max_delta": max(map(abs, numeric), default=0.0),
                                "signed_delta_sum": math.fsum(numeric), "absolute_delta_sum": math.fsum(map(abs, numeric))}
        return result
    missing, added = sorted(a.keys() - b.keys()), sorted(b.keys() - a.keys())
    return {"missing_count": len(missing), "added_count": len(added), "common_count": len(common),
            "first_missing_ids": [list(k) for k in missing[:SAMPLE_LIMIT]], "first_added_ids": [list(k) for k in added[:SAMPLE_LIMIT]],
            "variables": changes(common), "by_reference_parent_source": [
                {"reference_t_index": t, "source_label": s, "variables": changes([k for k in common if a[k]["t_index"] == t and k[0] == s])}
                for t, s in sorted({(a[k]["t_index"], k[0]) for k in common})],
            "grouping_policy": "matched rows grouped by reference parent; assignment changes reported separately",
            "reference_inventory": inventory(reference)}


def build_comparison(current, current_sha256, reference=None, reference_sha256=None, reviewed_f4=None, reviewed_f4_sha256=None):
    """Pure builder for already pinned JSON objects; validates independently."""
    rows = validate_f1(current)
    require((reference is None) == (reference_sha256 is None), "reference/hash pairing required")
    require((reviewed_f4 is None) == (reviewed_f4_sha256 is None), "F.4/hash pairing required")
    if reviewed_f4 is not None:
        require(reviewed_f4_sha256 == REVIEWED_F4_SHA, "reviewed F.4 authority SHA-256 mismatch")
    for digest in (current_sha256, reference_sha256, reviewed_f4_sha256):
        require(digest is None or re.fullmatch(r"[0-9a-f]{64}", digest), "invalid SHA-256")
    inv = inventory(rows)
    parents = [{"canonical_t_index": t, "application_event_count": p["count"], "baseline_parent_sum": p["signed_baseline_sum"], "absolute_baseline_parent_sum": p["absolute_baseline_sum"]} for t, p in enumerate(inv["f4_consumed_inside_phi"]["by_t"])]
    result = {"schema_version": SCHEMA, "non_authoritative": True, "failed_owner_diagnostic_only": True,
              "production_objects_mutated": False, "runtime_acceptance_granted": False, "setting": dict(SETTING),
              "input_sha256": {"current_f1": current_sha256, "reference_f1": reference_sha256, "reviewed_f4": reviewed_f4_sha256},
              "source_verified_lineage_role": "persisted F.1 b=c*w0; F.4 consumes inside_phi canonical-t rows; no factor reconstruction",
              "parser_checks": {"kind": "forensic parser, not canonical F.1/F.4 validator", "serialized_fingerprints": "SOURCE VERIFIED", "event_identity_geometry": "SOURCE VERIFIED", "signed_baseline_identity_formula": "F.4 _close, default relative tolerance 1e-12; exact comparisons never relaxed"},
              "inventory": inv, "recomputed_baseline_parents": parents,
              "event_level_source_attribution": "NOT VERIFIED", "reference_comparison": None,
              "reference_role": "absent", "reviewed_f4_parent_comparison": None,
              "data_availability": {"reference_f1": reference is not None, "reviewed_f4": reviewed_f4 is not None, "upstream_fit_templates_weight_array_bin_edges": False},
              "evidence_completeness": [
                  {"question": "persisted current coefficients/weights/MM/geometry and eligible sums", "availability": "observable from F.1 now", "classification": "SOURCE VERIFIED"},
                  {"question": "paired persisted event changes", "availability": "requires independently identified reference F.1", "classification": "SOURCE VERIFIED" if reference is not None else "NOT VERIFIED"},
                  {"question": "fit/template numerical drift, original MM weight-bin assignment and root cause", "availability": "requires upstream parent fits/templates/full weights/bin edges and historical inputs", "classification": "NOT VERIFIED"}],
              "scientific_limits": ["No canonical scientific-equivalence decision or runtime acceptance is made.", "Input bytes alone are not a farm/runtime receipt.", "F.1 does not expose upstream fit amplitudes, template objects, full weight arrays or MM bin edges.", "Persisted coefficient/weight/MM differences do not establish upstream numerical cause.", "Missing reviewed raw F.1 cannot be replaced by an alternative archive identity."]}
    if reference is not None:
        require(reference_sha256 != current_sha256, "reference bytes must be distinct")
        result["reference_comparison"] = reference_comparison(validate_f1(reference), rows)
        result["reference_role"] = "reviewed_raw_F1" if reference_sha256 == REVIEWED_F1_SHA else "historical_comparison_not_reviewed_authority"
        if reference_sha256 == REVIEWED_F1_SHA:
            result["event_level_source_attribution"] = "SOURCE VERIFIED for paired persisted fields only; upstream cause NOT VERIFIED"
        else:
            result["scientific_limits"].append("Historical comparison cannot resolve drift relative to missing reviewed F.1.")
    if reviewed_f4 is not None:
        reviewed = validate_f4(reviewed_f4, current["contract"]["t_edges"])
        result["reviewed_f4_parent_comparison"] = [{"canonical_t_index": p["canonical_t_index"], "metrics": {key: delta(reviewed[p["canonical_t_index"]][key], p[key]) for key in ("application_event_count", "baseline_parent_sum", "absolute_baseline_parent_sum")}} for p in parents]
    finite_tree(result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-f1", required=True, type=Path)
    parser.add_argument("--expected-current-sha256", required=True)
    parser.add_argument("--reference-f1", type=Path)
    parser.add_argument("--expected-reference-sha256")
    parser.add_argument("--reviewed-f4", type=Path)
    parser.add_argument("--expected-reviewed-f4-sha256")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    temporary = None
    try:
        require(bool(args.reference_f1) == bool(args.expected_reference_sha256), "reference path/hash pairing required")
        require(bool(args.reviewed_f4) == bool(args.expected_reviewed_f4_sha256), "F.4 path/hash pairing required")
        if args.reviewed_f4 is not None:
            require(args.expected_reviewed_f4_sha256 == REVIEWED_F4_SHA, "reviewed F.4 authority SHA-256 mismatch")
        inputs = [p for p in (args.current_f1, args.reference_f1, args.reviewed_f4) if p is not None]
        require(all(not p.is_symlink() for p in inputs), "input symlinks are refused")
        require(not args.output.is_symlink() and not args.output.exists(), "output must be fresh and not a symlink")
        require(args.output.parent.is_dir(), "output parent directory must exist")
        resolved = [p.resolve(strict=True) for p in inputs]
        require(all(p.is_file() for p in resolved), "inputs must be regular files")
        require(args.output.resolve() not in resolved, "output aliases input")
        for i, path in enumerate(resolved):
            require(all(not path.samefile(other) for other in resolved[:i]), "input paths/inodes must be distinct")
        loaded = [read_pinned(args.current_f1, args.expected_current_sha256)]
        reference = read_pinned(args.reference_f1, args.expected_reference_sha256) if args.reference_f1 else (None, None)
        reviewed = read_pinned(args.reviewed_f4, args.expected_reviewed_f4_sha256) if args.reviewed_f4 else (None, None)
        result = build_comparison(*loaded[0], *reference, *reviewed)
        raw = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")
        # Recheck pinned inputs before publishing; no timestamps/ambient paths.
        for path, expected in zip(inputs, [args.expected_current_sha256] + ([args.expected_reference_sha256] if args.reference_f1 else []) + ([args.expected_reviewed_f4_sha256] if args.reviewed_f4 else [])):
            require(hashlib.sha256(path.read_bytes()).hexdigest() == expected, "input changed during diagnosis")
        with tempfile.NamedTemporaryFile(dir=args.output.parent, prefix=".f1-diagnostic-", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
        # Atomic no-clobber publication: link fails if a concurrent target exists.
        os.link(temporary, args.output)
        return 0
    except (ValueError, KeyError, TypeError, OverflowError, OSError, f1.MethodAAcceptanceContractUnavailable) as exc:
        print("F.1 diagnostic failed: " + str(exc), file=sys.stderr)
        return 2
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
