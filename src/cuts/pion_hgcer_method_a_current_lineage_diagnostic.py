"""Aggregate-only, detached Left/lowe current-F.6.3 diagnostic.

Measurement A describes observed positive HGC response, without truth calibration.
Measurement B describes the accepted F.4 signed measure; N_abs is never applied.
The complete five-setting input set is needed to reproduce the pinned F.4, but
only Left/lowe is measured. No ROOT, production objects, or event rows are written.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Mapping

import numpy as np

import pion_hgcer_method_a_acceptance_map as f3
import pion_hgcer_method_a_parallel_full_procedure as lineage
import pion_hgcer_method_a_parent_preserving_correction as f4
import pion_hgcer_method_a_tphi_propagation as f5

SCHEMA = "pion_hgcer_method_a_current_lineage_diagnostic/v1"
ARTIFACT_SCHEMA = "pion_hgcer_method_a_current_lineage_diagnostic_artifact/v1"
FEATURES = f3.ACCEPTED_FEATURES
KINEMATIC = "Q4p4W2p74"
SETTING = "Left-lowe"
STEM = KINEMATIC + "_Left_lowe_kaon_pion-background_hgcer-method-a-current-lineage-diagnostic"
LIMITATIONS = {
    "operational_kaon_pid_npe_zero": True,
    "direct_hgc_free_true_pion_tag_present_in_consumed_artifacts": False,
    "hgc_free_pion_tag_available_in_consumed_artifacts": False,
    "direct_pion_to_kaon_misid_calibration_performed": False,
    "absolute_misid_probability_constructed": False,
    "proxy_validity_established": False,
    "proxy_validity_evidence_label": "NOT VERIFIED",
}
PID_SEMANTICS = {
    "kaon_pid_category": "P_hgcer_npeSum == 0",
    "pion_tree": "P_hgcer_npeSum > 0",
    "physical_pion_control": "P_hgcer_npeSum > 2",
    "weak_positive": "0 < P_hgcer_npeSum <= 2",
    "conditional_true_pion_zero_response":
        "P(NPE=0 | true pion,x) is pion-to-kaon HGC mis-ID into the kaon-selected sample",
}
# Fixed display bins only. Complete support is retained in under/overflow bins.
DISPLAY_EDGES = {
    "analysis_MM": np.linspace(0.8, 1.4, 41).tolist(),
    "SHMS_delta": np.linspace(-15.0, 25.0, 41).tolist(),
    "P_hgcer_xAtCer": np.linspace(-80.0, 80.0, 41).tolist(),
    "P_hgcer_yAtCer": np.linspace(-80.0, 80.0, 41).tolist(),
}
_FORBIDDEN = {
    "entry_index", "identities", "identity_table", "event_rows", "raw_event_rows",
    "application_records", "method_a_training_records", "feature_rows",
    "correction_factors", "raw_shape_factors", "alternative_factors", "factors",
    "per_event_correction", "per_event_factor", "production_weight", "response_matrix",
}


class CurrentLineageDiagnosticError(ValueError):
    """Required authority or descriptive measurement is unavailable."""


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def fingerprint(value):
    return hashlib.sha256(canonical_json(value).encode("utf-8")).hexdigest()


def guard_aggregate_persistence(value):
    """Reject event-identity/factor payloads recursively, including nested arrays."""
    if isinstance(value, Mapping):
        for key, child in value.items():
            if not isinstance(key, str) or key.lower() in _FORBIDDEN:
                raise CurrentLineageDiagnosticError("event_level_persistence_forbidden")
            guard_aggregate_persistence(child)
    elif isinstance(value, (list, tuple)):
        for child in value:
            guard_aggregate_persistence(child)
    elif isinstance(value, np.ndarray):
        raise CurrentLineageDiagnosticError("transient_array_persistence_forbidden")


def distribution(values):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)):
        raise CurrentLineageDiagnosticError("distribution_nonfinite")
    keys = ("min", "p01", "p10", "p25", "p50", "p75", "p90", "p99", "max")
    result = dict(zip(keys, map(float, np.percentile(values, [0, 1, 10, 25, 50, 75, 90, 99, 100])))) if len(values) else dict.fromkeys(keys)
    return {"count": len(values), **result}


def _positive(value, name):
    if not math.isfinite(value) or value <= 0:
        raise CurrentLineageDiagnosticError(name + "_invalid")
    return value


def signed_measure(b, r):
    """Descriptive signed/absolute sums only; never returns alternative factors."""
    b, r = np.asarray(b, dtype=float), np.asarray(r, dtype=float)
    if b.ndim != 1 or b.shape != r.shape or not len(b) or not np.all(np.isfinite(b)) or not np.all(np.isfinite(r)) or np.any(r <= 0):
        raise CurrentLineageDiagnosticError("signed_measure_input_invalid")
    B = _positive(math.fsum(b), "B_signed")
    A = _positive(math.fsum(abs(b)), "A_abs")
    U = _positive(math.fsum(b * r), "U_signed")
    V = _positive(math.fsum(abs(b) * r), "V_abs")
    ns = _positive(U / B, "N_signed")
    na = _positive(V / A, "N_abs")
    ratio = _positive(ns / na, "normalization_ratio")
    difference = (ns - na) / na
    if not math.isfinite(difference):
        raise CurrentLineageDiagnosticError("normalization_difference_invalid")
    return {"count": len(b), "B_signed": B, "A_abs": A,
            "cancellation_ratio": abs(B) / A, "U_signed": U, "V_abs": V,
            "N_signed": ns, "N_abs": na, "normalization_ratio": ratio,
            "relative_normalization_difference": difference,
            "N_abs_role": "descriptive_positive_support_comparator_only"}


def source_decomposition(rows, b, r, totals):
    result = []
    for source in sorted({row["source_label"] for row in rows}):
        mask = np.asarray([row["source_label"] == source for row in rows])
        sb, sr = b[mask], r[mask]
        B, A = math.fsum(sb), math.fsum(abs(sb))
        result.append({"source_label": source, "count": int(sum(mask)),
                       "B_signed": B, "A_abs": A, "U_signed": math.fsum(sb * sr),
                       "V_abs": math.fsum(abs(sb) * sr),
                       "signed_fraction": B / totals["B_signed"],
                       "absolute_support_fraction": A / totals["A_abs"],
                       "cancellation_ratio": abs(B) / A if A > 0 else None,
                       "cancellation_available": A > 0})
    return result


def _bin_summary(indices, size, b=None, c=None):
    result = []
    for index in range(size):
        mask = indices == index
        row = {"count": int(sum(mask))}
        if b is not None:
            before = math.fsum(b[mask]); after = math.fsum((b * c)[mask])
            row.update(baseline_signed_sum=before, adjusted_signed_sum=after,
                       signed_delta=after - before, absolute_baseline_support=math.fsum(abs(b[mask])))
        result.append(row)
    return result


def binned(values, edges, b=None, c=None):
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise CurrentLineageDiagnosticError("display_coordinate_nonfinite")
    indices = np.searchsorted(edges, values, side="right")
    indices[values == edges[-1]] = len(edges) - 1
    return {"edges": list(edges), "underflow_and_overflow_included": True,
            "bins": _bin_summary(indices, len(edges) + 1, b, c)}


def _features(rows):
    return np.asarray([[row[name] for name in FEATURES] for row in rows], dtype=float).reshape(-1, 3)


def _response(model, rows, training):
    if not rows:
        return np.asarray([]), np.asarray([], dtype=bool)
    return f3.evaluate_relative_response_grid(model, _features(rows), _features(training))


def population_summary(rows, model, training):
    r, inside = _response(model, rows, training)
    return {"count": len(rows), "in_support_count": int(sum(inside)),
            "ood_count": int(len(rows) - sum(inside)),
            "raw_relative_response": distribution(r[inside]),
            "response_summary_policy": "in_support_only; OOD reported separately",
            "features": {name: distribution([row[name] for row in rows])
                         for name in ("P_hgcer_npeSum",) + FEATURES},
            "coordinate_bins": {name: binned([row[name] for row in rows], DISPLAY_EDGES[name])
                                for name in FEATURES}}


def _band(npe):
    if not math.isfinite(npe) or npe <= 0:
        raise CurrentLineageDiagnosticError("training_not_positive_npe")
    return 0 if npe <= 1 else 1 if npe <= 2 else 2 if npe <= 4 else 3


def _proxy_parent(training, application, model):
    classes = {"weak_positive": [r for r in training if 0 < r["P_hgcer_npeSum"] <= 2],
               "control_response": [r for r in training if r["P_hgcer_npeSum"] > 2]}
    bands = [[row for row in training if _band(row["P_hgcer_npeSum"]) == i] for i in range(4)]
    controls = classes["control_response"]
    physical = [row for row in application if row["source_label"] == "prompt"]
    identity = lambda row: (row["source_label"], row["entry_index"])
    control_ids, physical_ids = set(map(identity, controls)), set(map(identity, physical))
    matched = control_ids & physical_ids
    shift = {"matched_count": len(matched), "training_control_only_count": len(control_ids - physical_ids),
             "physical_control_only_count": len(physical_ids - control_ids),
             "training_control": population_summary(controls, model, training),
             "physical_control": population_summary(physical, model, training),
             "matched_training_control": population_summary([r for r in controls if identity(r) in matched], model, training),
             "matched_physical_control": population_summary([r for r in physical if identity(r) in matched], model, training)}
    differences = {}
    for name in FEATURES + ("raw_relative_response",):
        before = (shift["training_control"][name] if name == "raw_relative_response"
                  else shift["training_control"]["features"][name])
        after = (shift["physical_control"][name] if name == "raw_relative_response"
                 else shift["physical_control"]["features"][name])
        differences[name] = {key: after[key] - before[key] if after[key] is not None and before[key] is not None else None
                             for key in before if key != "count"}
    shift["physical_minus_training_control_quantile_differences"] = differences
    return {"classes": {name: population_summary(rows, model, training) for name, rows in classes.items()},
            "positive_npe_bands": [{"definition": name, **population_summary(rows, model, training)}
                                   for name, rows in zip(("0<NPE<=1", "1<NPE<=2", "2<NPE<=4", "NPE>4"), bands)],
            "training_physical_population_shift": shift}


def _sensitivity_parent(rows, raw_rows, review, parent, phi_edges):
    b = np.asarray([row["signed_baseline_event_contribution"] for row in rows], dtype=float)
    r = np.asarray(review["raw_shape_factors"], dtype=float)
    c = np.asarray(review["correction_factors"], dtype=float)
    inside = np.asarray(review["in_support_mask"], dtype=bool)
    if b.shape != r.shape or b.shape != c.shape or b.shape != inside.shape:
        raise CurrentLineageDiagnosticError("review_row_alignment_invalid")
    if not np.all(np.isfinite(c)) or np.any(c <= 0):
        raise CurrentLineageDiagnosticError("accepted_correction_invalid")
    totals = signed_measure(b, r)
    if totals["B_signed"] != parent["baseline_parent_sum"] or totals["U_signed"] != parent["raw_shape_parent_sum"] or totals["N_signed"] != parent["parent_normalization"]:
        raise CurrentLineageDiagnosticError("accepted_parent_summary_mismatch")
    adjusted = math.fsum(b * c)
    if adjusted != parent["adjusted_parent_sum"]:
        raise CurrentLineageDiagnosticError("accepted_adjusted_sum_mismatch")
    phi = np.asarray([row["phi_index"] for row in rows], dtype=int)
    coords = {name: [raw_rows[(row["source_label"], row["entry_index"])][name] for row in rows]
              for name in DISPLAY_EDGES}
    return {"normalization": totals, "accepted_adjusted_sum": adjusted,
            "accepted_closure_residual": parent["closure_residual"],
            "source_decomposition": source_decomposition(rows, b, r, totals),
            "distributions": {"r": distribution(r), "C_signed": distribution(c),
                              "w0": distribution([row["baseline_pion_weight_w0"] for row in rows]),
                              "b": distribution(b), "abs_b": distribution(abs(b))},
            "support": {"in_support_count": int(sum(inside)), "ood_count": int(len(b) - sum(inside)),
                        "ood_fraction": float(1 - sum(inside) / len(b))},
            "phi": {"edges": phi_edges, "bins": _bin_summary(phi, len(phi_edges) - 1, b, c)},
            "coordinate_bins": {name: binned(values, DISPLAY_EDGES[name], b, c) for name, values in coords.items()}}


def build_diagnostic(f1_artifacts, f3_artifact, f4_artifact, *, f1_input_file_hashes,
                     f3_input_file_sha256, f4_input_file_sha256,
                     kinematic=KINEMATIC, setting_id=SETTING):
    """Fail closed before either measurement; authority cannot be supplied by CLI."""
    if kinematic != KINEMATIC or setting_id != SETTING:
        raise CurrentLineageDiagnosticError("diagnostic_scope_unsupported")
    if dict(f1_input_file_hashes) != lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC[KINEMATIC]["f1_source_file_sha256"]:
        raise CurrentLineageDiagnosticError("current_f1_hash_inventory_mismatch")
    # Same shared validation/reproduction path as F.6.3, without a branch fill.
    try:
        _, persisted, authority = f5._validate_f4_artifact(
            f4_artifact, f4_input_file_sha256, lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC)
        computed, reviews = f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(
            f1_artifacts, f3_artifact, f1_input_file_hashes=f1_input_file_hashes,
            f3_input_file_sha256=f3_input_file_sha256,
            accepted_f3_runtime_authority_by_kinematic=lineage.F6_3_CANDIDATE_F3_RECONSTRUCTION_AUTHORITY_BY_KINEMATIC)
        if computed != persisted or computed.get("fingerprint") != persisted.get("fingerprint"):
            raise CurrentLineageDiagnosticError("f4_exact_reproduction_mismatch")
        parsed = f3._validate_f1_artifacts(f1_artifacts)
        _, models = f4._validate_f3_artifact(f3_artifact, parsed, dict(f1_input_file_hashes), 1e-12)
        rows_by_parent = f4._raw_application_rows(f1_artifacts, parsed, 1e-12)
        f5._common_geometry(f1_artifacts)
    except (f4.MethodAParentPreservingCorrectionError, f5.MethodATPhiPropagationError, f3.MethodAAcceptanceMapError) as exc:
        raise CurrentLineageDiagnosticError("current_lineage_authority_failed:" + str(exc)) from exc
    selected = next(item for item in parsed if item["setting_id"] == SETTING)
    if selected["setting"]["kinematic_token"] != KINEMATIC:
        raise CurrentLineageDiagnosticError("input_kinematic_mismatch")
    raw = next(item for item in f1_artifacts if item["setting"] == selected["setting"])["contract"]
    # Source audit: consumed F.1 fields contain detector response, not a truth tag.
    tag_fields = {"true_species", "truth_pdg", "hgc_free_pion_tag", "true_pion"}
    field_inventory = sorted(set().union(*(set(row) for row in raw["method_a_training_records"] + raw["application_records"])))
    if tag_fields.intersection(field_inventory):
        raise CurrentLineageDiagnosticError("unexpected_truth_tag_schema")
    raw_by_id = {(row["source_label"], row["entry_index"]): row for row in raw["application_records"]}
    review_by_t = {int(row["canonical_t_index"]): row for row in reviews if row["setting_id"] == SETTING}
    parent_by_t = {int(row["canonical_t_index"]): row for row in persisted["parents"] if row["setting_id"] == SETTING}
    if set(review_by_t) != {0, 1, 2} or set(parent_by_t) != {0, 1, 2}:
        raise CurrentLineageDiagnosticError("selected_parent_inventory_invalid")
    proxy, sensitivity = [], []
    for index in range(3):
        training = [row for row in selected["training"] if row["t_index"] == index]
        application = [row for row in selected["application"] if row["t_index"] == index]
        proxy.append({"canonical_t_index": index, **_proxy_parent(training, application, models[(SETTING, index)])})
        sensitivity.append({"canonical_t_index": index, **_sensitivity_parent(
            rows_by_parent[(SETTING, index)], raw_by_id, review_by_t[index], parent_by_t[index], raw["phi_edges"])})
    result = {
        "schema_version": SCHEMA, "status": "available", "available": True, "reason": None,
        "non_authoritative": True, "production_objects_mutated": False,
        "method_b_numerical_dependency": False, "alternative_correction_constructed": False,
        "scope": {"kinematic": KINEMATIC, "setting_id": SETTING, "parents": [0, 1, 2], "primary_parent": 1},
        "pid_semantics": dict(PID_SEMANTICS),
        "input_authority": authority,
        "input_file_hashes": {"f1": dict(f1_input_file_hashes),
                              "f3": f3_input_file_sha256, "f4": f4_input_file_sha256},
        "f1_contract_fingerprints": {item["setting_id"]: item["fingerprints"]["fingerprint"] for item in parsed},
        "f3_map_fingerprint": persisted["f3_map_fingerprint"],
        "f3_artifact_fingerprint": persisted["f3_artifact_fingerprint"],
        "f3_reconstruction_role": "candidate_construction_sentinel_only; not farm authority",
        "candidate_validation_source_head": lineage.F6_3_CANDIDATE_VALIDATION_SOURCE_HEAD,
        "f4_exact_reproduction": {"payload_identical": True, "fingerprint": computed["fingerprint"]},
        "measurement_a": {"role": "weak_positive_relative_proxy_evidence_only",
                          "limitations": dict(LIMITATIONS), "audited_f1_record_fields": field_inventory,
                          "tag_audit_scope": "consumed_F1_F3_F4_only; not a global data-availability claim",
                          "parents": proxy},
        "measurement_b": {"role": "accepted_F4_signed_normalization_sensitivity", "parents": sensitivity},
        "display_policy": {"fixed_edges": DISPLAY_EDGES, "phi_edges": raw["phi_edges"],
                           "underflow_overflow": "retained", "clipping": False, "winsorization": False,
                           "adaptive_optimization": False, "new_physics_cuts": False,
                           "percentiles": "unweighted NumPy linear; complete finite population"},
    }
    # Field names are audit metadata only, not event records.
    guard_aggregate_persistence(result)
    result["fingerprint"] = fingerprint(result)
    return json.loads(canonical_json(result))


def build_artifact(diagnostic, *, provenance):
    guard_aggregate_persistence(diagnostic)
    core = dict(diagnostic); observed = core.pop("fingerprint", None)
    if observed != fingerprint(core):
        raise CurrentLineageDiagnosticError("diagnostic_fingerprint_mismatch")
    artifact = {"schema_version": ARTIFACT_SCHEMA, "diagnostic": diagnostic,
                "provenance": provenance, "non_authoritative": True}
    guard_aggregate_persistence(artifact)
    artifact["artifact_fingerprint"] = fingerprint(artifact)
    return json.loads(canonical_json(artifact))
