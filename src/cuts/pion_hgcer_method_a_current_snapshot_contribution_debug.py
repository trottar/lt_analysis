"""Read-only internal snapshot reproduction; never grants external authority.

Only the stored F.3 reproduction record is passed to the frozen F.4 builder,
after independent raw/content binding checks. No historical acceptance pins
or canonical-owner equivalence decision is consumed or replaced here.
"""

from __future__ import annotations

import copy

import pion_hgcer_method_a_reweighted_contribution_debug as aggregate
import pion_hgcer_method_a_parallel_full_procedure as lineage

f3, f4, f5 = aggregate.f3, aggregate.f4, aggregate.f5
require = aggregate.require
fingerprint = aggregate.fingerprint
SCHEMA = "pion_hgcer_method_a_current_snapshot_contribution_debug/v1"
SNAPSHOT_SCHEMA = "pion_hgcer_method_a_current_snapshot/v1"
KINEMATIC = "Q4p4W2p74"
SETTING = "Left-lowe"
SETTINGS = ("Left-lowe", "Left-highe", "Center-lowe", "Center-highe", "Right-highe")
ROLES = (*SETTINGS, "f3", "f4")
ContributionDebugError = aggregate.ContributionDebugError
seal = aggregate.seal
canonical_json = aggregate.canonical_json
SOURCES = aggregate.SOURCES
ROLE = "current-snapshot candidate — detached diagnostic; descriptive background contributions, no final yields"


def digest(value):
    require(type(value) is str and len(value) == 64 and
            all(c in "0123456789abcdef" for c in value), "snapshot_digest_invalid")
    return value


def check_hashes(hashes):
    require(set(hashes) == set(ROLES), "snapshot_input_role_inventory_mismatch")
    for value in hashes.values():
        digest(value)


def validate_f4_wrapper(artifact):
    """Intrinsic checks only; deliberately does not call f5._authority.

    Shared exact F.4 reproduction below checks the entire scientific body,
    including algorithm, fingerprints, links, roles and all fifteen parents.
    """
    require(artifact["schema_version"] == f5._F4_ARTIFACT_SCHEMA,
            "f4_artifact_schema_invalid")
    body = artifact["correction"]
    require(body["schema_version"] == f5._F4_SCHEMA and
            body["fingerprint_schema_version"] == f5._F4_FINGERPRINT_SCHEMA,
            "f4_body_schema_invalid")
    require(body["fingerprint"] == fingerprint(body["fingerprint_inputs"]),
            "f4_correction_fingerprint_mismatch")
    require(artifact["artifact_fingerprint"] == fingerprint({
        "schema_version": f5._F4_ARTIFACT_SCHEMA,
        "correction_fingerprint": body["fingerprint"],
        "input_paths": artifact["provenance"]["input_paths"]}),
        "f4_artifact_fingerprint_mismatch")
    flags = {
        "non_authoritative": True, "relative_map_consumed": True,
        "relative_map_modified": False, "parent_normalization_constructed": True,
        "correction_constructed": True, "correction_applied_to_production": False,
        "event_correction_evaluated_for_detached_diagnostic": True,
        "event_correction_persisted": False, "child_renormalization_performed": False,
        "downstream_template_application_performed": False,
        "absolute_probability_constructed": False, "production_application_performed": False,
        "production_objects_mutated": False, "method_b_numerical_dependency": False,
        "basis_frozen": True, "manual_review_required": True,
    }
    for name, expected in flags.items():
        require(artifact.get(name) is expected, "f4_wrapper_role_invalid:" + name)
    require(body["available"] is True and body["status"] == "available",
            "f4_unavailable")
    require(len(body["parents"]) == 15, "f4_parent_inventory_invalid")
    return body


def validate_snapshot_record(record):
    require(set(record) == {"schema_version", "input_file_hashes", "scientific_fingerprints",
                            "fingerprint"}, "snapshot_record_inventory_invalid")
    require(record["schema_version"] == SNAPSHOT_SCHEMA, "snapshot_schema_invalid")
    check_hashes(record["input_file_hashes"])
    fps = record["scientific_fingerprints"]
    require(set(fps) == set(ROLES), "snapshot_fingerprint_roles_invalid")
    for role, item in fps.items():
        expected = {"stable_content", "contract"} if role in SETTINGS else (
            {"body", "algorithm", "artifact"} if role == "f3" else {"body", "artifact"})
        require(set(item) == expected, "snapshot_fingerprint_inventory_invalid:" + role)
        for value in item.values():
            digest(value)
    require(record["fingerprint"] == fingerprint({k: v for k, v in record.items()
                                                 if k != "fingerprint"}),
            "snapshot_fingerprint_mismatch")


def reproduce_snapshot(f1_artifacts, f3_artifact, f4_artifact, *, hashes, _paired_parents=None):
    """Validate transient pairs; return aggregates after exact, bound reproduction."""
    require(_paired_parents is None or _paired_parents == {}, "paired_selection_invalid")
    check_hashes(hashes)
    before = fingerprint([f1_artifacts, f3_artifact, f4_artifact])
    try:
        parsed = f3._validate_f1_artifacts(f1_artifacts)
        require(all(p["setting"]["kinematic_token"] == KINEMATIC for p in parsed),
                "snapshot_kinematic_invalid")
        # Role assignment is an independent invariant, not merely set membership.
        require([p["setting"]["phi_setting"] + "-" + p["setting"]["epsilon_filename_token"]
                 for p in f1_artifacts] == list(SETTINGS), "snapshot_f1_role_assignment_invalid")
        f5._common_geometry(f1_artifacts)
        f1_hashes = {s: hashes[s] for s in SETTINGS}
        map_value, _ = f4._validate_f3_artifact(f3_artifact, parsed, f1_hashes, 1e-12)
        persisted = validate_f4_wrapper(f4_artifact)
        embedded = persisted["f3_runtime_authority"]
        require(embedded["kinematic_token"] == KINEMATIC,
                "snapshot_f3_embedded_kinematic_invalid")
        # This is stored internal reproduction input, NOT a new authority map.
        stored = embedded["accepted"]
        checked = f4._f3_runtime_authority(KINEMATIC, hashes["f3"], map_value,
                                         f3_artifact["artifact_fingerprint"], {KINEMATIC: stored})
        require(lineage.first_mismatch(checked, embedded) is None,
                "snapshot_f3_embedded_record_mismatch")
        expected_inputs = [{"setting": p["setting"], "setting_id": p["setting_id"],
                            "source_file_sha256": f1_hashes[p["setting_id"]],
                            "stable_f1_content_fingerprint": p["stable_content_fingerprint"],
                            "f1_contract_fingerprint": p["fingerprints"]["fingerprint"]}
                           for p in parsed]
        require(lineage.first_mismatch(expected_inputs, persisted["input_fingerprints"]) is None,
                "snapshot_f4_f1_binding_mismatch")
        for key, value in {"f3_source_file_sha256": hashes["f3"],
                           "f3_map_fingerprint": map_value["fingerprint"],
                           "f3_algorithm_fingerprint": map_value["algorithm_fingerprint"],
                           "f3_artifact_fingerprint": f3_artifact["artifact_fingerprint"]}.items():
            require(persisted[key] == value, "snapshot_f4_link_mismatch:" + key)
        computed, reviews = f4.build_pion_hgcer_method_a_parent_preserving_correction_with_review_data(
            f1_artifacts, f3_artifact, f1_input_file_hashes=f1_hashes,
            f3_input_file_sha256=hashes["f3"],
            accepted_f3_runtime_authority_by_kinematic={KINEMATIC: copy.deepcopy(stored)})
        mismatch = lineage.first_mismatch(computed, persisted)
        require(mismatch is None, "snapshot_f4_exact_reproduction_mismatch:" +
                (mismatch["path"] if mismatch else ""))
        rows = f4._raw_application_rows(f1_artifacts, parsed, 1e-12)
        parents = {(p["setting_id"], p["canonical_t_index"]): p for p in persisted["parents"]}
        review_map = {(p["setting_id"], p["canonical_t_index"]): p for p in reviews}
        inventory = {(s, t) for s in SETTINGS for t in range(3)}
        require(len(parents) == len(reviews) == len(review_map) == 15 and
                set(rows) == set(parents) == set(review_map) == inventory,
                "snapshot_parent_review_inventory_mismatch")
        # Validate every setting's row/factor/source/phi/MM closure, not just plotted Left-lowe.
        measured = {}
        for setting, artifact in zip(SETTINGS, f1_artifacts):
            raw = artifact["contract"]
            raw_by_id = {(r["source_label"], r["entry_index"]): r
                         for r in raw["application_records"]}
            require(len(raw_by_id) == len(raw["application_records"]), "raw_identity_duplicate")
            accounted = set()
            for t in range(3):
                key = (setting, t)
                for row in rows[key]:
                    identity = (row["source_label"], row["entry_index"])
                    require(row["t_index"] == t and identity not in accounted,
                            "snapshot_row_parent_alignment_invalid")
                    accounted.add(identity)
                p = aggregate.measure_parent(rows[key], raw_by_id, review_map[key], parents[key], raw["phi_edges"])
                p.pop("historical_1d_check")
                p["shared_1d_check"] = "unchanged aggregate producer on this exact snapshot"
                p["support"] = {k: parents[key][k] for k in
                                ("in_support_count", "ood_count", "ood_fraction", "ood_final_correction")}
                validate_parent(p)
                if _paired_parents is not None and setting in PAIRED_SETTINGS:
                    _paired_parents.setdefault(setting, []).append(copy.deepcopy(p))
                if setting == SETTING:
                    measured[t] = p
            require(accounted == set(raw_by_id), "snapshot_original_records_not_accounted:" + setting)
        scientific = {p["setting_id"]: {"stable_content": p["stable_content_fingerprint"],
                                       "contract": p["fingerprints"]["fingerprint"]} for p in parsed}
        scientific["f3"] = {"body": map_value["fingerprint"], "algorithm": map_value["algorithm_fingerprint"],
                            "artifact": f3_artifact["artifact_fingerprint"]}
        scientific["f4"] = {"body": persisted["fingerprint"], "artifact": f4_artifact["artifact_fingerprint"]}
        record = seal({"schema_version": SNAPSHOT_SCHEMA, "input_file_hashes": dict(hashes),
                       "scientific_fingerprints": scientific})
        validate_snapshot_record(record)
        require(before == fingerprint([f1_artifacts, f3_artifact, f4_artifact]), "production_objects_changed")
        return record, [measured[t] for t in range(3)]
    except (KeyError, TypeError, ValueError, IndexError, OverflowError) as exc:
        raise ContributionDebugError("snapshot_failed:" + str(exc)) from exc


def build_diagnostic(snapshot, parents, *, provenance):
    validate_snapshot_record(snapshot)
    value = seal({"schema_version": SCHEMA, "status": "non_authoritative", "non_authoritative": True,
        "scope": {"kinematic": KINEMATIC, "setting_id": SETTING, "parents": [0, 1, 2],
                  "primary_parents": [0, 2], "reference_parent": 1},
        "role": ROLE,
        "canonical_owner_equivalence_not_granted": True,
        "canonical_owner_equivalence": "BLOCKED",
        "canonical_owner_equivalence_evaluated": False,
        "internal_snapshot_consistency": "PASS",
        "production_objects_mutated": False, "method_b_numerical_dependency": False,
        "child_renormalization_performed": False, "correction_modified": False,
        "snapshot": snapshot, "provenance": copy.deepcopy(provenance), "parents": parents,
        "display_geometry": {"phi_edges": list(f5._PHI_EDGES), "mm_edges": list(aggregate.MM_EDGES),
                             "shape": [3, 9, 42]},
        "derivation": {"baseline": "b=s*w0", "reweighted": "b_rw=s*w0*C; w_rw=w0*C",
                       "neff": "(sum nonnegative measure)^2/sum measure^2; descriptive, not uncertainty",
                       "tails": "ceil(q*N), q=0.01,0.05; descending C,abs_b,abs_b_rw; ties by transient source/entry",
                       "moments": "math.fsum; finite records including zeros; signed cancellations retained"},
        "checks": {"f4_exact_reproduction": True, "all_15_parent_closures": True,
                   "original_records_accounted": True, "objects_preserved": True}})
    validate_diagnostic(value)
    return value


def validate_diagnostic(value):
    """Persisted aggregate validation only; never consumes F.1 or factors."""
    aggregate.guard(value)
    require(set(value) == {"schema_version", "status", "non_authoritative", "scope", "role",
            "canonical_owner_equivalence_not_granted", "canonical_owner_equivalence",
            "canonical_owner_equivalence_evaluated", "internal_snapshot_consistency",
            "production_objects_mutated", "method_b_numerical_dependency", "child_renormalization_performed",
            "correction_modified", "snapshot", "provenance", "parents", "display_geometry", "derivation",
            "checks", "fingerprint"}, "persisted_diagnostic_inventory_invalid")
    require(value["role"] == ROLE, "persisted_scientific_role_invalid")
    require(value["fingerprint"] == fingerprint({k: v for k, v in value.items() if k != "fingerprint"}),
            "diagnostic_fingerprint_mismatch")
    require(value["schema_version"] == SCHEMA and value["status"] == "non_authoritative" and
            value["non_authoritative"] is True, "diagnostic_schema_or_role_invalid")
    require(value["scope"] == {"kinematic": KINEMATIC, "setting_id": SETTING, "parents": [0, 1, 2],
                               "primary_parents": [0, 2], "reference_parent": 1}, "persisted_scope_invalid")
    require(value["canonical_owner_equivalence_not_granted"] is True and
            value["canonical_owner_equivalence"] == "BLOCKED" and
            value["canonical_owner_equivalence_evaluated"] is False and
            value["internal_snapshot_consistency"] == "PASS", "persisted_acceptance_invalid")
    for key in ("production_objects_mutated", "method_b_numerical_dependency",
                "child_renormalization_performed", "correction_modified"):
        require(value[key] is False, "frozen_role_invalid:" + key)
    validate_snapshot_record(value["snapshot"])
    require(value["display_geometry"] == {"phi_edges": list(f5._PHI_EDGES),
            "mm_edges": list(aggregate.MM_EDGES), "shape": [3, 9, 42]}, "persisted_geometry_invalid")
    require(value["checks"] == {"f4_exact_reproduction": True, "all_15_parent_closures": True,
            "original_records_accounted": True, "objects_preserved": True}, "persisted_checks_invalid")
    require([p["canonical_t_index"] for p in value["parents"]] == [0, 1, 2], "persisted_parent_inventory_invalid")
    provenance = value["provenance"]
    require(set(provenance) == {"git_head", "git_branch", "worktree_dirty", "input_basenames"},
            "persisted_provenance_inventory_invalid")
    require(type(provenance["worktree_dirty"]) is bool and type(provenance["git_head"]) is str and
            len(provenance["git_head"]) == 40 and all(c in "0123456789abcdef" for c in provenance["git_head"]),
            "persisted_git_identity_invalid")
    names = provenance["input_basenames"]
    require(set(names) == set(ROLES) and all(type(n) is str and n and not any(c in n for c in "/\\:")
                                         for n in names.values()), "persisted_basename_invalid")
    for parent in value["parents"]:
        validate_parent(parent)
    return True


def validate_parent(p):
    """Frozen moment/marginal checks, without historical counts or normalization pins."""
    require(set(p) == {"canonical_t_index", "totals", "f4_parent", "distributions", "histograms",
            "tails", "sources", "phi", "mm", "joint", "shared_1d_check", "support"},
            "persisted_parent_fields_invalid")
    require(set(p["distributions"]) == {"w0", "w_rw", "abs_b", "abs_b_rw", "C"} and
            set(p["histograms"]) == {"weights", "contributions"}, "persisted_distribution_inventory_invalid")
    total = p["totals"]
    aggregate.validate_summary(total)
    count = total["count"]
    require(count == p["f4_parent"]["application_event_count"] and count > 0, "persisted_parent_count_invalid")
    require(total["B"] > 0 and p["f4_parent"]["parent_normalization"] > 0, "persisted_normalization_invalid")
    require(p["f4_parent"]["closure_relative_scale"] == 1e-12*max(1.0, abs(total["B"])), "persisted_closure_tolerance_invalid")
    require(p["f4_parent"]["closure_residual"] == total["signed_delta"], "persisted_closure_residual_invalid")
    support = p["support"]
    require(type(support["ood_count"]) is int and type(support["in_support_count"]) is int and
            0 <= support["ood_count"] <= count and support["in_support_count"] + support["ood_count"] == count,
            "persisted_support_count_invalid")
    require(support["ood_fraction"] == support["ood_count"]/count and
            support["ood_final_correction"] == 1/p["f4_parent"]["parent_normalization"], "persisted_support_invalid")
    require(abs(total["signed_delta"]) <= p["f4_parent"]["closure_relative_scale"], "persisted_parent_closure_failed")
    require(total["B_rw"] == p["f4_parent"]["adjusted_parent_sum"] and total["B"] == p["f4_parent"]["baseline_parent_sum"], "persisted_f4_parent_sum_mismatch")
    require(len(p["joint"]) == len(p["phi"]) == 9 and len(p["mm"]) == 42 and all(len(row) == 42 for row in p["joint"]), "joint_shape_invalid")
    cells = [s for row in p["joint"] for s in row]
    for s in cells + p["phi"] + p["mm"]:
        aggregate.validate_summary(s)
    aggregate.reconstruct(cells, total, "joint_parent")
    aggregate.reconstruct(p["phi"], total, "phi_parent")
    aggregate.reconstruct(p["mm"], total, "mm_parent")
    for i in range(9):
        aggregate.reconstruct(p["joint"][i], p["phi"][i], "joint_phi")
    for i in range(42):
        aggregate.reconstruct([row[i] for row in p["joint"]], p["mm"][i], "joint_mm")
    require(set(p["sources"]) == set(SOURCES), "source_inventory_invalid")
    for cell in cells:
        require(set(cell["sources"]) == set(SOURCES), "joint_source_inventory_invalid")
        for s in cell["sources"].values():
            aggregate.validate_summary(s)
        aggregate.reconstruct(list(cell["sources"].values()), cell, "joint_sources")
    for name, source in p["sources"].items():
        aggregate.validate_summary(source["totals"])
        aggregate.reconstruct([cell["sources"][name] for cell in cells], source["totals"], "joint_source_parent")
        aggregate.validate_tails(source["tails"], source["totals"])
    aggregate.reconstruct([s["totals"] for s in p["sources"].values()], total, "source_parent")
    aggregate.validate_tails(p["tails"], total)
    for distribution_value in p["distributions"].values():
        require(distribution_value["count"] == count and 0 <= distribution_value["zero_count"] <= count, "distribution_count_invalid")
        quantiles = [distribution_value[k] for k in ("min", "p01", "p50", "p90", "p99", "max")]
        require(all(x >= 0 for x in quantiles) and all(a <= b for a, b in zip(quantiles, quantiles[1:])), "distribution_quantiles_invalid")
    for histogram in p["histograms"].values():
        require(len(histogram["edges"]) == 41 and all(a < b for a, b in zip(histogram["edges"], histogram["edges"][1:])), "histogram_edges_invalid")
        for name in ("baseline_counts", "reweighted_counts"):
            require(len(histogram[name]) == 42 and all(type(x) is int and x >= 0 for x in histogram[name]) and sum(histogram[name]) == count, "histogram_count_invalid")
    return True


PAIRED_SETTINGS = ("Left-lowe", "Center-lowe")
SETTING_SCHEMA = "pion_hgcer_method_a_current_snapshot_setting_contribution_debug/v1"
PAIRED_SCHEMA = "pion_hgcer_method_a_current_snapshot_paired_lowe_contribution_debug/v1"
PAIRED_ROLE = "paired low-epsilon current-snapshot candidate — detached diagnostic; separate setting-local background contributions, no final yields"


def reproduce_paired_snapshot(f1_artifacts, f3_artifact, f4_artifact, *, hashes):
    """Retain two of the already measured settings from one all-15-parent pass."""
    selected = {}
    snapshot, _ = reproduce_snapshot(f1_artifacts, f3_artifact, f4_artifact,
                                    hashes=hashes, _paired_parents=selected)
    require(tuple(selected) == PAIRED_SETTINGS, "paired_selection_inventory_invalid")
    return snapshot, selected


def build_setting_diagnostic(snapshot, parents, *, setting_id, provenance):
    """New leaf schema with explicit setting/t identity; v1 numerical body unchanged."""
    require(setting_id in PAIRED_SETTINGS, "paired_setting_invalid")
    value = build_diagnostic(snapshot, copy.deepcopy(parents), provenance=provenance)
    value["schema_version"] = SETTING_SCHEMA
    value["scope"]["setting_id"] = setting_id
    for parent in value["parents"]:
        parent["setting_id"] = setting_id
    value = seal({k: v for k, v in value.items() if k != "fingerprint"})
    validate_setting_diagnostic(value, setting_id=setting_id)
    return value


def validate_setting_diagnostic(value, *, setting_id):
    """Check leaf identities, then reuse the unchanged v1 aggregate validator."""
    require(setting_id in PAIRED_SETTINGS, "paired_setting_invalid")
    require(value["schema_version"] == SETTING_SCHEMA and
            value["scope"]["setting_id"] == setting_id, "paired_leaf_scope_invalid")
    require(value["fingerprint"] == fingerprint({k: v for k, v in value.items() if k != "fingerprint"}),
            "paired_leaf_fingerprint_invalid")
    projected = copy.deepcopy(value)
    for parent in projected["parents"]:
        require(parent.pop("setting_id") == setting_id, "paired_parent_setting_invalid")
    projected["schema_version"] = SCHEMA
    projected["scope"]["setting_id"] = SETTING
    projected = seal({k: v for k, v in projected.items() if k != "fingerprint"})
    validate_diagnostic(projected)
    return True


def build_paired_diagnostic(snapshot, parents_by_setting, *, provenance):
    require(tuple(parents_by_setting) == PAIRED_SETTINGS, "paired_parent_settings_invalid")
    leaves = [build_setting_diagnostic(snapshot, parents_by_setting[s], setting_id=s,
                                      provenance=provenance) for s in PAIRED_SETTINGS]
    value = seal({"schema_version": PAIRED_SCHEMA, "status": "non_authoritative",
        "non_authoritative": True, "role": PAIRED_ROLE,
        "scope": {"kinematic": KINEMATIC, "setting_ids": list(PAIRED_SETTINGS)},
        "snapshot": copy.deepcopy(snapshot), "provenance": copy.deepcopy(provenance),
        "settings": leaves, "canonical_owner_equivalence": "BLOCKED",
        "production_objects_mutated": False, "method_b_numerical_dependency": False,
        "cross_setting_sum_or_normalization_performed": False})
    validate_paired_diagnostic(value)
    return value


def validate_paired_diagnostic(value):
    aggregate.guard(value)
    require(set(value) == {"schema_version", "status", "non_authoritative", "role", "scope",
            "snapshot", "provenance", "settings", "canonical_owner_equivalence",
            "production_objects_mutated", "method_b_numerical_dependency",
            "cross_setting_sum_or_normalization_performed", "fingerprint"}, "paired_inventory_invalid")
    require(value["schema_version"] == PAIRED_SCHEMA and value["status"] == "non_authoritative"
            and value["non_authoritative"] is True and value["role"] == PAIRED_ROLE,
            "paired_schema_or_role_invalid")
    require(value["fingerprint"] == fingerprint({k: v for k, v in value.items() if k != "fingerprint"}),
            "paired_fingerprint_invalid")
    require(value["scope"] == {"kinematic": KINEMATIC, "setting_ids": list(PAIRED_SETTINGS)},
            "paired_scope_invalid")
    require(value["canonical_owner_equivalence"] == "BLOCKED" and
            all(value[k] is False for k in ("production_objects_mutated", "method_b_numerical_dependency",
                                           "cross_setting_sum_or_normalization_performed")),
            "paired_frozen_role_invalid")
    validate_snapshot_record(value["snapshot"])
    require(len(value["settings"]) == 2, "paired_leaf_inventory_invalid")
    for setting_id, leaf in zip(PAIRED_SETTINGS, value["settings"]):
        validate_setting_diagnostic(leaf, setting_id=setting_id)
        require(leaf["snapshot"] == value["snapshot"] and leaf["provenance"] == value["provenance"],
                "paired_common_snapshot_or_provenance_invalid")
    require(value["settings"][0]["display_geometry"] == value["settings"][1]["display_geometry"],
            "paired_geometry_mismatch")
    return True
