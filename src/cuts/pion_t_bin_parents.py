"""Particle-subtraction-stage ownership for authoritative per-t pion parents.

The cuts layer supplies the sole signed pion-control cache and direct final
proton-stage kaon objects, with RF restoration applied only when configured.
This module owns only parent fitting, freezing, diagnostics, and PDF
rendering; it never rebuilds parent inputs from averages.
"""

from __future__ import annotations

import textwrap

from background_config import (
    get_particle_subtraction_setting_key,
    resolve_particle_subtraction_mode,
    resolve_pion_subtraction_scope,
)
from pion_component_fits import (
    build_particle_subtraction_component_result,
    load_or_resolve_pion_component_alignment,
    print_particle_subtraction_component_application_pages,
    print_particle_subtraction_component_fit_pages,
    print_particle_subtraction_kaon_lambda_comparison_page,
    record_particle_subtraction_page,
    resolve_scope_component_shapes,
    resolve_scope_single_shape,
)
from pion_component_subtraction import (
    build_simc_shape_pion_control_weights,
    build_t_bin_pion_parent_identity,
    evaluate_particle_subtraction_component_fit_result,
    fingerprint_histogram_content_error,
    resolve_frozen_parent_application_policy,
    validate_authoritative_t_bin_pion_parent,
    validate_frozen_t_bin_pion_parent_collection,
)
from root_histogram_ownership import clone_root_histogram
from data_coordinates import validate_kaon_data_coordinate_contract


def _is_component_t_bin_production(inp_dict):
    return (
        str((inp_dict or {}).get("ParticleType", "")).strip().lower() == "kaon"
        and resolve_particle_subtraction_mode(inp_dict) == "simc_shape_components"
        and resolve_pion_subtraction_scope(inp_dict) == "t_bin"
    )


def _hist_integral(histogram):
    if histogram is None:
        return 0.0
    try:
        return float(histogram.Integral())
    except Exception:
        return 0.0


def _diagnostic_application_from_fit_model(fit_result, *, t_index):
    """Build a detached parent-level model application for diagnostics only.

    The fit remains authoritative even if these optional display clones cannot
    be built.  The payload is deliberately never eligible for child use.
    """
    try:
        before_source = (
            fit_result.get("H_kaon_nosub_input")
            or fit_result.get("H_kaon_input")
        )
        pion_model = fit_result.get("H_kaon_pion_model")
        if before_source is None or pion_model is None:
            raise RuntimeError("missing parent fit-model histogram")

        before = clone_root_histogram(
            before_source,
            scope="pion_parent_t{}".format(int(t_index) + 1),
            role="diagnostic_before_pion",
            name="H_MM_parent_t{}_before_pion".format(int(t_index) + 1),
        )
        template = clone_root_histogram(
            pion_model,
            scope="pion_parent_t{}".format(int(t_index) + 1),
            role="diagnostic_pion_model",
            name="H_MM_parent_t{}_pion_model".format(int(t_index) + 1),
        )
        after = clone_root_histogram(
            before,
            scope="pion_parent_t{}".format(int(t_index) + 1),
            role="diagnostic_after_pion",
            name="H_MM_parent_t{}_after_pion".format(int(t_index) + 1),
        )
        after.Add(template, -1.0)

        weights = build_simc_shape_pion_control_weights(fit_result)
        return {
            "accepted": True,
            "fallback_used": False,
            "fallback_reason": "",
            "diagnostic_only": True,
            "application_authoritative": False,
            "production_applied": False,
            "application_context": "authoritative_parent_t_integrated_diagnostic",
            "analysis_scope": fit_result.get("analysis_scope"),
            "H_pion_control_model": weights.get("H_pion_control_model"),
            "H_kaon_pion_model": weights.get("H_kaon_pion_model") or template,
            "H_pion_weight_vs_MM": weights.get("H_pion_weight_vs_MM"),
            "H_pion_subtraction_template_MM": template,
            "H_pion_subtraction_template_MM_nosub": template,
            "H_MM_before_pion_subtraction": before,
            "H_MM_after_pion_subtraction": after,
            "H_MM_nosub_before_pion_subtraction": before,
            "H_MM_nosub_after_pion_subtraction": after,
            "diagnostics": {
                **dict(weights.get("diagnostics") or {}),
                "event_template_closure": {
                    "status": "parent_fit_model_diagnostic",
                    "child_event_templates_remain_downstream_only": True,
                },
            },
        }, {
            "status": "available",
            "mode": "detached_parent_fit_model",
            "event_template_owned_by": "downstream_child_application",
        }
    except Exception as exc:
        return None, {
            "status": "unavailable",
            "reason": "diagnostic_application_clone_failed",
            "detail": str(exc),
        }


def _classify_diagnostic_exception(exc):
    detail = str(exc)
    normalized = detail.lower()
    if "missing" in normalized and "hist" in normalized:
        reason = "missing_source_histogram"
    elif "template" in normalized and "binning" in normalized:
        reason = "template_binning_mismatch"
    elif "template" in normalized:
        reason = "weight_build_failed"
    elif "tree" in normalized or "event" in normalized:
        reason = "event_template_failed"
    elif "clone" in normalized:
        reason = "histogram_clone_failed"
    else:
        reason = "proposal_build_failed"
    return {"status": "unavailable", "reason": reason, "detail": detail}


def _normalize_parent_diagnostic_outcome(outcome, fit_result, inp_dict):
    """Normalize legacy and proposal/final diagnostic callback results."""
    application_policy = resolve_frozen_parent_application_policy(
        {"fit_result": fit_result}, inp_dict
    )
    production_evaluation = dict(application_policy.get("evaluation") or {})
    if isinstance(outcome, dict) and "proposed_diagnostic_application_result" in outcome:
        proposed = outcome.get("proposed_diagnostic_application_result")
        proposed_status = dict(outcome.get("proposed_diagnostic_application_status") or {})
        final = outcome.get("final_diagnostic_application_result")
        final_status = dict(outcome.get("final_diagnostic_application_status") or {})
        production_evaluation = dict(outcome.get("production_evaluation") or production_evaluation)
    else:
        proposed = outcome if isinstance(outcome, dict) else None
        proposed_status = {
            "status": "available" if proposed is not None else "unavailable",
            "reason": None if proposed is not None else "proposal_unavailable",
            "detail": None,
        }
        final = proposed
        final_status = {
            "status": "available" if final is not None else "unavailable",
            "final_status": "applied_component" if final is not None else "unavailable",
            "final_reason": None if final is not None else "proposal_unavailable",
        }

    proposal_available = isinstance(proposed, dict) and bool(proposed.get("proposal_available", True))
    final_available = isinstance(final, dict)
    production_accepted = bool(production_evaluation.get("accepted"))
    rejection_reasons = [
        reason.strip()
        for reason in str(production_evaluation.get("reason") or "").split(";")
        if reason.strip()
    ]
    status = {
        "status": "available" if final_available else ("partial" if proposal_available else "unavailable"),
        "mode": "detached_parent_event_template",
        "event_template_owned_by": "parent_diagnostic_application",
        "proposal_status": proposed_status.get("status") or ("available" if proposal_available else "unavailable"),
        "proposal_reason": proposed_status.get("reason"),
        "proposal_detail": proposed_status.get("detail"),
        "production_evaluation": "accepted" if production_accepted else "rejected",
        "production_rejection_reasons": rejection_reasons,
        "fallback_mode": production_evaluation.get("fallback_mode"),
        "application_action": application_policy.get("action"),
        "application_policy": application_policy,
        "final_status": final_status.get("final_status") or ("applied_component" if final_available else "unavailable"),
        "final_reason": final_status.get("final_reason"),
        "detail": final_status.get("detail"),
    }
    return proposed, proposed_status, final, final_status, production_evaluation, status


def _parent_summary(parent):
    fit_result = parent.get("fit_result") or {}
    kaon_diagnostics = (fit_result.get("diagnostics") or {}).get("kaon") or {}
    return {
        "t_bin": int(parent["t_bin_index"]) + 1,
        "t_edges": list(parent["t_edges"]),
        "pion_parent_id": parent.get("pion_parent_id"),
        "A_n": fit_result.get("A_n"),
        "A_delta": fit_result.get("A_delta"),
        "A_sidis": fit_result.get("A_sidis"),
        "S_lambda": fit_result.get("S_lambda"),
        "chi2_ndf": fit_result.get("chi2_ndf_kaon"),
        "p_value": fit_result.get("fit_p_value_kaon"),
        "parent_fit_status": fit_result.get("fit_status_kaon"),
        "status": fit_result.get("fit_status_kaon"),
        "protected_fit_variant": kaon_diagnostics.get("fit_variant"),
        "protected_fit_attempted": kaon_diagnostics.get("protected_fit_attempted"),
        "protected_fit_succeeded": kaon_diagnostics.get("protected_fit_succeeded"),
        "parent_pi_delta_applied": kaon_diagnostics.get("pi_delta_applied"),
        "parent_failure_reason": kaon_diagnostics.get("failure_reason"),
        "proposal_status": (parent.get("proposed_diagnostic_application_status") or {}).get("status"),
        "production_gate_status": (parent.get("diagnostic_application_status") or {}).get("production_evaluation"),
        "production_gate_reason": "; ".join((parent.get("diagnostic_application_status") or {}).get("production_rejection_reasons") or []),
        "fallback_mode": (parent.get("diagnostic_application_status") or {}).get("fallback_mode"),
        "final_application_status": (parent.get("diagnostic_application_status") or {}).get("final_status"),
        "diagnostic_application_status": dict(parent.get("diagnostic_application_status") or {}),
    }


def _sum_parent_histograms(parents, getter, *, role, required=False):
    """Return a complete-only canonical parent sum and its audit record.

    A global is meaningful only when every resolved canonical parent supplied
    the requested product.  In particular, a rejected proposal/final state
    must never silently turn into a partial authoritative global curve.
    """
    parents = tuple(parents or ())
    candidates = [getter(parent) for parent in parents]
    missing = [
        int(parent.get("t_bin_index", index))
        for index, (parent, candidate) in enumerate(zip(parents, candidates))
        if candidate is None
    ]
    record = {
        "role": role,
        "expected_parent_count": len(parents),
        "available_parent_count": len(parents) - len(missing),
        "missing_t_bins": missing,
        "source_parent_ids": [parent.get("pion_parent_id") for parent in parents],
        "complete": not missing,
        "histogram": None,
        "bin_error_closure": {
            "checked_bins_including_underflow_overflow": 0,
            "fingerprint": None,
            "passed": False,
        },
    }
    if required and missing:
        raise RuntimeError(
            "authoritative_canonical_t_global_missing_required_input:{}".format(
                ",".join("t{}".format(index + 1) for index in missing)
            )
        )
    if missing or not candidates:
        return record
    total = clone_root_histogram(
        candidates[0],
        scope="pion_canonical_t_global",
        role=role,
        reset=True,
        sumw2=True,
    )
    for candidate in candidates:
        total.Add(candidate)
    record["histogram"] = total
    record["bin_error_closure"] = {
        "checked_bins_including_underflow_overflow": int(total.GetNbinsX()) + 2,
        "fingerprint": fingerprint_histogram_content_error(total),
        "passed": True,
    }
    return record


def _build_authoritative_canonical_t_global(parents):
    """Summarize only sums of frozen canonical-t parent products.

    The setting-wide pion fit is useful for comparison, but it deliberately
    never enters this object.  This distinction keeps a diagnostic fit from
    being mistaken for the production canonical-t result.
    """
    before_record = _sum_parent_histograms(
        parents, lambda parent: parent.get("H_proton_cleaned_final_rf"), role="input", required=True
    )
    estimated_record = _sum_parent_histograms(
        parents,
        lambda parent: (parent.get("proposed_diagnostic_application_result") or {}).get(
            "H_pion_subtraction_template_MM_nosub"
        ),
        role="estimated_contamination",
    )
    proposed_record = _sum_parent_histograms(
        parents,
        lambda parent: (parent.get("proposed_diagnostic_application_result") or {}).get(
            "H_MM_nosub_after_pion_subtraction"
        ),
        role="proposed_clean",
    )
    final_record = _sum_parent_histograms(
        parents,
        lambda parent: (parent.get("final_diagnostic_application_result") or {}).get(
            "H_MM_nosub_after_pion_subtraction"
        ),
        role="final_clean",
    )
    complete = all(record["complete"] for record in (
        before_record, estimated_record, proposed_record, final_record
    ))
    missing_products = {
        record["role"]: list(record["missing_t_bins"])
        for record in (before_record, estimated_record, proposed_record, final_record)
        if not record["complete"]
    }
    before = before_record["histogram"]
    estimated = estimated_record["histogram"]
    proposed_after = proposed_record["histogram"]
    final_after = final_record["histogram"]
    return {
        "label": "AUTHORITATIVE CANONICAL-t GLOBAL",
        "source": "strict_sum_of_frozen_canonical_t_parent_products",
        "t_parent_count": len(parents),
        "source_target_states": sorted({
            str(parent.get("source_target_state") or "unknown")
            for parent in parents
        }),
        "complete": complete,
        "status": "complete" if complete else "incomplete",
        "missing_products": missing_products,
        "sum_records": {
            record["role"]: {
                key: value for key, value in record.items() if key != "histogram"
            }
            for record in (before_record, estimated_record, proposed_record, final_record)
        },
        "H_MM_input": before if complete else None,
        "H_MM_estimated_contamination": estimated if complete else None,
        "H_MM_proposed_clean": proposed_after if complete else None,
        "H_MM_final_clean": final_after if complete else None,
        "integrals": {
            "input": _hist_integral(before),
            "estimated_contamination": _hist_integral(estimated),
            "proposed_clean": _hist_integral(proposed_after),
            "final_clean": _hist_integral(final_after),
        },
    }


def _parent_input_contract(
    parent_input, *, inp_dict, phi_setting, t_index, t_edges, coordinate_contract
):
    """Validate the direct producer handoff without cloning or refitting it."""
    if not isinstance(parent_input, dict):
        raise RuntimeError("authoritative_t_bin_parent_input_missing:t{}".format(t_index + 1))
    if int(parent_input.get("t_index")) != int(t_index):
        raise RuntimeError("authoritative_t_bin_parent_input_index_mismatch:t{}".format(t_index + 1))
    input_edges = [float(value) for value in (parent_input.get("t_edges") or ())]
    if input_edges != list(t_edges):
        raise RuntimeError("authoritative_t_bin_parent_input_edges_mismatch:t{}".format(t_index + 1))
    proton_histogram = parent_input.get("H_proton_cleaned_final_rf")
    pion_control = parent_input.get("H_pion_control")
    if proton_histogram is None or pion_control is None:
        raise RuntimeError("authoritative_t_bin_parent_input_missing_histogram:t{}".format(t_index + 1))
    recomputed_proton_fingerprint = fingerprint_histogram_content_error(proton_histogram)
    producer_fingerprint = str(parent_input.get("proton_final_output_fingerprint") or "")
    if not producer_fingerprint or producer_fingerprint != recomputed_proton_fingerprint:
        raise RuntimeError(
            "authoritative_t_bin_proton_output_fingerprint_mismatch:t{}".format(t_index + 1)
        )
    recomputed_pion_control_fingerprint = fingerprint_histogram_content_error(pion_control)
    declared_pion_control_fingerprint = str(parent_input.get("pion_control_input_fingerprint") or "")
    if not declared_pion_control_fingerprint or declared_pion_control_fingerprint != recomputed_pion_control_fingerprint:
        raise RuntimeError(
            "authoritative_t_bin_pion_control_input_fingerprint_mismatch:t{}".format(t_index + 1)
        )
    expected_epsilon = str(inp_dict.get("EPSSET", "")).strip().lower()
    if str(parent_input.get("source_epsilon", "")).strip().lower() != expected_epsilon:
        raise RuntimeError("authoritative_t_bin_parent_source_epsilon_mismatch:t{}".format(t_index + 1))
    if str(parent_input.get("consumer_epsilon", "")).strip().lower() != expected_epsilon:
        raise RuntimeError("authoritative_t_bin_parent_consumer_epsilon_mismatch:t{}".format(t_index + 1))
    parent_coordinate = validate_kaon_data_coordinate_contract(
        parent_input.get("kaon_data_coordinate"), phi_setting=phi_setting
    )
    if parent_coordinate["coordinate_fingerprint"] != coordinate_contract["coordinate_fingerprint"]:
        raise RuntimeError("authoritative_t_bin_parent_coordinate_fingerprint_mismatch:t{}".format(t_index + 1))
    if str(parent_input.get("coordinate_fingerprint") or "") != parent_coordinate["coordinate_fingerprint"]:
        raise RuntimeError("authoritative_t_bin_parent_coordinate_contract_missing:t{}".format(t_index + 1))
    if str(parent_input.get("pion_control_coordinate_fingerprint") or "") != parent_coordinate["coordinate_fingerprint"]:
        raise RuntimeError("authoritative_t_bin_pion_control_coordinate_mismatch:t{}".format(t_index + 1))
    return {
        "identity": build_t_bin_pion_parent_identity(inp_dict, phi_setting, t_index, t_edges),
        "proton_histogram": proton_histogram,
        "pion_control": pion_control,
        "proton_producer_fingerprint": producer_fingerprint,
        "proton_input_fingerprint": recomputed_proton_fingerprint,
        "pion_control_input_fingerprint": recomputed_pion_control_fingerprint,
        "kaon_data_coordinate": parent_coordinate,
        "coordinate_fingerprint": parent_coordinate["coordinate_fingerprint"],
    }


def _finalize_frozen_parent_collection(hist, inp_dict, parents, t_bins, *, reuse=None):
    """Publish immutable parents and their complete-only global diagnostics."""
    frozen_parents = tuple(parents)
    hist.update({"_pion_t_bin_parent_results": frozen_parents})
    hist["_pion_t_bin_parent_source"] = "particle_subtraction_stage_rf_restored"
    hist["pion_t_parent_collection_frozen"] = True
    hist["pion_t_amplitude_table"] = [_parent_summary(parent) for parent in frozen_parents]
    canonical_t_global = _build_authoritative_canonical_t_global(frozen_parents)
    hist["_pion_authoritative_canonical_t_global"] = canonical_t_global
    hist["pion_authoritative_canonical_t_global_summary"] = {
        key: value for key, value in canonical_t_global.items() if not key.startswith("H_")
    }
    setting_fit = hist.get("_particle_subtraction_component_fit_setting") or {}
    hist["pion_t_parent_diagnostics"] = {
        "phi_setting": hist.get("phi_setting"),
        "epsilon": str(inp_dict.get("EPSSET", "")).strip().lower(),
        "reuse": dict(reuse or {"reused": False}),
        "setting_wide": {
            "scope": "setting_wide_diagnostic",
            "A_n": setting_fit.get("A_n"),
            "A_delta": setting_fit.get("A_delta"),
            "A_sidis": setting_fit.get("A_sidis"),
            "S_lambda": setting_fit.get("S_lambda"),
            "chi2_ndf": setting_fit.get("chi2_ndf_kaon"),
            "p_value": setting_fit.get("fit_p_value_kaon"),
            "status": setting_fit.get("fit_status_kaon"),
            "diagnostic_only": bool(setting_fit.get("diagnostic_only", True)),
            "label": "SETTING-WIDE DIAGNOSTIC FIT - NON-AUTHORITATIVE",
        },
        "authoritative_canonical_t_global": hist["pion_authoritative_canonical_t_global_summary"],
        "parents": list(hist["pion_t_amplitude_table"]),
    }
    validate_frozen_t_bin_pion_parent_collection(hist, inp_dict, t_bins)
    return frozen_parents


def build_setting_t_bin_pion_parents(
    hist,
    inp_dict,
    *,
    t_bins,
    parent_inputs,
    pion_component_shape_payload,
    kaon_signal_shape_payload,
    kaon_sigma0_shape_payload,
    parent_pion_alignment,
    alignment_outpath,
    mm_offset_data,
    coordinate_contract,
    diagnostic_application_builder=None,
):
    """Fit and freeze one parent from direct proton-stage products.

    ``parent_inputs`` is the ownership boundary for component ``t_bin``
    production.  It is intentionally a cuts-layer product, not an averages
    reconstruction from the original kaon trees.
    """
    if not _is_component_t_bin_production(inp_dict):
        return []
    phi_setting = hist.get("phi_setting")
    coordinate_contract = validate_kaon_data_coordinate_contract(
        coordinate_contract, phi_setting=phi_setting
    )
    manifest = hist.setdefault("pion_component_page_manifest", [])
    parent_inputs = tuple(parent_inputs or ())
    expected_count = max(0, len(t_bins) - 1)
    if len(parent_inputs) != expected_count:
        raise RuntimeError(
            "authoritative_t_bin_parent_input_count_mismatch:{}:{}".format(
                len(parent_inputs), expected_count
            )
        )
    if pion_component_shape_payload is None:
        raise RuntimeError("authoritative_t_bin_parents_require_component_shapes")

    contracts = []
    for t_index in range(expected_count):
        t_edges = [float(t_bins[t_index]), float(t_bins[t_index + 1])]
        contracts.append(_parent_input_contract(
            parent_inputs[t_index],
            inp_dict=inp_dict,
            phi_setting=phi_setting,
            t_index=t_index,
            t_edges=t_edges,
            coordinate_contract=coordinate_contract,
        ))

    # Reuse is deliberately in-memory only.  A phi-only adaptation can keep
    # a frozen parent fit, but no change in its physics identity or either
    # producer fingerprint is allowed to cross this boundary.
    existing_parents = hist.get("_pion_t_bin_parent_results")
    if isinstance(existing_parents, tuple) and len(existing_parents) == expected_count:
        reusable = True
        provenance_changes = []
        for t_index, (parent, contract) in enumerate(zip(existing_parents, contracts)):
            t_edges = [float(t_bins[t_index]), float(t_bins[t_index + 1])]
            try:
                validate_authoritative_t_bin_pion_parent(
                    parent, inp_dict, phi_setting, t_index, t_edges
                )
            except RuntimeError:
                reusable = False
                break
            if (
                parent.get("pion_parent_id") != contract["identity"]["pion_parent_id"]
                or parent.get("proton_producer_fingerprint", parent.get("proton_output_fingerprint"))
                != contract["proton_producer_fingerprint"]
                or parent.get("pion_control_input_fingerprint")
                != contract["pion_control_input_fingerprint"]
                or parent.get("coordinate_fingerprint")
                != contract["coordinate_fingerprint"]
            ):
                reusable = False
                break
            previous_provenance = {
                key: parent.get(key)
                for key in (
                    "canonical_interval_pair_id", "canonical_interval_pair_hash",
                    "analysis_runtime_config_hash", "canonical_phi_edges",
                    "requested_num_phi_bins", "actual_num_phi_bins",
                )
            }
            current_provenance = {
                key: contract["identity"].get(key)
                for key in previous_provenance
            }
            if previous_provenance != current_provenance:
                provenance_changes.append({
                    "t_bin_index": t_index,
                    "previous": previous_provenance,
                    "current": current_provenance,
                })
        if reusable:
            return _finalize_frozen_parent_collection(
                hist,
                inp_dict,
                existing_parents,
                t_bins,
                reuse={
                    "reused": True,
                    "reason": "matching_parent_physics_and_producer_fingerprints",
                    "provenance_only_changes": provenance_changes,
                },
            )

    parents = []
    for t_index in range(len(t_bins) - 1):
        parent_input = parent_inputs[t_index]
        t_edges = [float(t_bins[t_index]), float(t_bins[t_index + 1])]
        contract = contracts[t_index]
        identity = contract["identity"]
        kaon_input = contract["proton_histogram"]
        pion_control = contract["pion_control"]
        proton_fingerprint = contract["proton_input_fingerprint"]
        pion_input = kaon_input
        pion_fingerprint = contract["proton_input_fingerprint"]

        scope = "t_bin{}".format(t_index + 1)
        scope_component_shapes = resolve_scope_component_shapes(
            pion_component_shape_payload,
            analysis_scope=scope,
            t_bin_index=t_index,
        )
        alignment_bin_key = {
            "kinematic_setting": get_particle_subtraction_setting_key(inp_dict),
            "epsilon": str(inp_dict.get("EPSSET", "")).strip().lower(),
            "phi_setting": phi_setting,
            "analysis_scope": "authoritative_particle_stage_t_bin",
            "kaon_data_coordinate_fingerprint": coordinate_contract[
                "coordinate_fingerprint"
            ],
            "t_bin": {"index": int(t_index), "edges": t_edges},
            "phi_bin": None,
            "active_dimensions": {
                "particle_type": inp_dict.get("ParticleType"),
                "polarization": inp_dict.get("POL"),
                "target": inp_dict.get("target") or inp_dict.get("Target"),
            },
        }
        active_alignment, persistence_status, persistence_reasons, alignment_paths = (
            load_or_resolve_pion_component_alignment(
                alignment_outpath,
                get_particle_subtraction_setting_key(inp_dict),
                phi_setting,
                inp_dict.get("EPSSET"),
                "authoritative_particle_stage_t_bin",
                alignment_bin_key,
                pion_control,
                scope_component_shapes,
                parent_alignment=parent_pion_alignment,
                inp_dict=inp_dict,
                common_setting_shift_gev=float(coordinate_contract["mm_shift"]),
            )
        )
        fit_result = build_particle_subtraction_component_result(
            pion_control,
            pion_input,
            scope_component_shapes,
            inp_dict,
            analysis_scope=scope,
            kaon_signal_shape=resolve_scope_single_shape(
                kaon_signal_shape_payload, analysis_scope=scope, t_bin_index=t_index
            ),
            kaon_sigma0_shape=resolve_scope_single_shape(
                kaon_sigma0_shape_payload, analysis_scope=scope, t_bin_index=t_index
            ),
            kaon_sigma0_source_diagnostics=(
                (kaon_sigma0_shape_payload or {}).get("diagnostics") or {}
            ),
            mm_offset_data=float(mm_offset_data or 0.0),
            phi_setting=phi_setting,
            context="particle_stage_{}_t{}".format(phi_setting, t_index + 1),
            parent_alignment=parent_pion_alignment,
            pion_component_alignment=active_alignment,
            alignment_bin_key=alignment_bin_key,
        )
        alignment_payload = fit_result.get("pion_component_alignment")
        if isinstance(alignment_payload, dict):
            alignment_payload["persistence_status"] = persistence_status
            alignment_payload["persistence_rejection_reasons"] = list(persistence_reasons)
            alignment_payload["artifact_paths"] = list(alignment_paths)
            for path in alignment_paths:
                if path not in inp_dict.setdefault("pion_component_alignment_artifacts", []):
                    inp_dict["pion_component_alignment_artifacts"].append(path)
        parent = {
            **identity,
            "t_bin_index": int(t_index),
            "t_edges": t_edges,
            "analysis_scope": scope,
            "input_selection": str(
                parent_input.get("input_selection")
                or "no_rf_proton_cleaning_then_rf_restored"
            ),
            "source_target_state": str(
                parent_input.get("source_target_state") or "post_proton_post_rf"
            ),
            "fit_result": fit_result,
            "H_random_dummy_subtracted_kaon": parent_input.get("H_random_dummy_subtracted_kaon"),
            "H_proton_estimate": parent_input.get("H_proton_estimate"),
            "H_proton_cleaned": parent_input.get("H_proton_cleaned"),
            "H_proton_cleaned_final_rf": kaon_input,
            "H_pion_control": pion_control,
            "source_accounting": dict(parent_input.get("source_accounting") or {}),
            "pion_control_source_accounting": dict(
                parent_input.get("pion_control_source_accounting") or {}
            ),
            "proton_producer_fingerprint": contract["proton_producer_fingerprint"],
            "proton_output_fingerprint": proton_fingerprint,
            "pion_input_fingerprint": pion_fingerprint,
            "pion_control_input_fingerprint": contract["pion_control_input_fingerprint"],
            "kaon_data_coordinate": dict(contract["kaon_data_coordinate"]),
            "coordinate_fingerprint": contract["coordinate_fingerprint"],
            "handoff_match": bool(
                contract["proton_producer_fingerprint"] == proton_fingerprint
            ),
            "handoff_status": "matched_producer_fingerprint",
            "pion_fit_input_is_proton_output": True,
            "application_result": None,
            "fit_summary": None,
            "application_summary": None,
        }
        # Missing K-Lambda is discovered by the existing fitter and remains a
        # hard error.  Optional detached diagnostic products are independent.
        validate_authoritative_t_bin_pion_parent(
            parent, inp_dict, phi_setting, t_index, t_edges
        )
        if diagnostic_application_builder is None:
            diagnostic_payload, diagnostic_status = _diagnostic_application_from_fit_model(
                fit_result, t_index=t_index
            )
            outcome = diagnostic_payload
        else:
            try:
                outcome = diagnostic_application_builder(
                    fit_result=fit_result,
                    parent_input=parent_input,
                    t_index=t_index,
                    t_edges=t_edges,
                )
            except Exception as exc:
                outcome = None
                diagnostic_status = _classify_diagnostic_exception(exc)
                if (
                    bool(inp_dict.get("pion_parent_diagnostic_strict", False))
                    and diagnostic_status["reason"] == "proposal_build_failed"
                ):
                    raise
                production_evaluation = evaluate_particle_subtraction_component_fit_result(
                    fit_result,
                    inp_dict,
                )
                parent["proposed_diagnostic_application_result"] = None
                parent["proposed_diagnostic_application_status"] = dict(diagnostic_status)
                parent["final_diagnostic_application_result"] = None
                parent["final_diagnostic_application_status"] = {
                    "status": "unavailable",
                    "final_status": "unavailable",
                    "final_reason": diagnostic_status["reason"],
                    "detail": diagnostic_status["detail"],
                }
                parent["diagnostic_application_result"] = None
                parent["diagnostic_application_status"] = {
                    **diagnostic_status,
                    "proposal_status": "unavailable",
                    "production_evaluation": "accepted" if production_evaluation.get("accepted") else "rejected",
                    "production_rejection_reasons": [
                        reason.strip()
                        for reason in str(production_evaluation.get("reason") or "").split(";")
                        if reason.strip()
                    ],
                    "fallback_mode": production_evaluation.get("fallback_mode"),
                    "final_status": "unavailable",
                    "final_reason": diagnostic_status["reason"],
                }
                parents.append(parent)
                continue
        if diagnostic_application_builder is None:
            # Retain the legacy model-only fallback for callers that do not
            # provide the event-level parent builder.
            proposed = diagnostic_payload
            proposed_status = diagnostic_status
            final = diagnostic_payload
            final_status = {
                "status": diagnostic_status.get("status"),
                "final_status": "applied_component" if diagnostic_payload is not None else "unavailable",
                "final_reason": diagnostic_status.get("reason"),
                "detail": diagnostic_status.get("detail"),
            }
            production_evaluation = evaluate_particle_subtraction_component_fit_result(fit_result, inp_dict)
            diagnostic_status = {
                "status": diagnostic_status.get("status"),
                "mode": diagnostic_status.get("mode", "detached_parent_fit_model"),
                "proposal_status": diagnostic_status.get("status"),
                "proposal_reason": diagnostic_status.get("reason"),
                "proposal_detail": diagnostic_status.get("detail"),
                "production_evaluation": "accepted" if production_evaluation.get("accepted") else "rejected",
                "production_rejection_reasons": [
                    reason.strip()
                    for reason in str(production_evaluation.get("reason") or "").split(";")
                    if reason.strip()
                ],
                "fallback_mode": production_evaluation.get("fallback_mode"),
                "final_status": final_status["final_status"],
                "final_reason": final_status["final_reason"],
                "detail": final_status["detail"],
            }
        else:
            proposed, proposed_status, final, final_status, production_evaluation, diagnostic_status = (
                _normalize_parent_diagnostic_outcome(outcome, fit_result, inp_dict)
            )
        parent["proposed_diagnostic_application_result"] = proposed
        parent["proposed_diagnostic_application_status"] = proposed_status
        parent["final_diagnostic_application_result"] = final
        parent["final_diagnostic_application_status"] = final_status
        parent["production_application_policy"] = {
            "production_evaluation": "accepted" if production_evaluation.get("accepted") else "rejected",
            "fallback_mode": production_evaluation.get("fallback_mode"),
            "rejection_reasons": list(diagnostic_status.get("production_rejection_reasons") or []),
            "application_action": resolve_frozen_parent_application_policy(parent, inp_dict).get("action"),
        }
        parent["diagnostic_application_result"] = final
        parent["diagnostic_application_status"] = diagnostic_status
        parents.append(parent)

    # A tuple makes the authoritative parent *collection* immutable to normal
    # consumers.  The ROOT state inside each validated parent remains owned by
    # that parent for the lifetime of the setting, but Step 6/yield/average
    # code cannot append, remove, or replace parent slots.
    frozen_parents = tuple(parents)
    # Keep the producer-owned freeze visible at the particle-stage boundary;
    # the shared finalizer only derives scalar summaries and global products.
    hist["_pion_t_bin_parent_results"] = frozen_parents
    frozen_parents = _finalize_frozen_parent_collection(
        hist,
        inp_dict,
        frozen_parents,
        t_bins,
        reuse={"reused": False, "reason": "new_particle_stage_parent_fit"},
    )
    source_states = sorted({
        str(parent.get("source_target_state") or "unknown")
        for parent in parents
    })
    print("[SIMC PION PARENTS] {} {} -- {}".format(
        phi_setting, inp_dict.get("EPSSET", ""), ", ".join(source_states)
    ))
    print("  canonical t bins: {}".format(len(frozen_parents)))
    for row in hist["pion_t_amplitude_table"]:
        print(
            "[SIMC PION PARENT] {setting} {epsilon} t{t_bin} {t_edges} "
            "fit={status} protected={protected} proposal={proposal} "
            "production_gate={gate} final={final} gate_reason={reason} parent={parent}".format(
                setting=phi_setting,
                epsilon=inp_dict.get("EPSSET", ""),
                protected=row.get("protected_fit_variant"),
                proposal=row.get("proposal_status"),
                gate=row.get("production_gate_status"),
                final=row.get("final_application_status"),
                reason=(row.get("production_gate_reason") or "none")[:240],
                parent=str(row.get("pion_parent_id"))[-16:],
                **row
            )
        )
    print("  parents frozen before rand_sub PDF close")
    return frozen_parents


def _record_parent_page_if_absent(manifest, page_id, *, scope):
    if any(entry.get("page_id") == page_id for entry in manifest if isinstance(entry, dict)):
        return
    record_particle_subtraction_page(
        manifest, page_id, scope=scope, authoritative=True
    )


def _print_parent_overview_page(pdf_name, parent, title_prefix, manifest, page_prefix):
    """Emit a compact proton-analog parent overview without owning data hists."""
    try:
        import ROOT

        canvas = ROOT.TCanvas(
            "pion_parent_overview_{}".format(parent.get("pion_parent_id", "unknown")[-12:]),
            "Pion parent overview",
            900,
            850,
        )
        text = ROOT.TPaveText(0.07, 0.05, 0.93, 0.95, "NDC")
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.SetTextAlign(12)
        text.SetTextSize(0.024)
        fit_result = parent.get("fit_result") or {}
        status = parent.get("diagnostic_application_status") or {}
        protected = ((fit_result.get("diagnostics") or {}).get("kaon") or {}).get(
            "pi_delta_signal_protected_fit"
        ) or {}
        lambda_retention = protected.get("lambda_retention") or {}
        text.AddText("{} AUTHORITATIVE PION PARENT".format(title_prefix))
        text.AddText("setting/epsilon: {} / {}".format(
            parent.get("phi_setting"), parent.get("epsilon")
        ))
        text.AddText("t{} [{:.3f}, {:.3f}]".format(
            int(parent["t_bin_index"]) + 1,
            float(parent["t_edges"][0]),
            float(parent["t_edges"][1]),
        ))
        text.AddText("input: {}".format(
            "post-proton, RF-restored"
            if parent.get("source_target_state") == "post_proton_post_rf"
            else "post-proton, RF not restored"
        ))
        text.AddText("source target state: {}".format(parent.get("source_target_state")))
        text.AddText("parent: {}".format(parent.get("pion_parent_id")))
        coordinate = parent.get("kaon_data_coordinate") or {}
        text.AddText("kaon coordinate: dMM={:+.6g}, dt={:+.6g}, fp={}".format(
            float(coordinate.get("mm_shift", 0.0)),
            float(coordinate.get("t_shift", 0.0)),
            str(parent.get("coordinate_fingerprint") or "")[-16:],
        ))
        text.AddText("fit status pion/kaon: {} / {}".format(
            fit_result.get("fit_status_pion"), fit_result.get("fit_status_kaon")
        ))
        text.AddText("production gate: {}  fallback: {}".format(
            status.get("production_evaluation"), status.get("fallback_mode")
        ))
        text.AddText("A_n={}  A_delta proposed={}  A_delta applied={}  A_sidis={}  S_lambda={}".format(
            fit_result.get("A_n"), fit_result.get("A_delta"),
            ((protected.get("protected_applied_amplitudes") or {}).get("pi_delta", fit_result.get("A_delta"))),
            fit_result.get("A_sidis"), fit_result.get("S_lambda"),
        ))
        text.AddText("protected: {}  chi2/ndf={}  p={}".format(
            protected.get("fit_variant") or protected.get("status") or "unavailable",
            protected.get("chi2_ndf"), protected.get("p_value"),
        ))
        text.AddText("Lambda retention: {}".format(
            lambda_retention.get("status") or lambda_retention.get("reason") or "unavailable"
        ))
        text.AddText("proposal: {}  final: {}".format(
            status.get("proposal_status"), status.get("final_status")
        ))
        text.AddText("diagnostic status: {}".format(status.get("status")))
        for reason in status.get("production_rejection_reasons") or []:
            for line in textwrap.wrap("gate reason: {}".format(reason), width=95):
                text.AddText(line)
        for line in textwrap.wrap("diagnostic reason: {}".format(
            status.get("final_reason") or status.get("proposal_reason") or "none"
        ), width=95):
            text.AddText(line)
        for line in textwrap.wrap("detail: {}".format(status.get("detail") or status.get("proposal_detail") or "none"), width=95):
            text.AddText(line)
        text.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(
            manifest, "{}.overview".format(page_prefix), scope=parent["analysis_scope"]
        )
        return True
    except Exception as exc:
        parent["diagnostic_application_status"] = {
            **dict(parent.get("diagnostic_application_status") or {}),
            "overview_status": "unavailable",
            "overview_reason": str(exc),
        }
        return False


def _print_parent_application_gate_page(pdf_name, parent, title_prefix, manifest, page_prefix):
    """Print the production-gate decision without conflating it with proposal availability."""
    try:
        import ROOT

        fit_result = parent.get("fit_result") or {}
        status = parent.get("diagnostic_application_status") or {}
        canvas = ROOT.TCanvas(
            "pion_parent_gate_{}".format(parent.get("pion_parent_id", "unknown")[-12:]),
            "Pion parent application gate",
            900,
            650,
        )
        text = ROOT.TPaveText(0.08, 0.08, 0.92, 0.92, "NDC")
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.SetTextAlign(12)
        text.SetTextSize(0.026)
        text.AddText("{} application gate".format(title_prefix))
        text.AddText("component production evaluation: {}".format(status.get("production_evaluation")))
        text.AddText("fit status pion/kaon: {} / {}".format(
            fit_result.get("fit_status_pion"), fit_result.get("fit_status_kaon")
        ))
        text.AddText("proposal: {}  final: {}  fallback: {}".format(
            status.get("proposal_status"), status.get("final_status"), status.get("fallback_mode")
        ))
        text.AddText("pion validation: {}  kaon validation: {}".format(
            (fit_result.get("diagnostics") or {}).get("pion", {}).get("validation", {}).get("accepted"),
            (fit_result.get("diagnostics") or {}).get("kaon", {}).get("validation", {}).get("accepted"),
        ))
        reasons = status.get("production_rejection_reasons") or ["none"]
        for reason in reasons:
            for line in textwrap.wrap("reason: {}".format(reason), width=100):
                text.AddText(line)
        text.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(
            manifest, "{}.application_gate".format(page_prefix), scope=parent["analysis_scope"]
        )
        return True
    except Exception as exc:
        parent["diagnostic_application_status"] = {
            **dict(parent.get("diagnostic_application_status") or {}),
            "application_gate_page_status": "unavailable",
            "application_gate_page_detail": str(exc),
        }
        return False


def _print_parent_proposal_final_pages(pdf_name, parent, title_prefix, manifest, page_prefix):
    """Additive proposed/final weight and before-after diagnostic pages."""
    try:
        import ROOT

        proposal = parent.get("proposed_diagnostic_application_result") or {}
        final = parent.get("final_diagnostic_application_result") or {}
        before = proposal.get("H_MM_nosub_before_pion_subtraction") or final.get("H_MM_nosub_before_pion_subtraction")
        proposed_after = proposal.get("H_MM_nosub_after_pion_subtraction")
        final_after = final.get("H_MM_nosub_after_pion_subtraction")
        proposed_weight = proposal.get("H_pion_weight_vs_MM")
        final_weight = final.get("H_pion_weight_vs_MM")
        if before is None and proposed_weight is None and final_weight is None:
            return False

        canvas = ROOT.TCanvas(
            "pion_parent_transition_{}".format(parent.get("pion_parent_id", "unknown")[-12:]),
            "Pion parent proposal/final",
            950,
            800,
        )
        canvas.Divide(1, 2)
        canvas.cd(1)
        legend_weight = ROOT.TLegend(0.56, 0.70, 0.90, 0.90)
        display_hists = []
        weight_drawn = False
        for histogram, label, color, style in (
            (proposed_weight, "proposed component weight", ROOT.kRed + 1, 2),
            (final_weight, "final applied weight", ROOT.kGreen + 2, 1),
        ):
            if histogram is None:
                continue
            display = clone_root_histogram(
                histogram,
                scope="{}".format(page_prefix),
                role="proposal_final_weight_display",
            )
            display.SetLineColor(color)
            display.SetLineStyle(style)
            display.SetTitle("{} proposed versus final pion weight;MM [GeV];weight".format(title_prefix))
            display.Draw("hist" if not weight_drawn else "hist same")
            display_hists.append(display)
            legend_weight.AddEntry(display, label, "l")
            weight_drawn = True
        if weight_drawn:
            legend_weight.Draw()

        canvas.cd(2)
        legend_spectrum = ROOT.TLegend(0.52, 0.62, 0.90, 0.90)
        spectrum_drawn = False
        input_label = (
            "post-proton/RF input"
            if parent.get("source_target_state") == "post_proton_post_rf"
            else "post-proton input (RF not restored)"
        )
        for histogram, label, color, style in (
            (before, input_label, ROOT.kBlack, 1),
            (
                proposal.get("H_pion_subtraction_template_MM_nosub"),
                "estimated pion contamination",
                ROOT.kRed + 1,
                3,
            ),
            (proposed_after, "proposed component-cleaned", ROOT.kRed + 1, 2),
            (final_after, "final applied/fallback-cleaned", ROOT.kGreen + 2, 1),
        ):
            if histogram is None:
                continue
            display = clone_root_histogram(
                histogram,
                scope="{}".format(page_prefix),
                role="proposal_final_spectrum_display",
            )
            display.SetLineColor(color)
            display.SetLineStyle(style)
            display.SetTitle("{} before / proposed / final;MM [GeV];yield".format(title_prefix))
            display.Draw("hist" if not spectrum_drawn else "hist same")
            display_hists.append(display)
            legend_spectrum.AddEntry(display, label, "l")
            spectrum_drawn = True
        if spectrum_drawn:
            legend_spectrum.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(
            manifest, "{}.proposal_final_transition".format(page_prefix), scope=parent["analysis_scope"]
        )
        return True
    except Exception as exc:
        parent["diagnostic_application_status"] = {
            **dict(parent.get("diagnostic_application_status") or {}),
            "proposal_final_page_status": "unavailable",
            "proposal_final_page_detail": str(exc),
        }
        return False


def _print_authoritative_canonical_t_global_page(pdf_name, global_payload, manifest, title_prefix):
    """Render the canonical-t sum separately from the setting-wide fit."""
    if not isinstance(global_payload, dict):
        return False
    if not bool(global_payload.get("complete")):
        return _print_parent_status_page(
            pdf_name,
            manifest,
            "pion.canonical_t_global.incomplete_status",
            scope="canonical_t_global",
            title="{} AUTHORITATIVE CANONICAL-t GLOBAL INCOMPLETE".format(title_prefix),
            detail="missing products: {}".format(global_payload.get("missing_products") or {}),
        )
    try:
        import ROOT

        canvas = ROOT.TCanvas("pion_canonical_t_global", "Canonical-t pion global", 950, 720)
        legend = ROOT.TLegend(0.50, 0.62, 0.90, 0.90)
        drawn = False
        display_hists = []
        source_states = set(global_payload.get("source_target_states") or [])
        input_label = (
            "input: post-proton/RF"
            if source_states == {"post_proton_post_rf"}
            else "input: post-proton (RF not restored)"
        )
        for histogram, label, color, style in (
            (global_payload.get("H_MM_input"), input_label, ROOT.kBlack, 1),
            (global_payload.get("H_MM_estimated_contamination"), "estimated pion contamination", ROOT.kRed + 1, 3),
            (global_payload.get("H_MM_proposed_clean"), "proposed clean result", ROOT.kOrange + 7, 2),
            (global_payload.get("H_MM_final_clean"), "final applied result", ROOT.kGreen + 2, 1),
        ):
            if histogram is None:
                continue
            display = clone_root_histogram(
                histogram,
                scope="pion_canonical_t_global",
                role="display_{}".format(label.replace(" ", "_")),
            )
            display.SetLineColor(color)
            display.SetLineStyle(style)
            display.SetTitle(
                "{} AUTHORITATIVE CANONICAL-t GLOBAL;MM [GeV];yield".format(title_prefix)
            )
            display.Draw("hist" if not drawn else "hist same")
            legend.AddEntry(display, label, "l")
            display_hists.append(display)
            drawn = True
        if not drawn:
            return False
        legend.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(
            manifest,
            "pion.canonical_t_global.before_estimated_proposed_final",
            scope="canonical_t_global",
        )
        return True
    except Exception:
        return False


def _print_parent_status_page(pdf_name, manifest, page_id, *, scope, title, detail):
    """Emit a real PDF status page for an unavailable required plot slot."""
    try:
        import ROOT

        canvas = ROOT.TCanvas(
            "pion_status_{}".format(page_id.replace(".", "_")[-48:]),
            "Pion diagnostic status",
            900,
            600,
        )
        text = ROOT.TPaveText(0.08, 0.12, 0.92, 0.88, "NDC")
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.SetTextAlign(12)
        text.SetTextSize(0.030)
        text.AddText(str(title))
        for line in textwrap.wrap(str(detail or "unavailable"), width=90):
            text.AddText(line)
        text.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(manifest, page_id, scope=scope)
        return True
    except Exception:
        return False


def _print_parent_summary_page(
    pdf_name, parents, setting_wide_summary, canonical_t_global, manifest, title_prefix
):
    """Add a compact setting-wide versus canonical-t parent table."""
    try:
        import ROOT

        canvas = ROOT.TCanvas("pion_parent_summary", "Pion parent summary", 1100, 700)
        text = ROOT.TPaveText(0.03, 0.05, 0.97, 0.95, "NDC")
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.SetTextAlign(12)
        text.SetTextSize(0.020)
        text.AddText("{} SETTING-WIDE DIAGNOSTIC FIT - NON-AUTHORITATIVE".format(title_prefix))
        text.AddText("AUTHORITATIVE CANONICAL-t GLOBAL = strict sum of frozen canonical-t products")
        if isinstance(canonical_t_global, dict):
            integrals = canonical_t_global.get("integrals") or {}
            text.AddText("canonical-t global integrals: input={} estimate={} proposed={} final={}".format(
                integrals.get("input"),
                integrals.get("estimated_contamination"),
                integrals.get("proposed_clean"),
                integrals.get("final_clean"),
            ))
        text.AddText("scope       A_n      A_delta    A_sidis    S_lambda   chi2/ndf   p       gate/final")
        rows = []
        if isinstance(setting_wide_summary, dict):
            rows.append(("setting", setting_wide_summary))
        rows.extend(("t{}".format(int(parent["t_bin_index"]) + 1), _parent_summary(parent)) for parent in parents)
        for label, row in rows:
            text.AddText("{:<10} {:>8} {:>10} {:>10} {:>10} {:>10} {:>7}  {}/{}".format(
                label,
                str(row.get("A_n")), str(row.get("A_delta")), str(row.get("A_sidis")),
                str(row.get("S_lambda")), str(row.get("chi2_ndf")), str(row.get("p_value")),
                row.get("production_gate_status", "diagnostic"), row.get("final_application_status", row.get("status")),
            ))
        text.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(manifest, "pion.t_bin.summary", scope="setting")
        return True
    except Exception:
        return False


def _print_coordinate_closure_page(
    pdf_name,
    coordinate_audit,
    coordinate_diagnostics,
    manifest,
    title_prefix,
    *,
    axis_label,
    raw_hist_key,
    analysis_hist_key,
    page_id,
):
    """Render the shared kaon-frame closure before parent-specific pages."""
    try:
        import ROOT

        audit = dict(coordinate_audit or {})
        contract = dict(audit.get("kaon_data_coordinate") or {})
        source_accounting = dict(audit.get("source_accounting") or {})
        canvas = ROOT.TCanvas(
            "pion_coordinate_closure_{}".format(
                "{}_{}".format(
                    page_id.replace(".", "_"),
                    str(contract.get("coordinate_fingerprint") or "unknown")[-12:],
                )
            ),
            "Pion coordinate closure",
            1050,
            800,
        )
        canvas.Divide(2, 2)
        for pad_index, source_label in enumerate(
            ("prompt", "rand", "dummy", "dummy_rand"), start=1
        ):
            canvas.cd(pad_index)
            hists = (coordinate_diagnostics or {}).get(source_label) or {}
            raw_hist = hists.get(raw_hist_key)
            analysis_hist = hists.get(analysis_hist_key)
            if raw_hist is not None:
                raw_hist.SetLineColor(ROOT.kBlack)
                raw_hist.SetTitle(
                    "{} {} {} coordinate closure;{};signed yield".format(
                        title_prefix, source_label, axis_label, axis_label
                    )
                )
                raw_hist.Draw("hist")
            if analysis_hist is not None:
                analysis_hist.SetLineColor(ROOT.kBlue + 1)
                analysis_hist.Draw("hist same")
            legend = ROOT.TLegend(0.48, 0.70, 0.88, 0.88)
            if raw_hist is not None:
                legend.AddEntry(raw_hist, "raw experimental {}".format(axis_label), "l")
            if analysis_hist is not None:
                legend.AddEntry(analysis_hist, "kaon analysis {}".format(axis_label), "l")
            closure = (source_accounting.get(source_label) or {}).get(
                "coordinate_closure"
            ) or {}
            legend.SetHeader(
                "shift={:+.4g}, pass={}, migration={}".format(
                    float(contract.get(
                        "mm_shift" if raw_hist_key == "H_MM_raw" else "t_shift", 0.0
                    )),
                    closure.get("passed"),
                    ", ".join(
                        "{}:{}".format(key, value)
                        for key, value in sorted(
                            ((source_accounting.get(source_label) or {}).get(
                                "t_bin_migration"
                            ) or {}).items()
                        )
                    ) or "none",
                ),
                "C",
            )
            legend.Draw()
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(
            manifest, page_id, scope="setting"
        )
        return True
    except Exception as exc:
        return _print_parent_status_page(
            pdf_name,
            manifest,
            page_id,
            scope="setting",
            title="{} pion {} coordinate closure unavailable".format(
                title_prefix, axis_label
            ),
            detail=str(exc),
        )


def _parent_plot_contract(
    page_manifest,
    parents,
    *,
    setting_wide_enabled,
    canonical_parent_k_lambda_render=(),
):
    # Legacy ``protected_fit_or_status`` remains represented by the compact
    # control/protection page's preserved ``pion_control_fit`` manifest slot.
    page_ids = [entry.get("page_id") for entry in page_manifest if isinstance(entry, dict)]
    duplicate_ids = sorted({page_id for page_id in page_ids if page_ids.count(page_id) > 1})
    required = [
        "pion.canonical_t_global.before_estimated_proposed_final",
        "pion.canonical_t_global.incomplete_status",
    ]
    for parent in parents:
        prefix = "pion.t_bin.{}".format(int(parent["t_bin_index"]))
        required.extend((
            "{}.{}".format(prefix, slot)
            for slot in (
                "pion_control_fit", "before_after", "lambda_comparison",
            )
        ))
    # Global is a terminal either/or slot, never both.
    global_complete = "pion.canonical_t_global.before_estimated_proposed_final" in page_ids
    global_incomplete = "pion.canonical_t_global.incomplete_status" in page_ids
    required = [page_id for page_id in required if not page_id.startswith("pion.canonical_t_global")]
    if not (global_complete ^ global_incomplete):
        required.append("pion.canonical_t_global.terminal_exactly_once")
    missing = [page_id for page_id in required if page_id not in page_ids]
    by_scope = {}
    for page_id in page_ids:
        if not isinstance(page_id, str):
            continue
        scope = "setting_wide" if page_id.startswith("pion.setting_wide") else ".".join(page_id.split(".")[:3])
        by_scope[scope] = by_scope.get(scope, 0) + 1
    report = {
        "page_counts": by_scope,
        "mandatory_lambda_pages": len([page_id for page_id in page_ids if str(page_id).endswith(".lambda_comparison")]),
        "mandatory_lambda_expected": len(parents) + (1 if setting_wide_enabled else 0),
        "duplicate_page_ids": duplicate_ids,
        "missing_mandatory_page_ids": missing,
        # These scalar records distinguish an actual comparison renderer
        # success from an explicit unavailable-status page in the same slot.
        "canonical_parent_k_lambda_render": sorted(
            (
                {
                    "t_index": int(record.get("t_index")),
                    "status": str(record.get("status") or "unavailable"),
                    "reason": record.get("reason"),
                    "page_recorded": bool(record.get("page_recorded")),
                }
                for record in canonical_parent_k_lambda_render or ()
                if isinstance(record, dict) and record.get("t_index") is not None
            ),
            key=lambda record: record["t_index"],
        ),
    }
    print("[PLOT CONTRACT]")
    for scope, count in sorted(by_scope.items()):
        print("  {} pion pages: {}".format(scope, count))
    print("  mandatory Lambda pages: {} / expected {}".format(
        report["mandatory_lambda_pages"], report["mandatory_lambda_expected"]
    ))
    print("  duplicate page IDs: {}".format(len(duplicate_ids)))
    print("  missing mandatory page IDs: {}".format(missing))
    return report


def _ensure_required_parent_plot_slots(
    pdf_name, manifest, parents, *, title_prefix, canonical_t_global, setting_wide_enabled,
):
    """Backfill every required slot with a rendered status page if necessary."""
    existing = {
        entry.get("page_id") for entry in manifest if isinstance(entry, dict)
    }
    desired_global = (
        "pion.canonical_t_global.before_estimated_proposed_final"
        if bool((canonical_t_global or {}).get("complete"))
        else "pion.canonical_t_global.incomplete_status"
    )
    desired = [(desired_global, "canonical_t_global")]
    for parent in parents:
        prefix = "pion.t_bin.{}".format(int(parent["t_bin_index"]))
        desired.extend((
            ("{}.{}".format(prefix, slot), parent["analysis_scope"])
            for slot in (
                "pion_control_fit", "before_after", "lambda_comparison",
            )
        ))
    for page_id, scope in desired:
        if page_id not in existing:
            _print_parent_status_page(
                pdf_name,
                manifest,
                page_id,
                scope=scope,
                title="{} diagnostic status".format(title_prefix),
                detail="required plot slot unavailable: {}".format(page_id),
            )
            existing.add(page_id)


def _clone_parent_display(histogram, parent, role):
    """Clone a parent histogram for display without touching its source object."""
    if histogram is None:
        return None
    try:
        source = (
            (parent.get("fit_result") or {}).get(histogram)
            if isinstance(histogram, str)
            else histogram
        )
        if source is None:
            return None
        return clone_root_histogram(
            source,
            scope="pion_parent_t{}".format(int(parent.get("t_bin_index", -1)) + 1),
            role=role,
        )
    except Exception:
        return None


def _draw_parent_overlay(pad, parent, specs, title):
    """Draw a compact display-only overlay and return its cloned objects."""
    try:
        import ROOT
    except Exception:
        return []
    if pad is not None and hasattr(pad, "cd"):
        pad.cd()
    legend = ROOT.TLegend(0.52, 0.64, 0.90, 0.90)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    retained = [legend]
    drawn = False
    for key, label, color, style in specs:
        display = _clone_parent_display((parent.get("fit_result") or {}).get(key), parent, key)
        if display is None:
            continue
        display.SetLineColor(color)
        display.SetLineStyle(style)
        display.SetTitle(title)
        display.Draw("hist" if not drawn else "hist same")
        legend.AddEntry(display, label, "l")
        retained.append(display)
        drawn = True
    if drawn:
        legend.Draw()
    else:
        ROOT.TLatex().DrawLatexNDC(0.12, 0.55, "display inputs unavailable")
    return retained


def _print_parent_control_protection_page(pdf_name, parent, title_prefix, manifest, page_prefix):
    """Combine parent control-fit and protected-fit/status evidence on one page."""
    try:
        import ROOT

        fit_result = parent.get("fit_result") or {}
        diagnostics = fit_result.get("diagnostics") or {}
        protected = diagnostics.get("pi_delta_signal_protected_fit") or (
            diagnostics.get("kaon") or {}
        ).get("pi_delta_signal_protected_fit") or {}
        canvas = ROOT.TCanvas(
            "pion_parent_control_protection_{}".format(parent.get("pion_parent_id", "unknown")[-12:]),
            "Pion parent control and protection",
            1400,
            900,
        )
        canvas.Divide(2, 2)
        retained = []
        retained.extend(_draw_parent_overlay(
            canvas.cd(1), parent,
            (
                ("H_pion_control_input", "pion-control data", ROOT.kBlack, 1),
                ("H_pion_fit_pi_n_scaled", "pi-n", ROOT.kRed + 1, 1),
                ("H_pion_fit_pi_sidis_scaled", "pi-SIDIS", ROOT.kMagenta + 2, 1),
                ("H_pion_fit_pi_delta_scaled", "pi-delta", ROOT.kAzure + 2, 1),
                ("H_pion_fit_total", "total fit", ROOT.kGreen + 2, 2),
            ),
            "{} t{} control fit;MM [GeV];yield".format(title_prefix, int(parent["t_bin_index"]) + 1),
        ))
        retained.extend(_draw_parent_overlay(
            canvas.cd(2), parent,
            (
                ("H_kaon_nosub_input", "kaon no-sub data", ROOT.kBlack, 1),
                ("H_kaon_fit_pi_n_scaled", "pi-n", ROOT.kRed + 1, 1),
                ("H_kaon_fit_pi_sidis_scaled", "pi-SIDIS", ROOT.kMagenta + 2, 1),
                ("H_kaon_fit_pi_delta_scaled", "pi-delta", ROOT.kAzure + 2, 1),
                ("H_kaon_pion_bg_fit_total", "pion background", ROOT.kOrange + 7, 2),
            ),
            "{} t{} protected kaon fit;MM [GeV];yield".format(title_prefix, int(parent["t_bin_index"]) + 1),
        ))
        canvas.cd(3)
        status = ROOT.TPaveText(0.06, 0.08, 0.94, 0.92, "NDC")
        status.SetFillStyle(0)
        status.SetBorderSize(0)
        status.SetTextAlign(12)
        status.SetTextSize(0.030)
        status.AddText("Protected pi-delta / K-Sigma0 status")
        status.AddText("status: {}".format(protected.get("status") or protected.get("fit_variant") or "unavailable"))
        status.AddText("Lambda retention: {}".format((protected.get("lambda_retention") or {}).get("status") or "not recorded"))
        status.AddText("K-Sigma0 template: {}".format(
            "available" if fit_result.get("H_kaon_fit_k_sigma0_scaled") is not None else "unavailable"
        ))
        status.AddText("No protected-signal bin is used as a pion normalization anchor.")
        for line in textwrap.wrap(str(protected.get("reason") or "none"), width=70):
            status.AddText("reason: {}".format(line))
        status.Draw()
        retained.append(status)
        canvas.cd(4)
        gate = parent.get("diagnostic_application_status") or {}
        application = ROOT.TPaveText(0.06, 0.08, 0.94, 0.92, "NDC")
        application.SetFillStyle(0)
        application.SetBorderSize(0)
        application.SetTextAlign(12)
        application.SetTextSize(0.030)
        application.AddText("Authoritative parent state")
        application.AddText("fit status pion/kaon: {} / {}".format(fit_result.get("fit_status_pion"), fit_result.get("fit_status_kaon")))
        application.AddText("production gate: {}".format(gate.get("production_evaluation") or "not recorded"))
        application.AddText("final application: {}".format(gate.get("final_status") or "not recorded"))
        application.Draw()
        retained.append(application)
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(manifest, "{}.pion_control_fit".format(page_prefix), scope=parent["analysis_scope"])
        return True
    except Exception as exc:
        return _print_parent_status_page(
            pdf_name, manifest, "{}.pion_control_fit".format(page_prefix),
            scope=parent["analysis_scope"],
            title="{} control/protection status".format(title_prefix), detail=str(exc),
        )


def _print_parent_application_closure_page(pdf_name, parent, title_prefix, manifest, page_prefix):
    """Combine committed application, MM closure, and parent decision on one page."""
    try:
        import ROOT

        payload = parent.get("final_diagnostic_application_result") or {}
        proposal = parent.get("proposed_diagnostic_application_result") or {}
        fit_result = parent.get("fit_result") or {}
        canvas = ROOT.TCanvas(
            "pion_parent_application_closure_{}".format(parent.get("pion_parent_id", "unknown")[-12:]),
            "Pion parent application and closure",
            1400,
            800,
        )
        canvas.Divide(2, 1)
        canvas.cd(1)
        legend = ROOT.TLegend(0.52, 0.64, 0.90, 0.90)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        retained = [legend]
        drawn = False
        sources = (
            (payload.get("H_MM_nosub_before_pion_subtraction") or proposal.get("H_MM_nosub_before_pion_subtraction") or fit_result.get("H_kaon_nosub_input"), "before pion subtraction", ROOT.kBlack, 1),
            (payload.get("H_pion_subtraction_template_MM_nosub") or proposal.get("H_pion_subtraction_template_MM_nosub"), "estimated pion contamination", ROOT.kOrange + 7, 2),
            (payload.get("H_MM_nosub_after_pion_subtraction") or proposal.get("H_MM_nosub_after_pion_subtraction"), "committed parent result", ROOT.kGreen + 2, 1),
        )
        for index, (histogram, label, color, style) in enumerate(sources):
            display = _clone_parent_display(histogram, parent, "application_{}".format(index))
            if display is None:
                continue
            display.SetLineColor(color)
            display.SetLineStyle(style)
            display.SetTitle("{} t{} application/MM closure;MM [GeV];yield".format(title_prefix, int(parent["t_bin_index"]) + 1))
            display.Draw("hist" if not drawn else "hist same")
            legend.AddEntry(display, label, "l")
            retained.append(display)
            drawn = True
        if drawn:
            legend.Draw()
        else:
            ROOT.TLatex().DrawLatexNDC(0.12, 0.55, "application display inputs unavailable")
        canvas.cd(2)
        diagnostics = payload.get("diagnostics") or {}
        text = ROOT.TPaveText(0.06, 0.08, 0.94, 0.92, "NDC")
        text.SetFillStyle(0)
        text.SetBorderSize(0)
        text.SetTextAlign(12)
        text.SetTextSize(0.030)
        text.AddText("Committed parent application and closure")
        text.AddText("accepted: {}; final status: {}".format(payload.get("accepted"), payload.get("final_application_status") or "not recorded"))
        text.AddText("proposal status: {}".format((parent.get("diagnostic_application_status") or {}).get("proposal_status") or "not recorded"))
        closure = diagnostics.get("model_closure") or {}
        text.AddText("model closure signature: {}; integral ratio: {}".format(closure.get("signature_match", "not recorded"), closure.get("integral_ratio", "not recorded")))
        text.AddText("source target: {}".format(parent.get("source_target_state") or "not recorded"))
        text.Draw()
        retained.append(text)
        canvas.Print(pdf_name)
        canvas.Close()
        _record_parent_page_if_absent(manifest, "{}.before_after".format(page_prefix), scope=parent["analysis_scope"])
        return True
    except Exception as exc:
        return _print_parent_status_page(
            pdf_name, manifest, "{}.before_after".format(page_prefix),
            scope=parent["analysis_scope"],
            title="{} application/MM closure status".format(title_prefix), detail=str(exc),
        )


def _print_parent_lambda_page(pdf_name, parent, inp_dict, title_prefix, manifest, page_prefix):
    """Render one K-Lambda slot and retain its actual presentation outcome."""
    t_index = int(parent["t_bin_index"])
    page_id = "{}.lambda_comparison".format(page_prefix)
    try:
        emitted = print_particle_subtraction_kaon_lambda_comparison_page(
            pdf_name,
            parent.get("fit_result") or {},
            parent.get("final_diagnostic_application_result"),
            title_prefix=title_prefix,
            cut_window=(float(inp_dict["mm_min"]), float(inp_dict["mm_max"])),
            page_manifest=manifest,
            page_id_prefix=page_prefix,
            authoritative=True,
            proposal_payload=parent.get("proposed_diagnostic_application_result"),
        )
    except Exception as exc:
        status_page_recorded = _print_parent_status_page(
            pdf_name, manifest, page_id,
            scope=parent["analysis_scope"],
            title="{} K-Lambda comparison unavailable".format(title_prefix),
            detail="K-Lambda source/reference/normalization unavailable: {}".format(exc),
        )
        return {
            "t_index": t_index,
            "status": "unavailable",
            "reason": str(exc),
            "page_recorded": bool(status_page_recorded),
        }
    if emitted:
        return {
            "t_index": t_index,
            "status": "available",
            "reason": None,
            "page_recorded": True,
        }
    status_page_recorded = _print_parent_status_page(
        pdf_name, manifest, page_id,
        scope=parent["analysis_scope"],
        title="{} K-Lambda comparison unavailable".format(title_prefix),
        detail="K-Lambda comparison renderer returned false.",
    )
    return {
        "t_index": t_index,
        "status": "unavailable",
        "reason": "comparison_renderer_returned_false",
        "page_recorded": bool(status_page_recorded),
    }


def render_setting_t_bin_pion_parent_pages(
    pdf_name,
    parents,
    inp_dict,
    *,
    title_prefix,
    page_manifest,
    setting_wide_summary=None,
    canonical_t_global=None,
    setting_wide_enabled=True,
    coordinate_audit=None,
    coordinate_diagnostics=None,
    coordinate_debug_pdf=None,
):
    """Render authoritative parent pages and route coordinate detail separately."""
    _print_parent_summary_page(
        pdf_name, parents, setting_wide_summary, canonical_t_global, page_manifest, title_prefix
    )
    _print_authoritative_canonical_t_global_page(
        pdf_name, canonical_t_global, page_manifest, title_prefix
    )
    coordinate_pdf = coordinate_debug_pdf or pdf_name
    _print_coordinate_closure_page(
        coordinate_pdf,
        coordinate_audit,
        coordinate_diagnostics,
        page_manifest,
        title_prefix,
        axis_label="MM [GeV]",
        raw_hist_key="H_MM_raw",
        analysis_hist_key="H_MM_analysis",
        page_id="pion.coordinate_mm_closure",
    )
    _print_coordinate_closure_page(
        coordinate_pdf,
        coordinate_audit,
        coordinate_diagnostics,
        page_manifest,
        title_prefix,
        axis_label="|t| [GeV^2]",
        raw_hist_key="H_t_raw",
        analysis_hist_key="H_t_analysis",
        page_id="pion.coordinate_t_closure",
    )
    canonical_parent_k_lambda_render = []
    for parent in parents:
        page_prefix = "pion.t_bin.{}".format(int(parent["t_bin_index"]))
        section_title = "{} t{} [{:.3f}, {:.3f}]".format(
            title_prefix,
            int(parent["t_bin_index"]) + 1,
            float(parent["t_edges"][0]),
            float(parent["t_edges"][1]),
        )
        _print_parent_control_protection_page(
            pdf_name, parent, section_title, page_manifest, page_prefix
        )
        _print_parent_application_closure_page(
            pdf_name, parent, section_title, page_manifest, page_prefix
        )
        canonical_parent_k_lambda_render.append(_print_parent_lambda_page(
            pdf_name, parent, inp_dict, section_title, page_manifest, page_prefix
        ))
    _ensure_required_parent_plot_slots(
        pdf_name,
        page_manifest,
        parents,
        title_prefix=title_prefix,
        canonical_t_global=canonical_t_global,
        setting_wide_enabled=bool(setting_wide_enabled),
    )
    # A failed direct status renderer can still be backfilled into the
    # mandatory parent slot.  This updates only the page-recorded flag; the
    # comparison remains unavailable unless its real renderer succeeded.
    recorded_page_ids = {
        entry.get("page_id") for entry in page_manifest if isinstance(entry, dict)
    }
    for record in canonical_parent_k_lambda_render:
        if record.get("status") != "unavailable":
            continue
        page_id = "pion.t_bin.{}.lambda_comparison".format(record["t_index"])
        record["page_recorded"] = bool(
            record.get("page_recorded") or page_id in recorded_page_ids
        )
    return _parent_plot_contract(
        page_manifest,
        parents,
        setting_wide_enabled=bool(setting_wide_enabled),
        canonical_parent_k_lambda_render=canonical_parent_k_lambda_render,
    )


__all__ = (
    "build_setting_t_bin_pion_parents",
    "render_setting_t_bin_pion_parent_pages",
    "validate_frozen_t_bin_pion_parent_collection",
)
