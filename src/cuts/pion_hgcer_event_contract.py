"""Detached Phase-A event contract for the canonical-t pion subtraction.

This module is deliberately a read-only consumer of already-resolved analysis
objects.  It does not open trees, fit a response, alter a parent, or apply an
HGCer refinement.  Its only job is to attach scalar event provenance to the
existing pion weights and to prove event-level histogram closure.
"""

from __future__ import annotations

import hashlib
import json
import math

from background_config import (
    resolve_particle_subtraction_weight_clip_bounds,
    resolve_particle_subtraction_weight_denominator_floor,
)
from canonical_binning import find_canonical_bin
from pion_component_subtraction import (
    build_simc_shape_pion_control_weights,
    resolve_frozen_parent_application_policy,
    simc_shape_pion_weight_from_value,
)
from root_histogram_ownership import (
    clone_root_histogram,
    fingerprint_histogram_content_error,
)


EVENT_CONTRACT_SCHEMA_VERSION = "pion_hgcer_event_contract/v1"
DEFAULT_CLOSURE_TOLERANCE = 1.0e-10


class EventContractUnavailable(RuntimeError):
    """Expected provenance/policy failure that must not affect production."""


def _json_ready(value):
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        return float(value) if math.isfinite(value) else None
    if isinstance(value, dict):
        return {str(key): _json_ready(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(child) for child in value]
    try:
        if hasattr(value, "item"):
            return _json_ready(value.item())
        if hasattr(value, "tolist"):
            return _json_ready(value.tolist())
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    return float(numeric) if math.isfinite(numeric) else None


def _payload_hash(payload):
    encoded = json.dumps(
        _json_ready(payload),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _finite_or_none(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _normalized_bin_index(value):
    """Normalize frozen bin membership to the public integer-or-None form."""
    if value is None:
        return None
    numeric = _finite_or_none(value)
    if numeric is None or numeric != math.floor(numeric):
        raise EventContractUnavailable("invalid_frozen_bin_index")
    index = int(numeric)
    return index if index >= 0 else None


def _strict_edges(edges, label):
    resolved = []
    for value in edges or ():
        numeric = _finite_or_none(value)
        if numeric is None:
            raise EventContractUnavailable("{}_edges_nonfinite".format(label))
        resolved.append(numeric)
    if len(resolved) < 2 or any(
        resolved[index] >= resolved[index + 1]
        for index in range(len(resolved) - 1)
    ):
        raise EventContractUnavailable("{}_edges_invalid".format(label))
    return resolved


def _geometry(value, edges, coordinate):
    numeric = _finite_or_none(value)
    if numeric is None:
        return None, None, None, "nonfinite_{}".format(coordinate)
    index = int(find_canonical_bin(numeric, edges))
    if index >= 0:
        return index, float(edges[index]), float(edges[index + 1]), "inside"
    if numeric < float(edges[0]):
        reason = "{}_underflow".format(coordinate)
    elif numeric > float(edges[-1]):
        reason = "{}_overflow".format(coordinate)
    else:
        reason = "{}_unassigned".format(coordinate)
    return None, None, None, reason


def _source_kind(label):
    return {
        "prompt": "data_prompt",
        "rand": "data_random",
        "dummy": "dummy_prompt",
        "dummy_rand": "dummy_random",
    }.get(str(label), str(label or "unknown"))


def _hist_snapshot(histogram):
    if histogram is None:
        raise EventContractUnavailable("missing_histogram")
    try:
        nbins = int(histogram.GetNbinsX())
        axis = histogram.GetXaxis()
        edges = [float(axis.GetBinLowEdge(index)) for index in range(1, nbins + 2)]
        bins = [
            {
                "index": int(index),
                "content": float(histogram.GetBinContent(index)),
                "error": float(histogram.GetBinError(index)),
            }
            for index in range(0, nbins + 2)
        ]
    except Exception as exc:
        raise EventContractUnavailable("invalid_histogram:{}".format(exc))
    scalar_values = edges + [
        value
        for entry in bins
        for value in (entry["content"], entry["error"])
    ]
    if not all(math.isfinite(value) for value in scalar_values):
        raise EventContractUnavailable("nonfinite_histogram")
    return {
        "nbins": nbins,
        "axis_min": float(axis.GetXmin()),
        "axis_max": float(axis.GetXmax()),
        "edges": edges,
        "bins": bins,
        "underflow": dict(bins[0]),
        "overflow": dict(bins[-1]),
        "integral": float(sum(entry["content"] for entry in bins[1:-1])),
        "full_integral": float(sum(entry["content"] for entry in bins)),
        "sum_error_squared": float(sum(entry["error"] ** 2 for entry in bins)),
        "fingerprint": fingerprint_histogram_content_error(histogram),
    }


def _histogram_closure(reconstructed, authoritative, tolerance):
    actual = _hist_snapshot(reconstructed)
    expected = _hist_snapshot(authoritative)
    binning_match = bool(
        actual["nbins"] == expected["nbins"]
        and len(actual["edges"]) == len(expected["edges"])
        and all(
            abs(left - right) <= tolerance
            for left, right in zip(actual["edges"], expected["edges"])
        )
        and abs(actual["axis_min"] - expected["axis_min"]) <= tolerance
        and abs(actual["axis_max"] - expected["axis_max"]) <= tolerance
    )
    comparable_bins = zip(actual["bins"], expected["bins"])
    content_differences = [
        abs(left["content"] - right["content"])
        for left, right in comparable_bins
    ] if actual["nbins"] == expected["nbins"] else [float("inf")]
    comparable_bins = zip(actual["bins"], expected["bins"])
    error_differences = [
        abs(left["error"] - right["error"])
        for left, right in comparable_bins
    ] if actual["nbins"] == expected["nbins"] else [float("inf")]
    maximum_content_difference = max(content_differences or (0.0,))
    maximum_error_difference = max(error_differences or (0.0,))
    return {
        "passed": bool(
            binning_match
            and maximum_content_difference <= tolerance
            and maximum_error_difference <= tolerance
            and abs(actual["full_integral"] - expected["full_integral"]) <= tolerance
        ),
        "tolerance": float(tolerance),
        "binning_match": binning_match,
        "maximum_absolute_content_difference": float(maximum_content_difference),
        "maximum_absolute_error_difference": float(maximum_error_difference),
        "integral_difference": float(actual["integral"] - expected["integral"]),
        "full_integral_difference": float(
            actual["full_integral"] - expected["full_integral"]
        ),
        "underflow_content_difference": float(
            actual["underflow"]["content"] - expected["underflow"]["content"]
        ),
        "underflow_error_difference": float(
            actual["underflow"]["error"] - expected["underflow"]["error"]
        ),
        "overflow_content_difference": float(
            actual["overflow"]["content"] - expected["overflow"]["content"]
        ),
        "overflow_error_difference": float(
            actual["overflow"]["error"] - expected["overflow"]["error"]
        ),
        "reconstructed": actual,
        "authoritative": expected,
    }


def _clone_reset(histogram, scope, role):
    return clone_root_histogram(
        histogram,
        scope=scope,
        role=role,
        reset=True,
        sumw2=True,
    )


def _sum_histograms(histograms, scope, role):
    histograms = list(histograms or ())
    if not histograms or any(histogram is None for histogram in histograms):
        raise EventContractUnavailable("missing_histogram_for_{}".format(role))
    total = _clone_reset(histograms[0], scope, role)
    for histogram in histograms:
        total.Add(histogram)
    return total


def _identity_mm_target_map(target_templates, scope, role):
    targets = {}
    for key in ("h_mm", "h_mm_nosub"):
        template = (target_templates or {}).get(key)
        if template is None:
            raise EventContractUnavailable(
                "identity_host_target_template_missing:{}".format(key)
            )
        targets[key] = _clone_reset(template, scope, "{}_{}".format(role, key))
    return targets


def _sum_identity_target_maps(products, field, scope):
    return {
        key: _sum_histograms(
            [(product.get(field) or {}).get(key) for product in products],
            scope,
            "{}_{}".format(field, key),
        )
        for key in ("h_mm", "h_mm_nosub")
    }


def _build_identity_no_proton_cleaning_application(
    *, proton_source_bundle, target_templates, t_edges, delta_edges,
    coordinate_fingerprint, proton_cleaning_result=None,
    closure_tolerance=DEFAULT_CLOSURE_TOLERANCE
):
    """Build the exact noRF identity host from already-prepared kaon records."""
    resolved_t_edges = _strict_edges(t_edges, "canonical_t")
    resolved_delta_edges = _strict_edges(delta_edges, "delta")
    tolerance = float(closure_tolerance)
    coordinate = str(coordinate_fingerprint or "")
    if not coordinate:
        raise EventContractUnavailable("identity_host_coordinate_fingerprint_missing")
    prepared_sources = (proton_source_bundle or {}).get("prepared_sources") or {}
    if not prepared_sources:
        raise EventContractUnavailable("identity_host_prepared_sources_missing")
    upstream_references = {
        key: (target_templates or {}).get(key)
        for key in ("h_mm", "h_mm_nosub")
    }
    if any(reference is None for reference in upstream_references.values()):
        raise EventContractUnavailable("identity_host_upstream_noRF_reference_missing")
    upstream_fingerprints_before = {
        key: fingerprint_histogram_content_error(reference)
        for key, reference in upstream_references.items()
    }

    products = []
    for t_index in range(len(resolved_t_edges) - 1):
        scope = "identity_no_proton_cleaning_t{}".format(t_index + 1)
        products.append({
            "t_index": int(t_index),
            "t_edges": [
                float(resolved_t_edges[t_index]),
                float(resolved_t_edges[t_index + 1]),
            ],
            "coordinate_fingerprint": coordinate,
            "raw_targets": _identity_mm_target_map(
                target_templates, scope, "raw"
            ),
            "proton_targets": _identity_mm_target_map(
                target_templates, scope, "proton"
            ),
            "cleaned_targets_pre_rf": _identity_mm_target_map(
                target_templates, scope, "cleaned"
            ),
            "final_targets": _identity_mm_target_map(
                target_templates, scope, "final"
            ),
            "source_accounting": {},
        })

    lookup = {}
    for source_label, source_spec in sorted(prepared_sources.items()):
        tree_name = str((source_spec or {}).get("tree_name") or "")
        if not tree_name.endswith("_noRF"):
            raise EventContractUnavailable(
                "identity_host_source_not_noRF:{}".format(source_label)
            )
        coefficient = _finite_or_none((source_spec or {}).get("coefficient"))
        if coefficient is None:
            raise EventContractUnavailable(
                "identity_host_source_coefficient_nonfinite:{}".format(source_label)
            )
        for entry_index, entry in sorted(
            ((source_spec or {}).get("entries") or {}).items()
        ):
            signature = "{}:{}".format(str(source_label), int(entry_index))
            t_index, _t_low, _t_high, _t_status = _geometry(
                (entry or {}).get("adj_t"), resolved_t_edges, "canonical_t"
            )
            delta_index, _d_low, _d_high, _d_status = _geometry(
                (entry or {}).get("delta_value"), resolved_delta_edges, "delta"
            )
            lookup[signature] = {
                "event_signature": signature,
                "t_index": t_index,
                "delta_index": delta_index,
                "proton_weight": 0.0,
                "applied_proton_probability": 0.0,
                "cleaned_factor": 1.0,
                "final_cleaned_factor": 1.0,
                "applied_cleaned_factor": 1.0,
                "applied_final_cleaned_factor": 1.0,
                "host_state": "identity_no_proton_cleaning",
            }
            allcuts = bool((entry or {}).get("allcuts"))
            nommcuts = bool((entry or {}).get("nommcuts"))
            if t_index is None or not (allcuts or nommcuts):
                continue
            analysis_mm = _finite_or_none((entry or {}).get("adj_mm"))
            if analysis_mm is None:
                raise EventContractUnavailable(
                    "identity_host_analysis_mm_nonfinite:{}".format(signature)
                )
            product = products[t_index]
            source_metrics = product["source_accounting"].setdefault(
                str(source_label),
                {
                    "coefficient": float(coefficient),
                    "tree_name": tree_name,
                    "selected_records": 0,
                    "allcuts_records": 0,
                    "nommcuts_records": 0,
                    "signed_allcuts_sum": 0.0,
                    "signed_nommcuts_sum": 0.0,
                    "absolute_weight_support": 0.0,
                },
            )
            source_metrics["selected_records"] += 1
            source_metrics["allcuts_records"] += int(allcuts)
            source_metrics["nommcuts_records"] += int(nommcuts)
            source_metrics["absolute_weight_support"] += abs(float(coefficient))
            if nommcuts:
                source_metrics["signed_nommcuts_sum"] += float(coefficient)
            if allcuts:
                source_metrics["signed_allcuts_sum"] += float(coefficient)
            for field in ("raw_targets", "cleaned_targets_pre_rf", "final_targets"):
                targets = product[field]
                if nommcuts:
                    targets["h_mm_nosub"].Fill(analysis_mm, coefficient)
                if allcuts:
                    targets["h_mm"].Fill(analysis_mm, coefficient)

    per_t_closure = []
    for product in products:
        full_closure = _histogram_closure(
            product["final_targets"]["h_mm_nosub"],
            product["raw_targets"]["h_mm_nosub"],
            tolerance,
        )
        cut_closure = _histogram_closure(
            product["final_targets"]["h_mm"],
            product["raw_targets"]["h_mm"],
            tolerance,
        )
        product["final_output_fingerprint"] = (
            fingerprint_histogram_content_error(
                product["final_targets"]["h_mm_nosub"]
            )
        )
        product["identity_host_closure"] = {
            "passed": bool(full_closure["passed"] and cut_closure["passed"]),
            "full": full_closure,
            "cut": cut_closure,
        }
        per_t_closure.append(product["identity_host_closure"])

    global_targets = {
        field: _sum_identity_target_maps(
            products, field, "identity_no_proton_cleaning_global"
        )
        for field in (
            "raw_targets", "proton_targets", "cleaned_targets_pre_rf",
            "final_targets",
        )
    }
    transform_global_full = _histogram_closure(
        global_targets["final_targets"]["h_mm_nosub"],
        global_targets["raw_targets"]["h_mm_nosub"],
        tolerance,
    )
    transform_global_cut = _histogram_closure(
        global_targets["final_targets"]["h_mm"],
        global_targets["raw_targets"]["h_mm"],
        tolerance,
    )
    upstream_noRF_closure = {
        "full": {
            "raw_vs_upstream": _histogram_closure(
                global_targets["raw_targets"]["h_mm_nosub"],
                upstream_references["h_mm_nosub"],
                tolerance,
            ),
            "final_vs_upstream": _histogram_closure(
                global_targets["final_targets"]["h_mm_nosub"],
                upstream_references["h_mm_nosub"],
                tolerance,
            ),
        },
        "cut": {
            "raw_vs_upstream": _histogram_closure(
                global_targets["raw_targets"]["h_mm"],
                upstream_references["h_mm"],
                tolerance,
            ),
            "final_vs_upstream": _histogram_closure(
                global_targets["final_targets"]["h_mm"],
                upstream_references["h_mm"],
                tolerance,
            ),
        },
    }
    upstream_failures = [
        "{}/{}".format(selection, comparison)
        for selection, comparisons in upstream_noRF_closure.items()
        for comparison, closure in comparisons.items()
        if not closure["passed"]
    ]
    upstream_fingerprints_after = {
        key: fingerprint_histogram_content_error(reference)
        for key, reference in upstream_references.items()
    }
    upstream_references_unchanged = bool(
        upstream_fingerprints_before == upstream_fingerprints_after
    )
    upstream_noRF_closure.update({
        "passed": not upstream_failures and upstream_references_unchanged,
        "reference_source": "caller_supplied_preexisting_target_templates",
        "reference_mapping": {
            "h_mm_nosub": "nommcuts",
            "h_mm": "allcuts",
        },
        "upstream_reference_fingerprints_before": upstream_fingerprints_before,
        "upstream_reference_fingerprints_after": upstream_fingerprints_after,
        "upstream_references_unchanged": upstream_references_unchanged,
        "failed_comparisons": upstream_failures,
    })
    transform_closure = {
        "passed": bool(
            all(entry["passed"] for entry in per_t_closure)
            and transform_global_full["passed"]
            and transform_global_cut["passed"]
        ),
        "per_t": per_t_closure,
        "global_full": transform_global_full,
        "global_cut": transform_global_cut,
    }
    result = proton_cleaning_result if isinstance(proton_cleaning_result, dict) else {}
    result_diagnostics = result.get("diagnostics") or {}
    lambda_gate = result_diagnostics.get("lambda_preservation_gate") or {}
    diagnostics = {
        "host_state": "identity_no_proton_cleaning",
        "identity_reason": (
            result.get("fallback_reason")
            or lambda_gate.get("status")
            or "upstream_proton_cleaning_not_committed"
        ),
        "upstream_proton_method": result.get("method"),
        "upstream_proton_result_accepted": bool(result.get("accepted", False)),
        "upstream_fallback_reason": result.get("fallback_reason"),
        "lambda_gate_status": lambda_gate.get("status"),
        "lambda_gate_production_action": lambda_gate.get("production_action"),
        "proton_cleaning_committed": False,
        "rf_applied": False,
        "source_target_state": "post_proton_noRF",
        "event_weight_source": "identity_no_proton_cleaning_prepared_noRF_lookup",
        "coordinate_fingerprint": coordinate,
        "identity_host_closure": {
            "passed": bool(
                transform_closure["passed"]
                and upstream_noRF_closure["passed"]
            ),
            "tolerance": tolerance,
            "identity_transform_closure": transform_closure,
            "upstream_noRF_closure": upstream_noRF_closure,
            "global_constructed_strictly_from_per_t": True,
        },
    }
    if not upstream_references_unchanged:
        raise EventContractUnavailable(
            "identity_host_upstream_noRF_reference_mutated"
        )
    if not upstream_noRF_closure["passed"]:
        raise EventContractUnavailable(
            "identity_host_upstream_noRF_closure_failed:{}".format(
                ",".join(upstream_failures)
            )
        )
    if not transform_closure["passed"]:
        raise EventContractUnavailable("identity_host_histogram_closure_failed")
    return {
        "accepted": True,
        "host_state": "identity_no_proton_cleaning",
        "source_target_state": "post_proton_noRF",
        "rf_restoration_applied": False,
        "coordinate_fingerprint": coordinate,
        "canonical_t_products": tuple(products),
        "raw_targets": global_targets["raw_targets"],
        "proton_targets": global_targets["proton_targets"],
        "cleaned_targets_pre_rf": global_targets["cleaned_targets_pre_rf"],
        "final_targets": global_targets["final_targets"],
        "_prepared_event_weight_lookup": lookup,
        "immutable_record_contract": True,
        "diagnostics": diagnostics,
    }


def _new_metrics():
    return {
        "selected_record_count": 0,
        "record_count": 0,
        "allcuts_record_count": 0,
        "signed_coefficient_sum": 0.0,
        "signed_weighted_sum": 0.0,
        "absolute_weighted_support": 0.0,
        "sumw2": 0.0,
        "effective_entries": 0.0,
    }


def _record_metrics(metrics, record):
    metrics["selected_record_count"] += 1
    metrics["allcuts_record_count"] += int(bool(record.get("allcuts")))
    if not bool(record.get("nommcuts")):
        return
    coefficient = float(record["signed_source_coefficient"])
    contribution = float(record["signed_baseline_event_contribution"])
    metrics["record_count"] += 1
    metrics["signed_coefficient_sum"] += coefficient
    metrics["signed_weighted_sum"] += contribution
    metrics["absolute_weighted_support"] += abs(contribution)
    metrics["sumw2"] += contribution * contribution


def _finalize_metrics(metrics):
    result = dict(metrics)
    sumw2 = float(result["sumw2"])
    result["effective_entries"] = (
        float(result["absolute_weighted_support"]) ** 2 / sumw2
        if sumw2 > 0.0 else 0.0
    )
    return result


def _unavailable(reason, *, stage="validation", exception=None):
    payload = {
        "schema_version": EVENT_CONTRACT_SCHEMA_VERSION,
        "status": "unavailable",
        "available": False,
        "reason": str(reason),
        "diagnostic_stage": str(stage),
        "contract_fingerprint": None,
        "pion_records": [],
        "kaon_host_records": [],
        "pion_closure": {"passed": False},
        "host_closure": {"passed": False},
    }
    if exception is not None:
        payload["exception_type"] = type(exception).__name__
        payload["exception_message"] = str(exception)
    return payload


def _weight_contract(parent, inp_dict, tolerance):
    policy = resolve_frozen_parent_application_policy(parent, inp_dict)
    action = str(policy.get("action") or "error")
    if action not in ("component_weight", "single_scale", "zero"):
        raise EventContractUnavailable(
            "baseline_parent_policy_unavailable:t{}:{}".format(
                int(parent.get("t_bin_index", -1)) + 1, action
            )
        )
    final_payload = parent.get("final_diagnostic_application_result")
    if not isinstance(final_payload, dict):
        raise EventContractUnavailable(
            "baseline_parent_final_payload_missing:t{}".format(
                int(parent.get("t_bin_index", -1)) + 1
            )
        )
    final_status = str(
        (parent.get("final_diagnostic_application_status") or {}).get(
            "final_status"
        )
        or final_payload.get("final_application_status")
        or ""
    )
    expected_status = {
        "component_weight": "applied_component",
        "single_scale": "applied_fallback",
        "zero": "zero",
    }[action]
    if final_status != expected_status:
        raise EventContractUnavailable(
            "baseline_parent_final_status_mismatch:t{}:{}:{}".format(
                int(parent.get("t_bin_index", -1)) + 1,
                action,
                final_status or "missing",
            )
        )
    reference = final_payload.get("H_pion_control_model")
    weights = final_payload.get("weights")
    if reference is None or weights is None:
        raise EventContractUnavailable(
            "baseline_parent_weight_inputs_missing:t{}".format(
                int(parent.get("t_bin_index", -1)) + 1
            )
        )
    weight_values = [float(value) for value in weights]
    rebuild_audit = {
        "performed": False,
        "passed": None,
        "rebuilt_reference_histogram_fingerprint": None,
        "rebuilt_weight_values_fingerprint": None,
    }
    if action == "component_weight":
        clip_min, clip_max = resolve_particle_subtraction_weight_clip_bounds(inp_dict)
        rebuilt_payload = build_simc_shape_pion_control_weights(
            parent.get("fit_result"),
            clip_min=clip_min,
            clip_max=clip_max,
            denom_floor=resolve_particle_subtraction_weight_denominator_floor(inp_dict),
        )
        rebuilt_reference = rebuilt_payload.get("H_pion_control_model")
        rebuilt_weight_values = rebuilt_payload.get("weights")
        rebuilt_weights = list(
            () if rebuilt_weight_values is None else rebuilt_weight_values
        )
        if rebuilt_reference is None or len(weight_values) != len(rebuilt_weights) or any(
            abs(float(left) - float(right)) > tolerance
            for left, right in zip(weight_values, rebuilt_weights)
        ):
            raise EventContractUnavailable(
                "baseline_parent_weight_payload_mismatch:t{}".format(
                    int(parent.get("t_bin_index", -1)) + 1
                )
            )
        frozen_snapshot = _hist_snapshot(reference)
        rebuilt_snapshot = _hist_snapshot(rebuilt_reference)
        if (
            frozen_snapshot["nbins"] != rebuilt_snapshot["nbins"]
            or frozen_snapshot["edges"] != rebuilt_snapshot["edges"]
        ):
            raise EventContractUnavailable(
                "baseline_parent_weight_reference_binning_mismatch:t{}".format(
                    int(parent.get("t_bin_index", -1)) + 1
                )
            )
        rebuild_audit = {
            "performed": True,
            "passed": True,
            "rebuilt_reference_histogram_fingerprint": (
                fingerprint_histogram_content_error(rebuilt_reference)
            ),
            "rebuilt_weight_values_fingerprint": _payload_hash(
                [float(value) for value in rebuilt_weights]
            ),
        }
    if action == "zero" and any(abs(value) > tolerance for value in weight_values):
        raise EventContractUnavailable(
            "baseline_parent_zero_weight_mismatch:t{}".format(
                int(parent.get("t_bin_index", -1)) + 1
            )
        )
    if action == "single_scale":
        scale = _finite_or_none(
            final_payload.get("particle_subtraction_effective_scale")
        )
        if scale is None or any(
            abs(value - scale) > tolerance for value in weight_values
        ):
            raise EventContractUnavailable(
                "baseline_parent_single_scale_weight_mismatch:t{}".format(
                    int(parent.get("t_bin_index", -1)) + 1
                )
            )
    return {
        "action": action,
        "policy": policy,
        "reference": reference,
        "weights": weights,
        "provenance": {
            "pion_parent_id": parent.get("pion_parent_id"),
            "parent_fit_configuration_hash": parent.get(
                "parent_fit_configuration_hash"
            ),
            "application_action": action,
            "fallback_mode": policy.get("fallback_mode"),
            "final_application_status": final_status,
            "authoritative_weight_source": "frozen_final_diagnostic_application_result",
            "reference_histogram_fingerprint": fingerprint_histogram_content_error(
                reference
            ),
            "weight_values_fingerprint": _payload_hash(weight_values),
            "weight_diagnostics": _json_ready(final_payload.get("diagnostics") or {}),
            "component_rebuild_audit": rebuild_audit,
        },
        "authoritative_cut": final_payload.get("H_pion_subtraction_template_MM"),
        "authoritative_full": final_payload.get(
            "H_pion_subtraction_template_MM_nosub"
        ),
    }


def _validate_parent_provenance(
    parent, cache_t, t_index, t_edges, coordinate_fingerprint, canonical_binning
):
    if int(parent.get("t_bin_index", -1)) != int(t_index):
        raise EventContractUnavailable("pion_parent_t_index_mismatch:t{}".format(t_index + 1))
    if [float(value) for value in (parent.get("t_edges") or ())] != list(t_edges):
        raise EventContractUnavailable("pion_parent_t_edges_mismatch:t{}".format(t_index + 1))
    if [float(value) for value in (cache_t.get("t_edges") or ())] != list(t_edges):
        raise EventContractUnavailable("pion_cache_t_edges_mismatch:t{}".format(t_index + 1))
    for owner, fingerprint in (
        ("parent", parent.get("coordinate_fingerprint")),
        ("pion_cache", cache_t.get("coordinate_fingerprint")),
    ):
        if str(fingerprint or "") != coordinate_fingerprint:
            raise EventContractUnavailable(
                "{}_coordinate_fingerprint_mismatch:t{}".format(owner, t_index + 1)
            )
    for key in ("canonical_interval_pair_id", "canonical_interval_pair_hash"):
        expected = canonical_binning.get(key)
        if expected is not None and parent.get(key) != expected:
            raise EventContractUnavailable(
                "pion_parent_{}_mismatch:t{}".format(key, t_index + 1)
            )
    if not str(parent.get("pion_parent_id") or ""):
        raise EventContractUnavailable("pion_parent_id_missing:t{}".format(t_index + 1))


def _build_pion_side(
    *, cache, parents, canonical_t_global, inp_dict, canonical_binning,
    t_edges, delta_edges, tolerance
):
    by_t_cache = tuple(cache.get("by_t") or ())
    if len(by_t_cache) != len(parents) or len(parents) != len(t_edges) - 1:
        raise EventContractUnavailable("pion_parent_cache_count_mismatch")
    coordinate_fingerprint = str(cache.get("coordinate_fingerprint") or "")
    mask_fingerprint = str(
        cache.get("physical_pion_control_mask_fingerprint") or ""
    )
    if not coordinate_fingerprint or not mask_fingerprint:
        raise EventContractUnavailable("pion_cache_provenance_missing")

    all_records = []
    by_t_results = []
    weight_provenance = []
    reconstructed_full = []
    reconstructed_cut = []
    authoritative_full = []
    authoritative_cut = []
    all_component_weight = True
    global_source_metrics = {}

    for t_index, (parent, cache_t) in enumerate(zip(parents, by_t_cache)):
        interval = [float(t_edges[t_index]), float(t_edges[t_index + 1])]
        _validate_parent_provenance(
            parent, cache_t, t_index, interval, coordinate_fingerprint,
            canonical_binning,
        )
        weight_contract = _weight_contract(parent, inp_dict, tolerance)
        all_component_weight = bool(
            all_component_weight and weight_contract["action"] == "component_weight"
        )
        weight_provenance.append(weight_contract["provenance"])
        reference_full = weight_contract["authoritative_full"]
        reference_cut = weight_contract["authoritative_cut"]
        if reference_full is None or reference_cut is None:
            raise EventContractUnavailable(
                "authoritative_pion_template_missing:t{}".format(t_index + 1)
            )
        full_hist = _clone_reset(reference_full, "phase_a_pion_t{}".format(t_index + 1), "full")
        cut_hist = _clone_reset(reference_cut, "phase_a_pion_t{}".format(t_index + 1), "cut")
        source_metrics = {}
        t_records = []
        for source_record in tuple(cache_t.get("records") or ()):
            analysis_t = _finite_or_none(source_record.get("adj_t"))
            resolved_t, t_low, t_high, t_status = _geometry(
                analysis_t, t_edges, "canonical_t"
            )
            if resolved_t != t_index or t_status != "inside":
                raise EventContractUnavailable(
                    "pion_event_analysis_t_assignment_mismatch:{}:{}".format(
                        source_record.get("source_label"), source_record.get("entry_index")
                    )
                )
            cached_t = _normalized_bin_index(source_record.get("t_index"))
            if cached_t != resolved_t:
                raise EventContractUnavailable(
                    "pion_event_cached_t_assignment_mismatch:{}:{}".format(
                        source_record.get("source_label"), source_record.get("entry_index")
                    )
                )
            delta_index, delta_low, delta_high, delta_status = _geometry(
                source_record.get("ssdelta"), delta_edges, "delta"
            )
            normalized_cached_delta = _normalized_bin_index(
                source_record.get("delta_index")
            )
            if normalized_cached_delta != delta_index:
                raise EventContractUnavailable(
                    "pion_event_delta_assignment_mismatch:{}:{}".format(
                        source_record.get("source_label"), source_record.get("entry_index")
                    )
                )
            analysis_mm = _finite_or_none(source_record.get("adj_MM"))
            if analysis_mm is None:
                raise EventContractUnavailable("pion_event_analysis_mm_nonfinite")
            w0 = float(simc_shape_pion_weight_from_value(
                analysis_mm, weight_contract["reference"], weight_contract["weights"]
            ))
            coefficient = float(source_record.get("coefficient", 0.0) or 0.0)
            contribution = coefficient * w0
            source_label = str(source_record.get("source_label") or "unknown")
            source_tree_name = str(source_record.get("source_tree_name") or "")
            source_rf_state = str(
                source_record.get("rf_state")
                or ("noRF" if source_tree_name.endswith("_noRF") else "RF_or_unknown")
            )
            if source_rf_state != "noRF":
                raise EventContractUnavailable(
                    "pion_event_source_not_noRF:{}:{}".format(
                        source_label, source_record.get("entry_index")
                    )
                )
            record = {
                "source_label": source_label,
                "source_kind": _source_kind(source_label),
                "source_tree_name": source_tree_name,
                "source_coefficient": abs(coefficient),
                "signed_source_coefficient": coefficient,
                "allcuts": bool(source_record.get("allcuts")),
                "nommcuts": bool(source_record.get("nommcuts")),
                "noRF_provenance": source_rf_state,
                "input_selection": parent.get("input_selection"),
                "source_target_state": parent.get("source_target_state"),
                "raw_MM": _finite_or_none(source_record.get("raw_MM")),
                "analysis_MM": analysis_mm,
                "raw_t": _finite_or_none(source_record.get("raw_t")),
                "analysis_abs_t": analysis_t,
                "canonical_t_index": int(t_index),
                "canonical_t_lower_edge": t_low,
                "canonical_t_upper_edge": t_high,
                "SHMS_delta": _finite_or_none(source_record.get("ssdelta")),
                "delta_index": delta_index,
                "delta_lower_edge": delta_low,
                "delta_upper_edge": delta_high,
                "refinement_geometry_status": delta_status,
                "P_hgcer_npeSum": _finite_or_none(source_record.get("P_hgcer_npeSum")),
                "P_hgcer_xAtCer": _finite_or_none(source_record.get("P_hgcer_xAtCer")),
                "P_hgcer_yAtCer": _finite_or_none(source_record.get("P_hgcer_yAtCer")),
                "baseline_pion_weight_w0": w0,
                "signed_baseline_event_contribution": contribution,
                "proton_cleaning_factor": None,
                "pion_parent_id": parent.get("pion_parent_id"),
                "parent_physics_identity": _json_ready(
                    parent.get("parent_physics_identity") or {}
                ),
                "parent_fit_configuration_hash": parent.get(
                    "parent_fit_configuration_hash"
                ),
                "baseline_weight_provenance_id": _payload_hash(
                    weight_contract["provenance"]
                ),
                "coordinate_fingerprint": coordinate_fingerprint,
                "physical_pion_control_mask_fingerprint": mask_fingerprint,
                "canonical_interval_pair_id": canonical_binning.get(
                    "canonical_interval_pair_id"
                ),
                "canonical_interval_pair_hash": canonical_binning.get(
                    "canonical_interval_pair_hash"
                ),
            }
            t_records.append(record)
            all_records.append(record)
            metrics = source_metrics.setdefault(source_label, _new_metrics())
            _record_metrics(metrics, record)
            global_metrics = global_source_metrics.setdefault(
                source_label, _new_metrics()
            )
            _record_metrics(global_metrics, record)
            if contribution != 0.0:
                if record["nommcuts"]:
                    full_hist.Fill(analysis_mm, contribution)
                if record["allcuts"]:
                    cut_hist.Fill(analysis_mm, contribution)
        full_closure = _histogram_closure(full_hist, reference_full, tolerance)
        cut_closure = _histogram_closure(cut_hist, reference_cut, tolerance)
        by_t_results.append({
            "canonical_t_index": int(t_index),
            "canonical_t_edges": interval,
            "pion_parent_id": parent.get("pion_parent_id"),
            "application_action": weight_contract["action"],
            "record_count": len(t_records),
            "records_outside_delta_support": sum(
                record["refinement_geometry_status"] != "inside"
                for record in t_records
            ),
            "source_accounting": {
                label: _finalize_metrics(metrics)
                for label, metrics in sorted(source_metrics.items())
            },
            "full_closure": full_closure,
            "cut_closure": cut_closure,
            "passed": bool(full_closure["passed"] and cut_closure["passed"]),
        })
        reconstructed_full.append(full_hist)
        reconstructed_cut.append(cut_hist)
        authoritative_full.append(reference_full)
        authoritative_cut.append(reference_cut)

    global_reconstructed_full = _sum_histograms(
        reconstructed_full, "phase_a_pion_global", "reconstructed_full"
    )
    global_reconstructed_cut = _sum_histograms(
        reconstructed_cut, "phase_a_pion_global", "reconstructed_cut"
    )
    global_authoritative_full = _sum_histograms(
        authoritative_full, "phase_a_pion_global", "authoritative_full"
    )
    global_authoritative_cut = _sum_histograms(
        authoritative_cut, "phase_a_pion_global", "authoritative_cut"
    )
    global_full_closure = _histogram_closure(
        global_reconstructed_full, global_authoritative_full, tolerance
    )
    global_cut_closure = _histogram_closure(
        global_reconstructed_cut, global_authoritative_cut, tolerance
    )
    global_source_accounting = {
        label: _finalize_metrics(metrics)
        for label, metrics in sorted(global_source_metrics.items())
    }
    source_weighted_sum = float(sum(
        metrics["signed_weighted_sum"]
        for metrics in global_source_accounting.values()
    ))
    source_accounting_difference = float(
        source_weighted_sum
        - global_full_closure["reconstructed"]["full_integral"]
    )
    source_accounting_closure = {
        "passed": bool(abs(source_accounting_difference) <= tolerance),
        "tolerance": float(tolerance),
        "signed_weighted_sum_across_sources": source_weighted_sum,
        "reconstructed_global_full_integral": global_full_closure[
            "reconstructed"
        ]["full_integral"],
        "difference": source_accounting_difference,
    }
    named_global_closure = None
    if all_component_weight and isinstance(canonical_t_global, dict):
        named_reference = canonical_t_global.get("H_MM_estimated_contamination")
        if named_reference is not None:
            named_global_closure = _histogram_closure(
                global_reconstructed_full, named_reference, tolerance
            )
    passed = bool(
        all(result["passed"] for result in by_t_results)
        and global_full_closure["passed"]
        and global_cut_closure["passed"]
        and source_accounting_closure["passed"]
        and (named_global_closure is None or named_global_closure["passed"])
    )
    return {
        "records": all_records,
        "by_t": by_t_results,
        "weight_provenance": weight_provenance,
        "source_accounting": global_source_accounting,
        "event_population_fingerprint": _payload_hash(all_records),
        "closure": {
            "passed": passed,
            "tolerance": tolerance,
            "per_t": by_t_results,
            "global_full": global_full_closure,
            "global_cut": global_cut_closure,
            "named_canonical_global": named_global_closure,
            "source_accounting": global_source_accounting,
            "source_accounting_closure": source_accounting_closure,
            "global_constructed_strictly_from_per_t": True,
        },
    }


def _classify_committed_host_state(
    proton_cleaning_result, proton_cleaning_application
):
    result = proton_cleaning_result if isinstance(proton_cleaning_result, dict) else {}
    application = (
        proton_cleaning_application
        if isinstance(proton_cleaning_application, dict) else {}
    )
    result_diagnostics = result.get("diagnostics") or {}
    application_diagnostics = application.get("diagnostics") or {}
    rf_states = [
        diagnostics.get("rf_applied")
        for diagnostics in (result_diagnostics, application_diagnostics)
        if "rf_applied" in diagnostics
    ]
    if not rf_states or any(value is not False for value in rf_states):
        raise EventContractUnavailable(
            "proton_host_rf_restoration_not_explicitly_disabled"
        )
    if (
        "rf_restoration_applied" in application
        and application.get("rf_restoration_applied") is not False
    ):
        raise EventContractUnavailable(
            "proton_host_rf_restoration_not_explicitly_disabled"
        )
    lambda_gate = (
        application_diagnostics.get("lambda_preservation_gate")
        or result_diagnostics.get("lambda_preservation_gate")
        or {}
    )
    committed = lambda_gate.get("proton_cleaning_committed")
    production_action = str(lambda_gate.get("production_action") or "").lower()
    if committed is not None and production_action:
        if bool(committed) != (production_action == "apply"):
            raise EventContractUnavailable(
                "proton_host_committed_state_contradiction"
            )
    if committed is not None:
        host_state = (
            "proton_cleaned"
            if bool(committed) else "identity_no_proton_cleaning"
        )
    elif production_action in ("apply", "bypass"):
        host_state = (
            "proton_cleaned"
            if production_action == "apply" else "identity_no_proton_cleaning"
        )
    else:
        explicit_state = str(application.get("host_state") or "")
        if explicit_state not in (
            "proton_cleaned", "identity_no_proton_cleaning"
        ):
            raise EventContractUnavailable(
                "proton_host_committed_state_unavailable"
            )
        host_state = explicit_state
    explicit_state = str(application.get("host_state") or "")
    if explicit_state and explicit_state != host_state:
        raise EventContractUnavailable(
            "proton_host_committed_state_contradiction"
        )
    return {
        "host_state": host_state,
        "source_target_state": "post_proton_noRF",
        "rf_restoration_applied": False,
        "lambda_gate_status": lambda_gate.get("status"),
        "lambda_gate_production_action": lambda_gate.get("production_action"),
        "proton_cleaning_committed": (
            host_state == "proton_cleaned"
        ),
    }


def _resolve_host_state(proton_cleaning_result, proton_cleaning_application):
    result = proton_cleaning_result if isinstance(proton_cleaning_result, dict) else {}
    application = (
        proton_cleaning_application
        if isinstance(proton_cleaning_application, dict) else {}
    )
    application_diagnostics = application.get("diagnostics") or {}
    result_diagnostics = result.get("diagnostics") or {}
    committed_state = _classify_committed_host_state(result, application)
    host_state = committed_state["host_state"]
    if host_state == "proton_cleaned":
        lookup = result.get("_prepared_event_weight_lookup") or {}
    else:
        lookup = (
            application.get("_prepared_event_weight_lookup")
            or result.get("_prepared_event_weight_lookup")
            or {}
        )
    if not lookup:
        raise EventContractUnavailable(
            "identity_no_proton_cleaning_host_provenance_unavailable"
            if host_state == "identity_no_proton_cleaning"
            else "proton_event_lookup_unavailable"
        )
    return {
        **committed_state,
        "lookup": lookup,
        "lookup_provenance": (
            application_diagnostics.get("event_weight_source")
            or result_diagnostics.get("event_weight_source")
        ),
    }


def _build_host_side(
    *, proton_source_bundle, proton_cleaning_result, proton_cleaning_application,
    pion_parents, coordinate_fingerprint, t_edges, delta_edges, tolerance
):
    host_contract = _resolve_host_state(
        proton_cleaning_result, proton_cleaning_application
    )
    lookup = host_contract["lookup"]
    prepared_sources = (proton_source_bundle or {}).get("prepared_sources") or {}
    products = tuple(
        (proton_cleaning_application or {}).get("canonical_t_products") or ()
    )
    parents = tuple(pion_parents or ())
    if len(products) != len(t_edges) - 1 or len(parents) != len(products):
        raise EventContractUnavailable("proton_canonical_t_product_count_mismatch")
    proton_coordinate = str(
        (proton_cleaning_application or {}).get("coordinate_fingerprint")
        or (proton_cleaning_result or {}).get("coordinate_fingerprint")
        or ""
    )
    if proton_coordinate and proton_coordinate != coordinate_fingerprint:
        raise EventContractUnavailable("proton_coordinate_fingerprint_mismatch")

    reconstructed_full = []
    reconstructed_cut = []
    reference_full = []
    reference_cut = []
    by_t_records = [[] for _ in products]
    all_records = []
    outside_geometry = []
    for source_label, source_spec in sorted(prepared_sources.items()):
        coefficient = float((source_spec or {}).get("coefficient", 0.0) or 0.0)
        tree_name = str((source_spec or {}).get("tree_name") or "")
        if not tree_name.endswith("_noRF"):
            raise EventContractUnavailable(
                "proton_host_source_not_noRF:{}".format(source_label)
            )
        for entry_index, entry in sorted(
            ((source_spec or {}).get("entries") or {}).items()
        ):
            if not bool((entry or {}).get("allcuts")) and not bool(
                (entry or {}).get("nommcuts")
            ):
                continue
            signature = "{}:{}".format(str(source_label), int(entry_index))
            frozen = lookup.get(signature)
            if not isinstance(frozen, dict):
                raise EventContractUnavailable(
                    "proton_event_lookup_missing:{}".format(signature)
                )
            t_index, t_low, t_high, t_status = _geometry(
                (entry or {}).get("adj_t"), t_edges, "canonical_t"
            )
            delta_index, delta_low, delta_high, delta_status = _geometry(
                (entry or {}).get("delta_value"), delta_edges, "delta"
            )
            if _normalized_bin_index(frozen.get("t_index")) != t_index:
                raise EventContractUnavailable(
                    "proton_event_t_assignment_mismatch:{}".format(signature)
                )
            if _normalized_bin_index(frozen.get("delta_index")) != delta_index:
                raise EventContractUnavailable(
                    "proton_event_delta_assignment_mismatch:{}".format(signature)
                )
            cleaned_factor = _finite_or_none(frozen.get("cleaned_factor"))
            final_factor = _finite_or_none(frozen.get("final_cleaned_factor"))
            if cleaned_factor is None or final_factor is None:
                raise EventContractUnavailable(
                    "proton_event_cleaning_factor_missing:{}".format(signature)
                )
            if host_contract["host_state"] == "proton_cleaned":
                if abs(cleaned_factor - final_factor) > tolerance:
                    raise EventContractUnavailable(
                        "proton_event_noRF_factor_alias_mismatch:{}".format(signature)
                    )
            elif (
                abs(cleaned_factor - 1.0) > tolerance
                or abs(final_factor - 1.0) > tolerance
            ):
                raise EventContractUnavailable(
                    "identity_no_proton_cleaning_factor_mismatch:{}".format(signature)
                )
            if host_contract["host_state"] == "identity_no_proton_cleaning":
                applied_probability = _finite_or_none(
                    frozen.get(
                        "applied_proton_probability",
                        frozen.get("proton_weight"),
                    )
                )
                if applied_probability is None or abs(applied_probability) > tolerance:
                    raise EventContractUnavailable(
                        "identity_no_proton_cleaning_probability_mismatch:{}".format(
                            signature
                        )
                    )
            else:
                applied_probability = _finite_or_none(
                    frozen.get(
                        "applied_proton_probability",
                        frozen.get("proton_weight"),
                    )
                )
            contribution = coefficient * final_factor
            record = {
                "source_label": str(source_label),
                "source_kind": _source_kind(source_label),
                "source_tree_name": tree_name,
                "source_coefficient": abs(coefficient),
                "signed_source_coefficient": coefficient,
                "allcuts": bool((entry or {}).get("allcuts")),
                "nommcuts": bool((entry or {}).get("nommcuts")),
                "noRF_provenance": "noRF",
                "input_selection": "no_rf_post_proton_host",
                "source_target_state": host_contract["source_target_state"],
                "host_state": host_contract["host_state"],
                "rf_restoration_applied": False,
                "analysis_MM": _finite_or_none((entry or {}).get("adj_mm")),
                "analysis_abs_t": _finite_or_none((entry or {}).get("adj_t")),
                "canonical_t_index": t_index,
                "canonical_t_lower_edge": t_low,
                "canonical_t_upper_edge": t_high,
                "SHMS_delta": _finite_or_none((entry or {}).get("delta_value")),
                "delta_index": delta_index,
                "delta_lower_edge": delta_low,
                "delta_upper_edge": delta_high,
                "refinement_geometry_status": (
                    delta_status if t_status == "inside" else t_status
                ),
                "proton_cleaning_factor": cleaned_factor,
                "final_cleaned_factor": final_factor,
                "applied_proton_probability": applied_probability,
                "signed_host_event_contribution": contribution,
                "pion_refinement_factor": None,
                "coordinate_fingerprint": coordinate_fingerprint,
                "proton_cleaning_method": (proton_cleaning_result or {}).get("method"),
                "proton_lookup_provenance": host_contract["lookup_provenance"],
                "lambda_gate_status": host_contract.get("lambda_gate_status"),
                "lambda_gate_production_action": host_contract.get(
                    "lambda_gate_production_action"
                ),
                "proton_cleaning_committed": host_contract.get(
                    "proton_cleaning_committed"
                ),
                "proton_final_output_fingerprint": (
                    products[t_index].get("final_output_fingerprint")
                    if t_index is not None else None
                ),
            }
            all_records.append(record)
            if t_index is None:
                outside_geometry.append(record)
            else:
                by_t_records[t_index].append(record)

    by_t_results = []
    for t_index, (product, parent) in enumerate(zip(products, parents)):
        targets = product.get("final_targets") or {}
        full_reference = targets.get("h_mm_nosub")
        cut_reference = targets.get("h_mm")
        parent_reference = parent.get("H_proton_cleaned_final_rf")
        if full_reference is None or cut_reference is None or parent_reference is None:
            raise EventContractUnavailable(
                "proton_final_target_missing:t{}".format(t_index + 1)
            )
        if _normalized_bin_index(product.get("t_index")) != t_index:
            raise EventContractUnavailable(
                "proton_product_t_index_mismatch:t{}".format(t_index + 1)
            )
        product_edges = [float(value) for value in (product.get("t_edges") or ())]
        if product_edges != [float(t_edges[t_index]), float(t_edges[t_index + 1])]:
            raise EventContractUnavailable(
                "proton_product_t_edges_mismatch:t{}".format(t_index + 1)
            )
        if str(product.get("coordinate_fingerprint") or coordinate_fingerprint) != coordinate_fingerprint:
            raise EventContractUnavailable(
                "proton_product_coordinate_fingerprint_mismatch:t{}".format(t_index + 1)
            )
        full_fingerprint = fingerprint_histogram_content_error(full_reference)
        for owner, expected in (
            ("product", product.get("final_output_fingerprint")),
            ("parent", parent.get("proton_output_fingerprint")),
        ):
            if expected is not None and str(expected) != str(full_fingerprint):
                raise EventContractUnavailable(
                    "proton_{}_output_fingerprint_mismatch:t{}".format(
                        owner, t_index + 1
                    )
                )
        parent_host_closure = _histogram_closure(
            full_reference, parent_reference, tolerance
        )
        if not parent_host_closure["passed"]:
            raise EventContractUnavailable(
                "proton_parent_host_mismatch:t{}".format(t_index + 1)
            )
        full_hist = _clone_reset(full_reference, "phase_a_host_t{}".format(t_index + 1), "full")
        cut_hist = _clone_reset(cut_reference, "phase_a_host_t{}".format(t_index + 1), "cut")
        for record in by_t_records[t_index]:
            mm_value = record["analysis_MM"]
            if mm_value is None:
                raise EventContractUnavailable("proton_event_analysis_mm_nonfinite")
            contribution = float(record["signed_host_event_contribution"])
            if contribution == 0.0:
                continue
            if record["nommcuts"]:
                full_hist.Fill(mm_value, contribution)
            if record["allcuts"]:
                cut_hist.Fill(mm_value, contribution)
        full_closure = _histogram_closure(full_hist, full_reference, tolerance)
        cut_closure = _histogram_closure(cut_hist, cut_reference, tolerance)
        by_t_results.append({
            "canonical_t_index": int(t_index),
            "canonical_t_edges": [float(t_edges[t_index]), float(t_edges[t_index + 1])],
            "record_count": len(by_t_records[t_index]),
            "parent_host_closure": parent_host_closure,
            "full_closure": full_closure,
            "cut_closure": cut_closure,
            "passed": bool(full_closure["passed"] and cut_closure["passed"]),
        })
        reconstructed_full.append(full_hist)
        reconstructed_cut.append(cut_hist)
        reference_full.append(full_reference)
        reference_cut.append(cut_reference)

    global_reconstructed_full = _sum_histograms(
        reconstructed_full, "phase_a_host_global", "reconstructed_full"
    )
    global_reconstructed_cut = _sum_histograms(
        reconstructed_cut, "phase_a_host_global", "reconstructed_cut"
    )
    global_reference_full = _sum_histograms(
        reference_full, "phase_a_host_global", "reference_full"
    )
    global_reference_cut = _sum_histograms(
        reference_cut, "phase_a_host_global", "reference_cut"
    )
    global_full_closure = _histogram_closure(
        global_reconstructed_full, global_reference_full, tolerance
    )
    global_cut_closure = _histogram_closure(
        global_reconstructed_cut, global_reference_cut, tolerance
    )
    named_targets = (proton_cleaning_application or {}).get("final_targets") or {}
    named_global_closure = {
        "full": (
            _histogram_closure(global_reference_full, named_targets.get("h_mm_nosub"), tolerance)
            if named_targets.get("h_mm_nosub") is not None else None
        ),
        "cut": (
            _histogram_closure(global_reference_cut, named_targets.get("h_mm"), tolerance)
            if named_targets.get("h_mm") is not None else None
        ),
        "authoritative_for_canonical_sum": False,
    }
    return {
        "records": all_records,
        "outside_geometry": outside_geometry,
        "host_state": host_contract["host_state"],
        "source_target_state": host_contract["source_target_state"],
        "rf_restoration_applied": False,
        "lookup_provenance": host_contract["lookup_provenance"],
        "lambda_gate_status": host_contract.get("lambda_gate_status"),
        "lambda_gate_production_action": host_contract.get(
            "lambda_gate_production_action"
        ),
        "proton_cleaning_committed": host_contract.get(
            "proton_cleaning_committed"
        ),
        "event_population_fingerprint": _payload_hash(all_records),
        "closure": {
            "passed": bool(
                all(result["passed"] for result in by_t_results)
                and global_full_closure["passed"]
                and global_cut_closure["passed"]
            ),
            "tolerance": tolerance,
            "per_t": by_t_results,
            "global_full": global_full_closure,
            "global_cut": global_cut_closure,
            "named_application_global": named_global_closure,
            "global_constructed_strictly_from_per_t": True,
        },
    }


def build_pion_hgcer_event_contract(
    *,
    pion_control_cache,
    pion_parents,
    canonical_t_global,
    proton_source_bundle,
    proton_cleaning_result,
    proton_cleaning_application,
    inp_dict,
    canonical_binning,
    delta_edge_source=None,
    closure_tolerance=DEFAULT_CLOSURE_TOLERANCE,
):
    """Build and validate the detached Phase-A pion/host event contract."""
    try:
        tolerance = float(closure_tolerance)
        if not math.isfinite(tolerance) or tolerance <= 0.0:
            raise EventContractUnavailable("closure_tolerance_invalid")
        cache = pion_control_cache if isinstance(pion_control_cache, dict) else {}
        parents = tuple(pion_parents or ())
        canonical = dict(canonical_binning or {})
        configured_t_edges = canonical.get("t_edges")
        if configured_t_edges is not None:
            t_edges = _strict_edges(configured_t_edges, "canonical_t")
        else:
            derived_edges = []
            for entry in cache.get("by_t") or ():
                interval = list(entry.get("t_edges") or ())
                if len(interval) != 2:
                    raise EventContractUnavailable("canonical_t_edges_missing")
                if not derived_edges:
                    derived_edges.append(float(interval[0]))
                derived_edges.append(float(interval[1]))
            t_edges = _strict_edges(derived_edges, "canonical_t")
        if len(t_edges) != len(parents) + 1:
            raise EventContractUnavailable("canonical_t_parent_edge_count_mismatch")
        delta_edges = _strict_edges(cache.get("delta_edges"), "delta")
        coordinate_fingerprint = str(cache.get("coordinate_fingerprint") or "")

        pion_side = _build_pion_side(
            cache=cache,
            parents=parents,
            canonical_t_global=canonical_t_global,
            inp_dict=inp_dict or {},
            canonical_binning=canonical,
            t_edges=t_edges,
            delta_edges=delta_edges,
            tolerance=tolerance,
        )
        host_side = _build_host_side(
            proton_source_bundle=proton_source_bundle,
            proton_cleaning_result=proton_cleaning_result,
            proton_cleaning_application=proton_cleaning_application,
            pion_parents=parents,
            coordinate_fingerprint=coordinate_fingerprint,
            t_edges=t_edges,
            delta_edges=delta_edges,
            tolerance=tolerance,
        )
        fingerprint_inputs = {
            "schema_version": EVENT_CONTRACT_SCHEMA_VERSION,
            "canonical_t_parent_ids": [
                parent.get("pion_parent_id") for parent in parents
            ],
            "parent_physics_identity": [
                _json_ready(parent.get("parent_physics_identity") or {})
                for parent in parents
            ],
            "canonical_t_edges": t_edges,
            "delta_edges": delta_edges,
            "delta_edge_source": delta_edge_source,
            "coordinate_fingerprint": coordinate_fingerprint,
            "physical_pion_control_mask_fingerprint": cache.get(
                "physical_pion_control_mask_fingerprint"
            ),
            "canonical_interval_pair_id": canonical.get(
                "canonical_interval_pair_id"
            ),
            "canonical_interval_pair_hash": canonical.get(
                "canonical_interval_pair_hash"
            ),
            "source_algebra": [
                {
                    "source_label": label,
                    "coefficient": (entry or {}).get("coefficient"),
                    "tree_name": (entry or {}).get("tree_name"),
                    "rf_state": (entry or {}).get("rf_state"),
                }
                for label, entry in sorted(
                    (cache.get("source_accounting") or {}).items()
                )
            ],
            "baseline_weight_provenance": pion_side["weight_provenance"],
            "pion_source_accounting": pion_side["source_accounting"],
            "pion_event_population_fingerprint": pion_side[
                "event_population_fingerprint"
            ],
            "kaon_host_population_fingerprint": host_side[
                "event_population_fingerprint"
            ],
            "proton_host_state": host_side["host_state"],
            "proton_source_target_state": host_side["source_target_state"],
            "rf_restoration_applied": host_side["rf_restoration_applied"],
            "lambda_gate_status": host_side.get("lambda_gate_status"),
            "lambda_gate_production_action": host_side.get(
                "lambda_gate_production_action"
            ),
            "proton_cleaning_committed": host_side.get(
                "proton_cleaning_committed"
            ),
            "proton_cleaning_method": (proton_cleaning_result or {}).get("method"),
            "proton_lookup_provenance": host_side["lookup_provenance"],
        }
        fingerprint = _payload_hash(fingerprint_inputs)
        contract = {
            "schema_version": EVENT_CONTRACT_SCHEMA_VERSION,
            "status": "available" if (
                pion_side["closure"]["passed"] and host_side["closure"]["passed"]
            ) else "unavailable",
            "available": bool(
                pion_side["closure"]["passed"] and host_side["closure"]["passed"]
            ),
            "reason": None,
            "diagnostic_stage": "complete",
            "contract_fingerprint": fingerprint,
            "fingerprint_inputs": fingerprint_inputs,
            "canonical_t_edges": t_edges,
            "delta_edges": delta_edges,
            "delta_edge_source": delta_edge_source,
            "coordinate_fingerprint": coordinate_fingerprint,
            "physical_pion_control_mask_fingerprint": cache.get(
                "physical_pion_control_mask_fingerprint"
            ),
            "pion_records": pion_side["records"],
            "kaon_host_records": host_side["records"],
            "pion_event_population_fingerprint": pion_side[
                "event_population_fingerprint"
            ],
            "kaon_host_population_fingerprint": host_side[
                "event_population_fingerprint"
            ],
            "baseline_weight_provenance": pion_side["weight_provenance"],
            "pion_source_accounting": pion_side["source_accounting"],
            "pion_closure": pion_side["closure"],
            "host_closure": host_side["closure"],
            "host_records_outside_geometry": host_side["outside_geometry"],
            "host_state": host_side["host_state"],
            "source_target_state": host_side["source_target_state"],
            "rf_restoration_applied": host_side["rf_restoration_applied"],
            "lambda_gate_status": host_side.get("lambda_gate_status"),
            "lambda_gate_production_action": host_side.get(
                "lambda_gate_production_action"
            ),
            "proton_cleaning_committed": host_side.get(
                "proton_cleaning_committed"
            ),
            "immutable_record_contract": True,
            "production_objects_mutated": False,
            "refinement_applied": False,
        }
        if not contract["available"]:
            contract["reason"] = "phase_a_closure_failed"
            contract["diagnostic_stage"] = "histogram_closure"
        # Fail hard here only for accidental non-JSON state.  The returned
        # value itself is the public proof that no ROOT object escaped.
        json.dumps(_json_ready(contract), allow_nan=False)
        return _json_ready(contract)
    except EventContractUnavailable as exc:
        return _unavailable(str(exc), stage="provenance_validation", exception=exc)


def summarize_pion_hgcer_event_contract(contract):
    """Return a lightweight JSON-safe summary without duplicating records."""
    payload = contract if isinstance(contract, dict) else {}
    summary = {
        "schema_version": payload.get("schema_version", EVENT_CONTRACT_SCHEMA_VERSION),
        "status": payload.get("status", "unavailable"),
        "available": bool(payload.get("available", False)),
        "reason": payload.get("reason"),
        "diagnostic_stage": payload.get("diagnostic_stage"),
        "contract_fingerprint": payload.get("contract_fingerprint"),
        "coordinate_fingerprint": payload.get("coordinate_fingerprint"),
        "physical_pion_control_mask_fingerprint": payload.get(
            "physical_pion_control_mask_fingerprint"
        ),
        "canonical_t_edges": payload.get("canonical_t_edges") or [],
        "delta_edges": payload.get("delta_edges") or [],
        "delta_edge_source": payload.get("delta_edge_source"),
        "pion_record_count": len(payload.get("pion_records") or ()),
        "kaon_host_record_count": len(payload.get("kaon_host_records") or ()),
        "host_records_outside_geometry_count": len(
            payload.get("host_records_outside_geometry") or ()
        ),
        "host_state": payload.get("host_state"),
        "source_target_state": payload.get("source_target_state"),
        "rf_restoration_applied": bool(
            payload.get("rf_restoration_applied", False)
        ),
        "lambda_gate_status": payload.get("lambda_gate_status"),
        "lambda_gate_production_action": payload.get(
            "lambda_gate_production_action"
        ),
        "proton_cleaning_committed": payload.get(
            "proton_cleaning_committed"
        ),
        "pion_closure_passed": bool(
            (payload.get("pion_closure") or {}).get("passed", False)
        ),
        "host_closure_passed": bool(
            (payload.get("host_closure") or {}).get("passed", False)
        ),
        "production_objects_mutated": bool(
            payload.get("production_objects_mutated", False)
        ),
        "refinement_applied": bool(payload.get("refinement_applied", False)),
    }
    for key in ("exception_type", "exception_message"):
        if key in payload:
            summary[key] = payload[key]
    json.dumps(_json_ready(summary), allow_nan=False)
    return _json_ready(summary)


__all__ = (
    "DEFAULT_CLOSURE_TOLERANCE",
    "EVENT_CONTRACT_SCHEMA_VERSION",
    "build_pion_hgcer_event_contract",
    "summarize_pion_hgcer_event_contract",
)
