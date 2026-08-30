"""Detached Phase-C Method-B local pion-background closure diagnostics.

Method B compares the frozen Phase-A pion prediction with the actual noRF host
inside canonical (t, delta) cells.  It is diagnostic-only: no result produced
here is a correction factor and no production object is accepted or mutated.
"""

from __future__ import annotations

import bisect
import hashlib
import itertools
import json
import math
from copy import deepcopy
from types import MappingProxyType

from background_config import (
    get_particle_subtraction_component_fit_window_config,
    get_proton_contamination_cleaning_config,
    resolve_particle_subtraction_component_fit_windows,
)
from canonical_binning import find_canonical_bin


METHOD_B_SCHEMA_VERSION = "pion_hgcer_method_b/v1"
METHOD_B_METHOD = "local_pion_background_closure"

DEFAULT_METHOD_B_CONFIG = {
    "mm_regions": [],
    "protected_regions": [],
    "support": {
        "minimum_host_neff": 10.0,
        "minimum_baseline_neff": 10.0,
        "minimum_baseline_significance": 2.0,
        "denominator_absolute_epsilon": 1.0e-12,
    },
    "parent_reference": {
        "minimum_usable_delta_cells": 2,
        "weighting": "inverse_variance",
        "covariance_treatment": "independent_diagnostic_approximation",
    },
    "region_consistency": {
        "primary_sigma": 1.0,
        "marginal_sigma": 2.0,
    },
    "candidate": {
        "combination": "inverse_variance_log_space",
        "single_region_policy": "retain_null",
        "poor_shape_veto": True,
    },
    "shape": {
        "minimum_usable_bins": 2,
        "good_maximum_chi2_ndf": 2.0,
        "good_maximum_abs_pull": 3.0,
        "poor_minimum_chi2_ndf": 5.0,
        "poor_minimum_abs_pull": 5.0,
    },
}


class MethodBUnavailable(RuntimeError):
    """Internal marker for a detached Method-B provenance failure."""

    def __init__(self, reason, stage="validation"):
        super().__init__(str(reason))
        self.reason = str(reason)
        self.stage = str(stage)


def _json_ready(value):
    if isinstance(value, MappingProxyType):
        value = dict(value)
    if isinstance(value, dict):
        return {str(key): _json_ready(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(child) for child in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    try:
        if hasattr(value, "tolist"):
            return _json_ready(value.tolist())
        if hasattr(value, "item"):
            return _json_ready(value.item())
    except Exception:
        pass
    return value


def _payload_hash(value):
    encoded = json.dumps(
        _json_ready(value), sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _deep_merge(base, override):
    result = deepcopy(base)
    for key, value in dict(override or {}).items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = deepcopy(value)
    return result


def _contains_phi_key(value):
    if isinstance(value, dict):
        return any(
            "phi" in str(key).lower() or _contains_phi_key(child)
            for key, child in value.items()
        )
    if isinstance(value, (list, tuple)):
        return any(_contains_phi_key(child) for child in value)
    return False


def _finite(value, name):
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        raise MethodBUnavailable("{}_nonfinite".format(name), "record_validation")
    if not math.isfinite(scalar):
        raise MethodBUnavailable("{}_nonfinite".format(name), "record_validation")
    return scalar


def _strict_edges(values, name):
    try:
        edges = [float(value) for value in values]
    except (TypeError, ValueError):
        raise MethodBUnavailable("{}_edges_missing".format(name), "canonical_geometry")
    if len(edges) < 2 or any(not math.isfinite(value) for value in edges):
        raise MethodBUnavailable("{}_edges_invalid".format(name), "canonical_geometry")
    if any(edges[index] >= edges[index + 1] for index in range(len(edges) - 1)):
        raise MethodBUnavailable("{}_edges_invalid".format(name), "canonical_geometry")
    return edges


def _normalized_index(value):
    if value is None:
        return None
    if isinstance(value, bool):
        raise MethodBUnavailable("stored_assignment_invalid", "record_validation")
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        raise MethodBUnavailable("stored_assignment_invalid", "record_validation")
    integer = int(scalar)
    if not math.isfinite(scalar) or scalar != integer:
        raise MethodBUnavailable("stored_assignment_invalid", "record_validation")
    return integer


def _canonical_index(value, edges):
    scalar = _finite(value, "canonical_coordinate")
    index = int(find_canonical_bin(scalar, tuple(edges)))
    return index if index >= 0 else None


def _window_collection(value, context):
    if not isinstance(value, (list, tuple)):
        raise MethodBUnavailable("{}_invalid".format(context), "configuration")
    if len(value) == 2:
        try:
            low, high = float(value[0]), float(value[1])
        except (TypeError, ValueError):
            pass
        else:
            if math.isfinite(low) and math.isfinite(high) and low < high:
                return [(low, high)]
    windows = []
    for entry in value:
        if not isinstance(entry, (list, tuple)) or len(entry) != 2:
            raise MethodBUnavailable("{}_invalid".format(context), "configuration")
        low, high = float(entry[0]), float(entry[1])
        if not math.isfinite(low) or not math.isfinite(high) or low >= high:
            raise MethodBUnavailable("{}_invalid".format(context), "configuration")
        windows.append((low, high))
    return windows


def resolve_pion_hgcer_method_b_config(
    *, inp_dict, phi_setting, mm_offset_data=0.0
):
    """Resolve Method-B windows from their existing configuration owners.

    ``phi_setting`` is used only by the repository's existing override resolver;
    it is deliberately absent from the returned configuration and fingerprint.
    """
    pion_config = get_particle_subtraction_component_fit_window_config(
        "pion_control", inp_dict=inp_dict, phi_setting=phi_setting
    ) or {}
    resolved_windows = resolve_particle_subtraction_component_fit_windows(
        "pion_control",
        mm_offset_data=mm_offset_data,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
    )
    offset = (
        float(mm_offset_data)
        if bool(pion_config.get("apply_mm_offset_data", False)) else 0.0
    )
    regions = []
    order = ("pi_n", "pi_sidis", "pi_delta")
    for base_name in order:
        windows = list(resolved_windows.get(base_name) or ())
        for index, bounds in enumerate(windows):
            if len(windows) == 1:
                region_name = base_name
            else:
                suffixes = ("low", "high")
                suffix = suffixes[index] if index < len(suffixes) else str(index + 1)
                region_name = "{}_{}".format(base_name, suffix)
            regions.append({
                "region_name": region_name,
                "mm_low": float(bounds[0]),
                "mm_high": float(bounds[1]),
                "region_role": "pion_sensitive",
                "window_source": "pion_control.windows.{}".format(base_name),
                "mm_offset_applied": offset,
            })

    proton_config = get_proton_contamination_cleaning_config(
        inp_dict=inp_dict, phi_setting=phi_setting
    ) or {}
    lambda_bounds = (proton_config.get("validation_windows") or {}).get(
        "lambda_peak"
    )
    lambda_windows = _window_collection(lambda_bounds, "lambda_protected_window")
    if len(lambda_windows) != 1:
        raise MethodBUnavailable("lambda_protected_window_invalid", "configuration")

    kaon_config = get_particle_subtraction_component_fit_window_config(
        "kaon_nosub", inp_dict=inp_dict, phi_setting=phi_setting
    ) or {}
    sigma_bounds = (kaon_config.get("windows") or {}).get("k_sigma0_signal")
    sigma_windows = _window_collection(sigma_bounds, "sigma0_protected_window")
    if len(sigma_windows) != 1:
        raise MethodBUnavailable("sigma0_protected_window_invalid", "configuration")
    sigma_offset = (
        float(mm_offset_data)
        if bool(kaon_config.get("apply_mm_offset_data", False)) else 0.0
    )
    sigma_low, sigma_high = sigma_windows[0]
    protected = [
        {
            "region_name": "KLambda",
            "mm_low": float(lambda_windows[0][0]),
            "mm_high": float(lambda_windows[0][1]),
            "region_role": "protected_signal",
            "window_source": "proton_cleaning.validation_windows.lambda_peak",
            "mm_offset_applied": 0.0,
        },
        {
            "region_name": "KSigma0",
            "mm_low": float(sigma_low + sigma_offset),
            "mm_high": float(sigma_high + sigma_offset),
            "region_role": "protected_signal",
            "window_source": "kaon_nosub.windows.k_sigma0_signal",
            "mm_offset_applied": sigma_offset,
        },
    ]
    resolved = deepcopy(DEFAULT_METHOD_B_CONFIG)
    resolved["mm_regions"] = regions
    resolved["protected_regions"] = protected
    return _json_ready(resolved)


def _normalize_regions(entries, expected_role):
    normalized = []
    names = set()
    for entry in entries or ():
        if not isinstance(entry, dict):
            raise MethodBUnavailable("region_definition_invalid", "configuration")
        name = str(entry.get("region_name") or "")
        if not name or name in names:
            raise MethodBUnavailable("region_name_invalid", "configuration")
        names.add(name)
        low = _finite(entry.get("mm_low"), "region_mm_low")
        high = _finite(entry.get("mm_high"), "region_mm_high")
        if low >= high or str(entry.get("region_role") or "") != expected_role:
            raise MethodBUnavailable("region_definition_invalid", "configuration")
        normalized.append({
            "region_name": name,
            "mm_low": low,
            "mm_high": high,
            "region_role": expected_role,
            "window_source": str(entry.get("window_source") or "unspecified"),
            "mm_offset_applied": _finite(
                entry.get("mm_offset_applied", 0.0), "region_mm_offset"
            ),
        })
    return normalized


def _resolved_config(config):
    if not isinstance(config, dict):
        raise MethodBUnavailable("configuration_missing", "configuration")
    if _contains_phi_key(config):
        raise MethodBUnavailable("phi_dependent_configuration_forbidden", "configuration")
    unknown = set(config) - set(DEFAULT_METHOD_B_CONFIG)
    if unknown:
        raise MethodBUnavailable("configuration_key_unsupported", "configuration")
    resolved = _deep_merge(DEFAULT_METHOD_B_CONFIG, config)
    for section in (
        "support", "parent_reference", "region_consistency", "candidate", "shape"
    ):
        if not isinstance(resolved.get(section), dict):
            raise MethodBUnavailable("configuration_section_invalid", "configuration")
        if set(resolved[section]) - set(DEFAULT_METHOD_B_CONFIG[section]):
            raise MethodBUnavailable("configuration_key_unsupported", "configuration")
    resolved["mm_regions"] = _normalize_regions(
        resolved.get("mm_regions"), "pion_sensitive"
    )
    resolved["protected_regions"] = _normalize_regions(
        resolved.get("protected_regions"), "protected_signal"
    )
    if not resolved["mm_regions"]:
        raise MethodBUnavailable("region_definitions_missing", "configuration")
    protected_names = {entry["region_name"] for entry in resolved["protected_regions"]}
    if not {"KLambda", "KSigma0"}.issubset(protected_names):
        raise MethodBUnavailable("protected_region_definitions_missing", "configuration")

    support = resolved["support"]
    for key in (
        "minimum_host_neff", "minimum_baseline_neff",
        "minimum_baseline_significance", "denominator_absolute_epsilon",
    ):
        support[key] = float(support.get(key, 0.0))
        if not math.isfinite(support[key]) or support[key] <= 0.0:
            raise MethodBUnavailable("support_configuration_invalid", "configuration")
    parent = resolved["parent_reference"]
    parent["minimum_usable_delta_cells"] = int(
        parent.get("minimum_usable_delta_cells", 0)
    )
    if parent["minimum_usable_delta_cells"] < 2:
        raise MethodBUnavailable("parent_reference_configuration_invalid", "configuration")
    if parent.get("weighting") != "inverse_variance" or parent.get(
        "covariance_treatment"
    ) != "independent_diagnostic_approximation":
        raise MethodBUnavailable("parent_reference_configuration_invalid", "configuration")
    consistency = resolved["region_consistency"]
    for key in ("primary_sigma", "marginal_sigma"):
        consistency[key] = float(consistency.get(key, 0.0))
    if not 0.0 < consistency["primary_sigma"] < consistency["marginal_sigma"]:
        raise MethodBUnavailable("region_consistency_configuration_invalid", "configuration")
    candidate = resolved["candidate"]
    if candidate.get("combination") != "inverse_variance_log_space" or candidate.get(
        "single_region_policy"
    ) != "retain_null" or candidate.get("poor_shape_veto") is not True:
        raise MethodBUnavailable("candidate_configuration_invalid", "configuration")
    shape = resolved["shape"]
    shape["minimum_usable_bins"] = int(shape.get("minimum_usable_bins", 0))
    for key in (
        "good_maximum_chi2_ndf", "good_maximum_abs_pull",
        "poor_minimum_chi2_ndf", "poor_minimum_abs_pull",
    ):
        shape[key] = float(shape.get(key, 0.0))
        if not math.isfinite(shape[key]) or shape[key] <= 0.0:
            raise MethodBUnavailable("shape_configuration_invalid", "configuration")
    if (
        shape["minimum_usable_bins"] < 2
        or shape["good_maximum_chi2_ndf"] >= shape["poor_minimum_chi2_ndf"]
        or shape["good_maximum_abs_pull"] >= shape["poor_minimum_abs_pull"]
    ):
        raise MethodBUnavailable("shape_configuration_invalid", "configuration")
    return resolved


def _positive_overlap(left, right):
    return max(float(left["mm_low"]), float(right["mm_low"])) < min(
        float(left["mm_high"]), float(right["mm_high"])
    )


def _annotate_regions(regions, protected):
    result = []
    for region in regions:
        overlaps = [
            entry["region_name"] for entry in protected
            if _positive_overlap(region, entry)
        ]
        result.append({
            **deepcopy(region),
            "protected_signal_overlap": bool(overlaps),
            "overlapping_protected_regions": overlaps,
            "available": not overlaps,
            "reason": "protected_signal_overlap" if overlaps else None,
        })
    return result


def _phase_a_mm_edges(phase):
    pion = (((phase.get("pion_closure") or {}).get("global_full") or {}).get(
        "authoritative"
    ) or {})
    host = (((phase.get("host_closure") or {}).get("global_full") or {}).get(
        "authoritative"
    ) or {})
    pion_edges = _strict_edges(pion.get("edges"), "pion_mm")
    host_edges = _strict_edges(host.get("edges"), "host_mm")
    if pion_edges != host_edges:
        raise MethodBUnavailable("phase_a_mm_binning_mismatch", "mm_binning")
    if int(pion.get("nbins", -1)) != len(pion_edges) - 1 or int(
        host.get("nbins", -1)
    ) != len(host_edges) - 1:
        raise MethodBUnavailable("phase_a_mm_binning_invalid", "mm_binning")
    return pion_edges


def _new_metric():
    return {
        "record_count": 0,
        "signed_yield": 0.0,
        "absolute_weight_support": 0.0,
        "sumw2": 0.0,
    }


def _add_metric(metric, weight):
    metric["record_count"] += 1
    metric["signed_yield"] += weight
    metric["absolute_weight_support"] += abs(weight)
    metric["sumw2"] += weight * weight


def _finish_metric(metric):
    result = dict(metric)
    result["sigma"] = math.sqrt(result["sumw2"])
    result["effective_entries"] = (
        result["absolute_weight_support"] ** 2 / result["sumw2"]
        if result["sumw2"] > 0.0 else 0.0
    )
    return result


def _new_cell(t_index, delta_index, t_edges, delta_edges, regions):
    return {
        "t_index": int(t_index),
        "t_low": float(t_edges[t_index]),
        "t_high": float(t_edges[t_index + 1]),
        "delta_index": int(delta_index),
        "delta_low": float(delta_edges[delta_index]),
        "delta_high": float(delta_edges[delta_index + 1]),
        "host_states": set(),
        "pion_events": [],
        "host_events": [],
        "region_metrics": {
            region["region_name"]: {
                "host": _new_metric(), "baseline": _new_metric()
            }
            for region in regions
        },
    }


def _in_region(mm_value, region):
    return region["mm_low"] <= mm_value < region["mm_high"]


def _record_population(
    records, population, cells, t_edges, delta_edges, regions, audit, host_state
):
    contribution_key = (
        "signed_baseline_event_contribution"
        if population == "pion" else "signed_host_event_contribution"
    )
    for record in records:
        if not isinstance(record, (dict, MappingProxyType)):
            raise MethodBUnavailable("{}_record_invalid".format(population), "record_validation")
        audit["{}_records_seen".format(population)] += 1
        recomputed_t = _canonical_index(record.get("analysis_abs_t"), t_edges)
        recomputed_delta = _canonical_index(record.get("SHMS_delta"), delta_edges)
        stored_t = _normalized_index(record.get("canonical_t_index"))
        stored_delta = _normalized_index(record.get("delta_index"))
        if stored_t != recomputed_t:
            raise MethodBUnavailable(
                "{}_canonical_t_assignment_mismatch".format(population),
                "record_validation",
            )
        if stored_delta != recomputed_delta:
            raise MethodBUnavailable(
                "{}_delta_assignment_mismatch".format(population),
                "record_validation",
            )
        if recomputed_t is None or recomputed_delta is None:
            audit["{}_records_outside_geometry".format(population)] += 1
            continue
        if not bool(record.get("nommcuts")):
            audit["{}_records_not_nommcuts".format(population)] += 1
            continue
        mm_value = _finite(record.get("analysis_MM"), "{}_analysis_mm".format(population))
        weight = _finite(record.get(contribution_key), contribution_key)
        cell = cells[(recomputed_t, recomputed_delta)]
        cell["{}_events".format(population)].append((mm_value, weight))
        if population == "host":
            record_host_state = str(record.get("host_state") or "")
            if record_host_state != host_state:
                raise MethodBUnavailable("host_state_record_mismatch", "record_validation")
            cell["host_states"].add(record_host_state)
        for region in regions:
            if region["available"] and _in_region(mm_value, region):
                key = "baseline" if population == "pion" else "host"
                _add_metric(cell["region_metrics"][region["region_name"]][key], weight)
        audit["{}_nommcuts_records_accepted".format(population)] += 1


def _support_reason(host, baseline, support):
    if baseline["effective_entries"] < support["minimum_baseline_neff"]:
        return "baseline_effective_entries_below_minimum"
    if baseline["signed_yield"] <= support["denominator_absolute_epsilon"]:
        return "baseline_signed_yield_nonpositive"
    significance = (
        baseline["signed_yield"] / baseline["sigma"]
        if baseline["sigma"] > 0.0 else None
    )
    if significance is None or significance < support["minimum_baseline_significance"]:
        return "baseline_cancellation_dominated"
    if host["effective_entries"] < support["minimum_host_neff"]:
        return "host_effective_entries_below_minimum"
    if not math.isfinite(host["signed_yield"]):
        return "host_signed_yield_nonfinite"
    return None


def _region_row(region, metrics, support):
    host = _finish_metric(metrics["host"])
    baseline = _finish_metric(metrics["baseline"])
    row = {
        **deepcopy(region),
        "host_record_count": host["record_count"],
        "host_yield": host["signed_yield"],
        "host_abs_support": host["absolute_weight_support"],
        "host_sumw2": host["sumw2"],
        "host_neff": host["effective_entries"],
        "host_sigma": host["sigma"],
        "baseline_record_count": baseline["record_count"],
        "baseline_pion_yield": baseline["signed_yield"],
        "baseline_pion_abs_support": baseline["absolute_weight_support"],
        "baseline_pion_sumw2": baseline["sumw2"],
        "baseline_pion_neff": baseline["effective_entries"],
        "baseline_pion_sigma": baseline["sigma"],
        "baseline_pion_significance": (
            baseline["signed_yield"] / baseline["sigma"]
            if baseline["sigma"] > 0.0 else None
        ),
        "residual": None,
        "residual_sigma": None,
        "fractional_residual": None,
        "raw_ratio": None,
        "raw_ratio_sigma": None,
        "parent_reference_ratio": None,
        "parent_reference_sigma": None,
        "parent_relative_ratio": None,
        "parent_relative_sigma": None,
        "parent_relative_status": "unavailable",
        "parent_relative_reason": "parent_reference_not_evaluated",
        "support_status": "unavailable",
        "support_reason": region.get("reason"),
    }
    if not region["available"]:
        return row
    reason = _support_reason(host, baseline, support)
    if reason is not None:
        row["support_reason"] = reason
        return row
    host_yield = host["signed_yield"]
    baseline_yield = baseline["signed_yield"]
    residual = host_yield - baseline_yield
    ratio = host_yield / baseline_yield
    ratio_variance = (
        host["sumw2"] / (baseline_yield ** 2)
        + (host_yield ** 2) * baseline["sumw2"] / (baseline_yield ** 4)
    )
    row.update({
        "residual": residual,
        "residual_sigma": math.sqrt(host["sumw2"] + baseline["sumw2"]),
        "fractional_residual": residual / baseline_yield,
        "raw_ratio": ratio,
        "raw_ratio_sigma": math.sqrt(max(0.0, ratio_variance)),
        "support_status": "usable",
        "support_reason": None,
    })
    return row


def _aggregate_rows(rows, prefix):
    return {
        "record_count": sum(row["{}_record_count".format(prefix)] for row in rows),
        "signed_yield": sum(row[
            "baseline_pion_yield" if prefix == "baseline" else "host_yield"
        ] for row in rows),
        "absolute_weight_support": sum(row[
            "baseline_pion_abs_support" if prefix == "baseline" else "host_abs_support"
        ] for row in rows),
        "sumw2": sum(row[
            "baseline_pion_sumw2" if prefix == "baseline" else "host_sumw2"
        ] for row in rows),
    }


def _parent_references(serialized_cells, regions, t_edges, support, parent_config):
    references = []
    by_key = {}
    for t_index in range(len(t_edges) - 1):
        for region in regions:
            region_name = region["region_name"]
            rows = [
                next(entry for entry in cell["regions"] if entry["region_name"] == region_name)
                for cell in serialized_cells if cell["t_index"] == t_index
            ]
            usable = [
                row for row in rows
                if row["support_status"] == "usable"
                and row["raw_ratio_sigma"] is not None
                and row["raw_ratio_sigma"] > 0.0
            ]
            host = _finish_metric(_aggregate_rows(usable, "host"))
            baseline = _finish_metric(_aggregate_rows(usable, "baseline"))
            reference = {
                "t_index": int(t_index),
                "t_low": float(t_edges[t_index]),
                "t_high": float(t_edges[t_index + 1]),
                "region_name": region_name,
                "usable_delta_cell_count": len(usable),
                "contributing_delta_indices": [
                    int(cell["delta_index"])
                    for cell in serialized_cells if cell["t_index"] == t_index
                    for row in cell["regions"]
                    if row["region_name"] == region_name
                    and row["support_status"] == "usable"
                    and row["raw_ratio_sigma"] is not None
                    and row["raw_ratio_sigma"] > 0.0
                ],
                "combined_host_abs_support": host["absolute_weight_support"],
                "combined_host_sumw2": host["sumw2"],
                "combined_host_neff": host["effective_entries"],
                "combined_baseline_abs_support": baseline["absolute_weight_support"],
                "combined_baseline_sumw2": baseline["sumw2"],
                "combined_baseline_neff": baseline["effective_entries"],
                "parent_reference_ratio": None,
                "parent_reference_uncertainty": None,
                "parent_reference_status": "unavailable",
                "parent_reference_reason": None,
                "weighting": parent_config["weighting"],
            }
            if not region["available"]:
                reference["parent_reference_reason"] = "protected_signal_overlap"
            elif len(usable) < parent_config["minimum_usable_delta_cells"]:
                reference["parent_reference_reason"] = "insufficient_usable_delta_cells"
            elif host["effective_entries"] < support["minimum_host_neff"]:
                reference["parent_reference_reason"] = "combined_host_support_below_minimum"
            elif baseline["effective_entries"] < support["minimum_baseline_neff"]:
                reference["parent_reference_reason"] = "combined_baseline_support_below_minimum"
            else:
                inverse_variances = [1.0 / (row["raw_ratio_sigma"] ** 2) for row in usable]
                total_inverse_variance = sum(inverse_variances)
                ratio = sum(
                    row["raw_ratio"] * weight
                    for row, weight in zip(usable, inverse_variances)
                ) / total_inverse_variance
                if math.isfinite(ratio) and ratio > 0.0:
                    reference.update({
                        "parent_reference_ratio": ratio,
                        "parent_reference_uncertainty": math.sqrt(
                            1.0 / total_inverse_variance
                        ),
                        "parent_reference_status": "available",
                        "parent_reference_reason": None,
                    })
                else:
                    reference["parent_reference_reason"] = "parent_reference_nonpositive"
            references.append(reference)
            by_key[(t_index, region_name)] = reference
    return references, by_key


def _apply_parent_references(cells, reference_lookup):
    for cell in cells:
        for row in cell["regions"]:
            reference = reference_lookup[(cell["t_index"], row["region_name"])]
            row["parent_reference_ratio"] = reference["parent_reference_ratio"]
            row["parent_reference_sigma"] = reference[
                "parent_reference_uncertainty"
            ]
            if row["support_status"] != "usable":
                row["parent_relative_reason"] = row["support_reason"]
                continue
            if reference["parent_reference_status"] != "available":
                row["parent_relative_reason"] = reference["parent_reference_reason"]
                continue
            ratio = row["raw_ratio"] / reference["parent_reference_ratio"]
            sigma = math.sqrt(
                (row["raw_ratio_sigma"] / reference["parent_reference_ratio"]) ** 2
                + (
                    row["raw_ratio"] * reference["parent_reference_uncertainty"]
                    / (reference["parent_reference_ratio"] ** 2)
                ) ** 2
            )
            if not math.isfinite(ratio) or not math.isfinite(sigma) or ratio <= 0.0:
                row["parent_relative_reason"] = "nonpositive_parent_relative_ratio"
                continue
            row.update({
                "parent_relative_ratio": ratio,
                "parent_relative_sigma": sigma,
                "parent_relative_status": "available",
                "parent_relative_reason": None,
            })


def _histogram(events, edges):
    nbins = len(edges) - 1
    contents = [0.0] * (nbins + 2)
    sumw2 = [0.0] * (nbins + 2)
    counts = [0] * (nbins + 2)
    for value, weight in events:
        if value < edges[0]:
            index = 0
        elif value > edges[-1]:
            index = nbins + 1
        elif value == edges[-1]:
            index = nbins
        else:
            index = bisect.bisect_right(edges, value)
        contents[index] += weight
        sumw2[index] += weight * weight
        counts[index] += 1
    return contents, sumw2, counts


def _bin_allowed(low, high, regions, protected):
    inside_region = any(
        region["available"]
        and low >= region["mm_low"] and high <= region["mm_high"]
        for region in regions
    )
    protected_overlap = any(
        max(low, entry["mm_low"]) < min(high, entry["mm_high"])
        for entry in protected
    )
    return inside_region and not protected_overlap


def _shape_payload(cell, edges, regions, protected, shape_config, epsilon):
    host, host_sumw2, host_count = _histogram(cell["host_events"], edges)
    baseline, baseline_sumw2, baseline_count = _histogram(cell["pion_events"], edges)
    bins = []
    pulls = []
    nbins = len(edges) - 1
    for index in range(nbins + 2):
        regular = 1 <= index <= nbins
        low = float(edges[index - 1]) if regular else None
        high = float(edges[index]) if regular else None
        allowed = bool(regular and _bin_allowed(low, high, regions, protected))
        variance = host_sumw2[index] + baseline_sumw2[index]
        usable = bool(allowed and math.isfinite(variance) and variance > epsilon)
        residual = host[index] - baseline[index]
        pull = residual / math.sqrt(variance) if usable else None
        if pull is not None:
            pulls.append(pull)
        bins.append({
            "index": int(index),
            "mm_low": low,
            "mm_high": high,
            "underflow": index == 0,
            "overflow": index == nbins + 1,
            "in_allowed_pion_sensitive_domain": allowed,
            "usable_for_shape": usable,
            "host_record_count": host_count[index],
            "host_yield": host[index],
            "host_sumw2": host_sumw2[index],
            "baseline_record_count": baseline_count[index],
            "baseline_yield": baseline[index],
            "baseline_sumw2": baseline_sumw2[index],
            "residual": residual,
            "pull": pull,
        })
    usable_count = len(pulls)
    chi2 = sum(value * value for value in pulls) if pulls else None
    chi2_ndf = chi2 / usable_count if usable_count else None
    maximum = max((abs(value) for value in pulls), default=None)
    if usable_count < shape_config["minimum_usable_bins"]:
        status = "unavailable"
        reason = "insufficient_usable_shape_bins"
    elif (
        chi2_ndf <= shape_config["good_maximum_chi2_ndf"]
        and maximum <= shape_config["good_maximum_abs_pull"]
    ):
        status, reason = "good", None
    elif (
        chi2_ndf > shape_config["poor_minimum_chi2_ndf"]
        or maximum > shape_config["poor_minimum_abs_pull"]
    ):
        status, reason = "poor", "shape_closure_threshold_exceeded"
    else:
        status, reason = "marginal", "shape_closure_marginal"
    return {
        "mm_edges": list(edges),
        "bins": bins,
        "underflow": dict(bins[0]),
        "overflow": dict(bins[-1]),
        "shape_chi2": chi2,
        "shape_ndf": usable_count,
        "shape_chi2_ndf": chi2_ndf,
        "shape_max_abs_pull": maximum,
        "shape_usable_bin_count": usable_count,
        "shape_status": status,
        "shape_reason": reason,
    }


def _region_consistency(rows, config):
    usable = [
        row for row in rows
        if row["parent_relative_status"] == "available"
    ]
    if len(usable) < 2:
        return "insufficient_regions", "fewer_than_two_parent_relative_regions", usable
    primary = config["primary_sigma"]
    marginal = config["marginal_sigma"]
    all_primary = True
    all_marginal = True
    for left, right in itertools.combinations(usable, 2):
        left_value, left_sigma = left["parent_relative_ratio"], left["parent_relative_sigma"]
        right_value, right_sigma = right["parent_relative_ratio"], right["parent_relative_sigma"]
        left_primary = (left_value - primary * left_sigma, left_value + primary * left_sigma)
        right_primary = (right_value - primary * right_sigma, right_value + primary * right_sigma)
        primary_match = (
            left_primary[0] <= right_value <= left_primary[1]
            or right_primary[0] <= left_value <= right_primary[1]
        )
        left_marginal = (left_value - marginal * left_sigma, left_value + marginal * left_sigma)
        right_marginal = (right_value - marginal * right_sigma, right_value + marginal * right_sigma)
        marginal_match = max(left_marginal[0], right_marginal[0]) <= min(
            left_marginal[1], right_marginal[1]
        )
        all_primary = all_primary and primary_match
        all_marginal = all_marginal and marginal_match
    if all_primary:
        return "region_consistent", None, usable
    if all_marginal:
        return "region_marginal", "regions_overlap_only_at_marginal_tolerance", usable
    return "region_inconsistent", "regional_parent_relative_intervals_disjoint", usable


def _candidate(rows, consistency_status, shape_status):
    if consistency_status != "region_consistent":
        return None, None, (
            "single_region_only" if consistency_status == "insufficient_regions"
            else consistency_status
        )
    if shape_status == "poor":
        return None, None, "shape_poor_veto"
    weighted = []
    for row in rows:
        value = row["parent_relative_ratio"]
        sigma = row["parent_relative_sigma"]
        if value is None or sigma is None or value <= 0.0 or sigma <= 0.0:
            continue
        sigma_log = sigma / value
        if math.isfinite(sigma_log) and sigma_log > 0.0:
            weighted.append((math.log(value), 1.0 / (sigma_log * sigma_log)))
    if len(weighted) < 2:
        return None, None, "single_region_only"
    total_weight = sum(weight for _, weight in weighted)
    log_value = sum(value * weight for value, weight in weighted) / total_weight
    log_sigma = math.sqrt(1.0 / total_weight)
    value = math.exp(log_value)
    return value, value * log_sigma, "available_multi_region"


def _validate_phase_a(phase):
    if not phase.get("available") or phase.get("status") != "available":
        raise MethodBUnavailable("phase_a_contract_unavailable", "phase_a_provenance")
    if not str(phase.get("contract_fingerprint") or ""):
        raise MethodBUnavailable("phase_a_contract_fingerprint_missing", "phase_a_provenance")
    if not str(phase.get("coordinate_fingerprint") or ""):
        raise MethodBUnavailable("coordinate_fingerprint_missing", "phase_a_provenance")
    if not (phase.get("pion_closure") or {}).get("passed"):
        raise MethodBUnavailable("phase_a_pion_closure_unavailable", "phase_a_provenance")
    if not (phase.get("host_closure") or {}).get("passed"):
        raise MethodBUnavailable("phase_a_host_closure_unavailable", "phase_a_provenance")
    if not phase.get("baseline_weight_provenance"):
        raise MethodBUnavailable("baseline_weight_provenance_missing", "phase_a_provenance")
    host_state = str(phase.get("host_state") or "")
    if host_state not in ("proton_cleaned", "identity_no_proton_cleaning"):
        raise MethodBUnavailable("host_state_invalid", "phase_a_provenance")
    if phase.get("rf_restoration_applied") is not False:
        raise MethodBUnavailable("phase_a_host_not_norf", "phase_a_provenance")
    return host_state


def _unavailable(reason, stage, config=None, t_edges=None, delta_edges=None, exception=None):
    result = {
        "schema_version": METHOD_B_SCHEMA_VERSION,
        "method": METHOD_B_METHOD,
        "status": "unavailable",
        "available": False,
        "reason": str(reason),
        "diagnostic_stage": str(stage),
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "rf_ct_required": False,
        "interpolation_used": False,
        "phase_a_records_only": True,
        "method_a_numerical_dependency": False,
        "t_edges": list(t_edges or ()),
        "delta_edges": list(delta_edges or ()),
        "mm_regions": list((config or {}).get("mm_regions") or ()),
        "protected_regions": list((config or {}).get("protected_regions") or ()),
        "mm_binning": [],
        "configuration": _json_ready(config or {}),
        "fingerprint": None,
        "cells": [],
        "parent_region_references": [],
        "summary": {},
    }
    if exception is not None:
        result["exception_type"] = type(exception).__name__
        result["exception_message"] = str(exception)
    result = _json_ready(result)
    json.dumps(result, allow_nan=False)
    return result


def build_pion_hgcer_method_b(phase_a_contract, *, config):
    """Build the detached local pion-background closure diagnostic."""
    phase = phase_a_contract if isinstance(phase_a_contract, dict) else {}
    resolved = None
    t_edges = []
    delta_edges = []
    try:
        resolved = _resolved_config(config)
        host_state = _validate_phase_a(phase)
        t_edges = _strict_edges(phase.get("canonical_t_edges"), "canonical_t")
        delta_edges = _strict_edges(phase.get("delta_edges"), "delta")
        mm_edges = _phase_a_mm_edges(phase)
        pion_records = tuple(phase.get("pion_records") or ())
        host_records = tuple(phase.get("kaon_host_records") or ())
        if not pion_records:
            raise MethodBUnavailable("missing_pion_records", "event_population")
        if not host_records:
            raise MethodBUnavailable("missing_host_records", "event_population")
        regions = _annotate_regions(
            resolved["mm_regions"], resolved["protected_regions"]
        )
        cells = {
            (t_index, delta_index): _new_cell(
                t_index, delta_index, t_edges, delta_edges, regions
            )
            for t_index in range(len(t_edges) - 1)
            for delta_index in range(len(delta_edges) - 1)
        }
        audit = {
            "pion_records_seen": 0,
            "pion_records_outside_geometry": 0,
            "pion_records_not_nommcuts": 0,
            "pion_nommcuts_records_accepted": 0,
            "host_records_seen": 0,
            "host_records_outside_geometry": 0,
            "host_records_not_nommcuts": 0,
            "host_nommcuts_records_accepted": 0,
        }
        _record_population(
            pion_records, "pion", cells, t_edges, delta_edges, regions, audit, host_state
        )
        _record_population(
            host_records, "host", cells, t_edges, delta_edges, regions, audit, host_state
        )
        serialized_cells = []
        for key in sorted(cells):
            cell = cells[key]
            cell_host_states = sorted(cell["host_states"])
            if cell_host_states and cell_host_states != [host_state]:
                raise MethodBUnavailable("cell_host_state_mismatch", "record_validation")
            rows = [
                _region_row(
                    region,
                    cell["region_metrics"][region["region_name"]],
                    resolved["support"],
                )
                for region in regions
            ]
            serialized_cells.append({
                "t_index": cell["t_index"],
                "t_low": cell["t_low"],
                "t_high": cell["t_high"],
                "delta_index": cell["delta_index"],
                "delta_low": cell["delta_low"],
                "delta_high": cell["delta_high"],
                "host_state": host_state,
                "host_record_count": len(cell["host_events"]),
                "baseline_record_count": len(cell["pion_events"]),
                "regions": rows,
                "_shape_source": cell,
            })

        references, reference_lookup = _parent_references(
            serialized_cells, regions, t_edges,
            resolved["support"], resolved["parent_reference"],
        )
        _apply_parent_references(serialized_cells, reference_lookup)
        status_counts = {}
        shape_counts = {}
        candidate_counts = {}
        for cell in serialized_cells:
            shape = _shape_payload(
                cell.pop("_shape_source"), mm_edges, regions,
                resolved["protected_regions"], resolved["shape"],
                resolved["support"]["denominator_absolute_epsilon"],
            )
            consistency, consistency_reason, usable_rows = _region_consistency(
                cell["regions"], resolved["region_consistency"]
            )
            candidate, candidate_sigma, candidate_status = _candidate(
                usable_rows, consistency, shape["shape_status"]
            )
            if candidate_status == "available_multi_region":
                method_status, method_reason = "available", None
            elif consistency == "region_inconsistent":
                method_status, method_reason = "internally_inconsistent", consistency_reason
            elif shape["shape_status"] == "poor" and consistency == "region_consistent":
                method_status, method_reason = "shape_inconsistent", shape["shape_reason"]
            elif consistency == "region_marginal":
                method_status, method_reason = "marginal", consistency_reason
            else:
                method_status, method_reason = "unavailable", consistency_reason
            cell.update({
                "region_consistency_status": consistency,
                "region_consistency_reason": consistency_reason,
                **shape,
                "candidate_L_B": candidate,
                "candidate_L_B_uncertainty": candidate_sigma,
                "candidate_L_B_status": candidate_status,
                "method_B_status": method_status,
                "method_B_reason": method_reason,
            })
            status_counts[method_status] = status_counts.get(method_status, 0) + 1
            shape_counts[shape["shape_status"]] = shape_counts.get(
                shape["shape_status"], 0
            ) + 1
            candidate_counts[candidate_status] = candidate_counts.get(
                candidate_status, 0
            ) + 1

        fingerprint_inputs = {
            "schema_version": METHOD_B_SCHEMA_VERSION,
            "method": METHOD_B_METHOD,
            "phase_a_contract_fingerprint": phase.get("contract_fingerprint"),
            "canonical_t_edges": t_edges,
            "delta_edges": delta_edges,
            "coordinate_fingerprint": phase.get("coordinate_fingerprint"),
            "pion_event_population_fingerprint": phase.get(
                "pion_event_population_fingerprint"
            ),
            "host_event_population_fingerprint": phase.get(
                "kaon_host_population_fingerprint"
            ),
            "baseline_weight_provenance": phase.get("baseline_weight_provenance"),
            "host_state": host_state,
            "source_target_state": phase.get("source_target_state"),
            "mm_regions": regions,
            "protected_regions": resolved["protected_regions"],
            "mm_binning": mm_edges,
            "support": resolved["support"],
            "parent_reference": resolved["parent_reference"],
            "region_consistency": resolved["region_consistency"],
            "candidate": resolved["candidate"],
            "shape": resolved["shape"],
            "population_selection": "phase_a_nommcuts",
        }
        result = {
            "schema_version": METHOD_B_SCHEMA_VERSION,
            "method": METHOD_B_METHOD,
            "status": "available",
            "available": True,
            "reason": None,
            "diagnostic_stage": "complete",
            "non_authoritative": True,
            "production_objects_mutated": False,
            "refinement_applied": False,
            "rf_ct_required": False,
            "interpolation_used": False,
            "phase_a_records_only": True,
            "method_a_numerical_dependency": False,
            "phase_a_contract_fingerprint": phase.get("contract_fingerprint"),
            "coordinate_fingerprint": phase.get("coordinate_fingerprint"),
            "host_state": host_state,
            "source_target_state": phase.get("source_target_state"),
            "t_edges": t_edges,
            "delta_edges": delta_edges,
            "mm_regions": regions,
            "protected_regions": resolved["protected_regions"],
            "mm_binning": mm_edges,
            "configuration": resolved,
            "fingerprint_inputs": fingerprint_inputs,
            "fingerprint": _payload_hash(fingerprint_inputs),
            "cells": serialized_cells,
            "parent_region_references": references,
            "summary": {
                **audit,
                "cell_count": len(serialized_cells),
                "parent_region_reference_count": len(references),
                "method_B_status_counts": status_counts,
                "shape_status_counts": shape_counts,
                "candidate_status_counts": candidate_counts,
                "available_region_count": sum(region["available"] for region in regions),
                "protected_overlap_region_count": sum(
                    region["protected_signal_overlap"] for region in regions
                ),
            },
        }
        result = _json_ready(result)
        json.dumps(result, allow_nan=False)
        return result
    except MethodBUnavailable as exc:
        return _unavailable(
            exc.reason, exc.stage, resolved, t_edges, delta_edges, exception=exc
        )
    except Exception as exc:
        return _unavailable(
            "unexpected_method_b_build_failure", "unexpected_exception",
            resolved, t_edges, delta_edges, exception=exc,
        )


def summarize_pion_hgcer_method_b(result):
    """Return a lightweight detached Method-B summary."""
    payload = result if isinstance(result, dict) else {}
    summary = {
        "schema_version": payload.get("schema_version", METHOD_B_SCHEMA_VERSION),
        "method": payload.get("method", METHOD_B_METHOD),
        "status": payload.get("status", "unavailable"),
        "available": bool(payload.get("available", False)),
        "reason": payload.get("reason"),
        "diagnostic_stage": payload.get("diagnostic_stage"),
        "fingerprint": payload.get("fingerprint"),
        "phase_a_contract_fingerprint": payload.get("phase_a_contract_fingerprint"),
        "coordinate_fingerprint": payload.get("coordinate_fingerprint"),
        "host_state": payload.get("host_state"),
        "t_edges": list(payload.get("t_edges") or ()),
        "delta_edges": list(payload.get("delta_edges") or ()),
        "mm_binning": list(payload.get("mm_binning") or ()),
        "cell_count": len(payload.get("cells") or ()),
        "method_B_status_counts": dict(
            (payload.get("summary") or {}).get("method_B_status_counts") or {}
        ),
        "shape_status_counts": dict(
            (payload.get("summary") or {}).get("shape_status_counts") or {}
        ),
        "candidate_status_counts": dict(
            (payload.get("summary") or {}).get("candidate_status_counts") or {}
        ),
        "support_thresholds": dict(
            (payload.get("configuration") or {}).get("support") or {}
        ),
        "rf_ct_required": bool(payload.get("rf_ct_required", False)),
        "interpolation_used": bool(payload.get("interpolation_used", False)),
        "production_objects_mutated": bool(
            payload.get("production_objects_mutated", False)
        ),
        "refinement_applied": bool(payload.get("refinement_applied", False)),
    }
    for key in ("exception_type", "exception_message"):
        if key in payload:
            summary[key] = payload[key]
    summary = _json_ready(summary)
    json.dumps(summary, allow_nan=False)
    return summary


__all__ = [
    "METHOD_B_SCHEMA_VERSION",
    "build_pion_hgcer_method_b",
    "resolve_pion_hgcer_method_b_config",
    "summarize_pion_hgcer_method_b",
]
