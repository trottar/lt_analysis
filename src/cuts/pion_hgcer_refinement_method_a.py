"""Detached Phase-B Method-A pion HGCer response diagnostics.

Method A measures the observed positive HGCer response in canonical (t, delta)
cells.  It is deliberately non-authoritative: it never evaluates or modifies a
pion subtraction weight and it owns no ROOT objects or output writers.
"""

from __future__ import annotations

import hashlib
import json
import math
from copy import deepcopy
from types import MappingProxyType

from canonical_binning import find_canonical_bin


METHOD_A_SCHEMA_VERSION = "pion_hgcer_method_a/v1"
METHOD_A_METHOD = "observed_positive_hgcer_response"
WILSON_Z_95 = 1.959963984540054

DEFAULT_METHOD_A_CONFIG = {
    "positive_response_threshold": 0.0,
    "low_response_upper_threshold": 2.0,
    "uncertainty_method": "wilson_95_percent",
    "quantile_method": "type7_linear",
    "support": {
        "supported_positive_count": 25,
        "supported_low_count": 5,
        "supported_control_count": 5,
        "marginal_positive_count": 10,
        "minimum_control_count_for_ratio": 1,
        "minimum_low_count_for_ratio": 0,
    },
    "signed_support": {
        "minimum_effective_entries": 10.0,
        "minimum_control_significance": 2.0,
        "denominator_absolute_epsilon": 1.0e-12,
    },
    "comparison": {
        "consistent_sigma": 1.0,
        "marginal_sigma": 2.0,
    },
}


class MethodAUnavailable(RuntimeError):
    """Internal marker for a detached diagnostic provenance failure."""

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


def _resolved_config(config):
    if _contains_phi_key(config or {}):
        raise MethodAUnavailable("phi_dependent_configuration_forbidden", "configuration")
    if not isinstance(config, (dict, type(None))):
        raise MethodAUnavailable("configuration_invalid", "configuration")
    unknown_top_level = set(config or {}) - set(DEFAULT_METHOD_A_CONFIG)
    if unknown_top_level:
        raise MethodAUnavailable("configuration_key_unsupported", "configuration")
    for section in ("support", "signed_support", "comparison"):
        override = (config or {}).get(section)
        if override is not None:
            if not isinstance(override, dict) or set(override) - set(DEFAULT_METHOD_A_CONFIG[section]):
                raise MethodAUnavailable("configuration_key_unsupported", "configuration")
    resolved = _deep_merge(DEFAULT_METHOD_A_CONFIG, config)
    if _contains_phi_key(resolved):
        raise MethodAUnavailable("phi_dependent_configuration_forbidden", "configuration")
    if float(resolved.get("positive_response_threshold", -1.0)) != 0.0:
        raise MethodAUnavailable("positive_response_threshold_must_be_zero", "configuration")
    if float(resolved.get("low_response_upper_threshold", -1.0)) != 2.0:
        raise MethodAUnavailable("low_response_threshold_must_be_two", "configuration")
    if resolved.get("uncertainty_method") != "wilson_95_percent":
        raise MethodAUnavailable("uncertainty_method_must_be_wilson_95_percent", "configuration")
    if resolved.get("quantile_method") != "type7_linear":
        raise MethodAUnavailable("quantile_method_must_be_type7_linear", "configuration")

    support = resolved["support"]
    integer_keys = (
        "supported_positive_count",
        "supported_low_count",
        "supported_control_count",
        "marginal_positive_count",
        "minimum_control_count_for_ratio",
        "minimum_low_count_for_ratio",
    )
    for key in integer_keys:
        value = int(support.get(key, -1))
        if value < 0:
            raise MethodAUnavailable("invalid_support_configuration", "configuration")
        support[key] = value
    if support["supported_positive_count"] < support["marginal_positive_count"]:
        raise MethodAUnavailable("invalid_support_configuration", "configuration")
    signed = resolved["signed_support"]
    for key in (
        "minimum_effective_entries",
        "minimum_control_significance",
        "denominator_absolute_epsilon",
    ):
        signed[key] = float(signed.get(key, 0.0))
        if not math.isfinite(signed[key]) or signed[key] <= 0.0:
            raise MethodAUnavailable("invalid_signed_support_configuration", "configuration")
    comparison = resolved["comparison"]
    comparison["consistent_sigma"] = float(comparison.get("consistent_sigma", 0.0))
    comparison["marginal_sigma"] = float(comparison.get("marginal_sigma", 0.0))
    if not 0.0 < comparison["consistent_sigma"] <= comparison["marginal_sigma"]:
        raise MethodAUnavailable("invalid_comparison_configuration", "configuration")
    return resolved


def _strict_edges(values, name):
    try:
        edges = [float(value) for value in values]
    except (TypeError, ValueError):
        raise MethodAUnavailable("{}_edges_missing".format(name), "canonical_geometry")
    if len(edges) < 2 or any(not math.isfinite(value) for value in edges):
        raise MethodAUnavailable("{}_edges_invalid".format(name), "canonical_geometry")
    if any(edges[index] >= edges[index + 1] for index in range(len(edges) - 1)):
        raise MethodAUnavailable("{}_edges_invalid".format(name), "canonical_geometry")
    return edges


def _canonical_index(value, edges):
    scalar = _finite(value, "canonical_coordinate")
    index = int(find_canonical_bin(scalar, tuple(edges)))
    return index if index >= 0 else None


def _stored_index(value, name):
    if value is None:
        return None
    if isinstance(value, bool):
        raise MethodAUnavailable("{}_assignment_invalid".format(name), "record_validation")
    try:
        integer = int(value)
        scalar = float(value)
    except (TypeError, ValueError):
        raise MethodAUnavailable("{}_assignment_invalid".format(name), "record_validation")
    if not math.isfinite(scalar) or scalar != float(integer):
        raise MethodAUnavailable("{}_assignment_invalid".format(name), "record_validation")
    return integer


def _finite(value, name):
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        raise MethodAUnavailable("nonfinite_{}".format(name), "record_validation")
    if not math.isfinite(scalar):
        raise MethodAUnavailable("nonfinite_{}".format(name), "record_validation")
    return scalar


def _type7_quantile(values, probability):
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return None
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * float(probability)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def _distribution(values):
    samples = [float(value) for value in values]
    if not samples:
        return {
            "mean": None,
            "rms": None,
            "minimum": None,
            "maximum": None,
            "q10": None,
            "q25": None,
            "q50": None,
            "q75": None,
            "q90": None,
        }
    mean = sum(samples) / float(len(samples))
    rms = math.sqrt(sum((value - mean) ** 2 for value in samples) / float(len(samples)))
    return {
        "mean": mean,
        "rms": rms,
        "minimum": min(samples),
        "maximum": max(samples),
        "q10": _type7_quantile(samples, 0.10),
        "q25": _type7_quantile(samples, 0.25),
        "q50": _type7_quantile(samples, 0.50),
        "q75": _type7_quantile(samples, 0.75),
        "q90": _type7_quantile(samples, 0.90),
    }


def _wilson_interval(successes, total):
    successes = int(successes)
    total = int(total)
    if total <= 0 or successes < 0 or successes > total:
        return None, None
    fraction = float(successes) / float(total)
    z2 = WILSON_Z_95 * WILSON_Z_95
    denominator = 1.0 + z2 / total
    center = (fraction + z2 / (2.0 * total)) / denominator
    spread = (
        WILSON_Z_95
        * math.sqrt(fraction * (1.0 - fraction) / total + z2 / (4.0 * total * total))
        / denominator
    )
    return max(0.0, center - spread), min(1.0, center + spread)


def _ratio_from_fraction(value):
    if value is None or value >= 1.0:
        return None
    return float(value) / (1.0 - float(value))


def _support_class(positive, low, control, config):
    support = config["support"]
    if (
        positive >= support["supported_positive_count"]
        and low >= support["supported_low_count"]
        and control >= support["supported_control_count"]
    ):
        return "supported"
    if (
        positive >= support["marginal_positive_count"]
        and control >= support["minimum_control_count_for_ratio"]
        and low >= support["minimum_low_count_for_ratio"]
    ):
        return "marginal"
    return "unsupported"


def _prompt_metrics(records, config):
    positive = len(records)
    low_records = [record for record in records if record["npe"] <= 2.0]
    control_records = [record for record in records if record["npe"] > 2.0]
    low = len(low_records)
    control = len(control_records)
    support_class = _support_class(positive, low, control, config)
    available = support_class != "unsupported"
    fraction = float(low) / float(positive) if positive else None
    interval_low, interval_high = _wilson_interval(low, positive)
    ratio = float(low) / float(control) if control else None
    npe = _distribution([record["npe"] for record in records])
    x = _distribution([record["x"] for record in records])
    y = _distribution([record["y"] for record in records])
    if not available:
        fraction = interval_low = interval_high = ratio = None
        npe = _distribution(())
        x = _distribution(())
        y = _distribution(())
    return {
        "positive": positive,
        "low": low,
        "control": control,
        "partition_closure_passed": positive == low + control,
        "support_class": support_class,
        "available": available,
        "f_low": fraction,
        "f_low_low": interval_low if available else None,
        "f_low_high": interval_high if available else None,
        "ratio": ratio,
        "ratio_low": _ratio_from_fraction(interval_low) if available else None,
        "ratio_high": _ratio_from_fraction(interval_high) if available else None,
        "npe": npe,
        "x": x,
        "y": y,
    }


def _empty_signed():
    return {"sumw": 0.0, "sumabs": 0.0, "sumw2": 0.0, "records": 0}


def _add_signed(target, weight):
    target["sumw"] += weight
    target["sumabs"] += abs(weight)
    target["sumw2"] += weight * weight
    target["records"] += 1


def _finish_signed(target):
    result = dict(target)
    result["neff"] = (
        result["sumabs"] ** 2 / result["sumw2"] if result["sumw2"] > 0.0 else 0.0
    )
    return result


def _signed_ratio(signed, config):
    positive = signed["positive"]
    low = signed["low"]
    control = signed["control"]
    thresholds = config["signed_support"]
    uncertainty = math.sqrt(control["sumw2"])
    denominator_floor = max(
        thresholds["denominator_absolute_epsilon"],
        thresholds["minimum_control_significance"] * uncertainty,
    )
    if (
        positive["neff"] < thresholds["minimum_effective_entries"]
        or control["neff"] < thresholds["minimum_effective_entries"]
        or control["sumw"] <= denominator_floor
    ):
        return None, None
    ratio = low["sumw"] / control["sumw"]
    variance = (
        low["sumw2"] / (control["sumw"] ** 2)
        + (low["sumw"] ** 2) * control["sumw2"] / (control["sumw"] ** 4)
    )
    sigma = math.sqrt(max(0.0, variance))
    return ratio, sigma


def _prompt_ratio_sigma(metrics):
    if not metrics["available"] or metrics["ratio"] is None:
        return None
    total = metrics["positive"]
    fraction = metrics["f_low"]
    if total <= 0 or fraction is None or fraction >= 1.0:
        return None
    fraction_variance = fraction * (1.0 - fraction) / float(total)
    return math.sqrt(max(0.0, fraction_variance)) / ((1.0 - fraction) ** 2)


def _prompt_vs_signed_status(prompt, signed_ratio, signed_sigma, config):
    if signed_ratio is None or signed_sigma is None or prompt["ratio"] is None:
        return "signed_unavailable"
    prompt_sigma = _prompt_ratio_sigma(prompt)
    if prompt_sigma is None:
        return "signed_unavailable"
    combined = math.sqrt(prompt_sigma * prompt_sigma + signed_sigma * signed_sigma)
    difference = abs(prompt["ratio"] - signed_ratio)
    if combined == 0.0:
        return "consistent" if difference == 0.0 else "inconsistent"
    significance = difference / combined
    comparison = config["comparison"]
    if significance <= comparison["consistent_sigma"]:
        return "consistent"
    if significance <= comparison["marginal_sigma"]:
        return "marginal"
    return "inconsistent"


def _nommcuts_vs_allcuts_status(nommcuts, allcuts):
    if not nommcuts["available"] or not allcuts["available"]:
        return "allcuts_unavailable"
    mutual_centers = (
        allcuts["f_low_low"] <= nommcuts["f_low"] <= allcuts["f_low_high"]
        and nommcuts["f_low_low"] <= allcuts["f_low"] <= nommcuts["f_low_high"]
    )
    if mutual_centers:
        return "consistent"
    if max(nommcuts["f_low_low"], allcuts["f_low_low"]) <= min(
        nommcuts["f_low_high"], allcuts["f_low_high"]
    ):
        return "marginal"
    return "inconsistent"


def _empty_cell(t_index, delta_index, t_edges, delta_edges):
    return {
        "t_index": int(t_index),
        "t_low": t_edges[t_index],
        "t_high": t_edges[t_index + 1],
        "delta_index": int(delta_index),
        "delta_low": delta_edges[delta_index],
        "delta_high": delta_edges[delta_index + 1],
        "prompt_nommcuts": [],
        "prompt_allcuts": [],
        "signed": {
            "positive": _empty_signed(),
            "low": _empty_signed(),
            "control": _empty_signed(),
        },
    }


def _validate_provenance(response, phase):
    if response.get("status") != "available":
        raise MethodAUnavailable("response_records_missing", "response_provenance")
    if not bool(phase.get("available")) or phase.get("status") != "available":
        raise MethodAUnavailable("phase_a_contract_unavailable", "phase_a_provenance")
    if bool(response.get("rf_restoration_applied", False)):
        raise MethodAUnavailable("response_not_noRF", "response_provenance")
    if not str(phase.get("contract_fingerprint") or ""):
        raise MethodAUnavailable("phase_a_contract_fingerprint_missing", "phase_a_provenance")
    if not str(phase.get("physical_pion_control_mask_fingerprint") or ""):
        raise MethodAUnavailable("physical_pion_control_mask_fingerprint_missing", "phase_a_provenance")
    pion_sources = ((response.get("source_provenance") or {}).get("pion") or {})
    if not pion_sources or "prompt" not in pion_sources:
        raise MethodAUnavailable("response_source_provenance_missing", "response_provenance")
    if any(
        entry.get("rf_state") != "noRF"
        or not str(entry.get("tree_name") or "").endswith("_noRF")
        or entry.get("proton_factor_scope") != "none"
        for entry in pion_sources.values()
    ):
        raise MethodAUnavailable("response_not_noRF", "response_provenance")


def _unavailable(reason, stage, config=None, t_edges=None, delta_edges=None, exception=None):
    result = {
        "schema_version": METHOD_A_SCHEMA_VERSION,
        "method": METHOD_A_METHOD,
        "status": "unavailable",
        "available": False,
        "reason": str(reason),
        "diagnostic_stage": str(stage),
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "rf_ct_required": False,
        "zerope_model_used": False,
        "t_edges": list(t_edges or ()),
        "delta_edges": list(delta_edges or ()),
        "configuration": _json_ready(config or {}),
        "fingerprint": None,
        "cells": [],
        "summary": {},
    }
    if exception is not None:
        result["exception_type"] = type(exception).__name__
        result["exception_message"] = str(exception)
    json.dumps(result, allow_nan=False)
    return result


def build_pion_hgcer_method_a(response_diagnostic, phase_a_contract, *, config=None):
    """Build the detached observed-positive-response Method-A diagnostic."""
    response = response_diagnostic if isinstance(response_diagnostic, dict) else {}
    phase = phase_a_contract if isinstance(phase_a_contract, dict) else {}
    resolved = None
    t_edges = []
    delta_edges = []
    try:
        resolved = _resolved_config(config)
        _validate_provenance(response, phase)
        t_edges = _strict_edges(phase.get("canonical_t_edges"), "canonical_t")
        delta_edges = _strict_edges(phase.get("delta_edges"), "delta")
        response_t_edges = _strict_edges(response.get("t_edges"), "response_t")
        response_delta_edges = _strict_edges(response.get("delta_edges"), "response_delta")
        if response_t_edges != t_edges:
            raise MethodAUnavailable("canonical_geometry_mismatch", "canonical_geometry")
        if response_delta_edges != delta_edges:
            raise MethodAUnavailable("delta_geometry_mismatch", "canonical_geometry")
        coordinate_fingerprint = str(phase.get("coordinate_fingerprint") or "")
        if not coordinate_fingerprint or str(response.get("coordinate_fingerprint") or "") != coordinate_fingerprint:
            raise MethodAUnavailable("coordinate_fingerprint_mismatch", "canonical_geometry")

        records = tuple(((response.get("records") or {}).get("pion") or ()))
        if not records:
            raise MethodAUnavailable("response_records_missing", "response_population")
        cells = {
            (t_index, delta_index): _empty_cell(
                t_index, delta_index, t_edges, delta_edges
            )
            for t_index in range(len(t_edges) - 1)
            for delta_index in range(len(delta_edges) - 1)
        }
        audit = {
            "records_seen": 0,
            "records_outside_t": 0,
            "records_outside_delta": 0,
            "records_nonpositive_response": 0,
            "records_not_nommcuts_or_allcuts": 0,
            "prompt_positive_nommcuts_records": 0,
            "signed_positive_nommcuts_records": 0,
        }
        fingerprint_records = []
        pion_source_provenance = (
            (response.get("source_provenance") or {}).get("pion") or {}
        )
        for source_record in records:
            if not isinstance(source_record, dict) and not isinstance(source_record, MappingProxyType):
                raise MethodAUnavailable("response_record_invalid", "record_validation")
            record = dict(source_record)
            audit["records_seen"] += 1
            if str(record.get("coordinate_fingerprint") or "") != coordinate_fingerprint:
                raise MethodAUnavailable("record_coordinate_fingerprint_mismatch", "record_validation")
            if bool(record.get("rf_applied_to_diagnostic", False)):
                raise MethodAUnavailable("record_rf_applied", "record_validation")
            if record.get("proton_cleaning_factor") is not None:
                raise MethodAUnavailable("pion_record_has_proton_factor", "record_validation")
            analysis_t = _finite(record.get("analysis_t"), "analysis_t")
            delta = _finite(record.get("ssdelta"), "delta")
            npe = _finite(record.get("P_hgcer_npeSum"), "response_value")
            x = _finite(record.get("P_hgcer_xAtCer"), "hgcer_x")
            y = _finite(record.get("P_hgcer_yAtCer"), "hgcer_y")
            stored_t = _stored_index(record.get("canonical_t_index"), "canonical_t")
            stored_delta = _stored_index(record.get("delta_index"), "delta")
            recomputed_t = _canonical_index(analysis_t, t_edges)
            recomputed_delta = _canonical_index(delta, delta_edges)
            if stored_t != recomputed_t:
                raise MethodAUnavailable("canonical_t_assignment_mismatch", "record_validation")
            if stored_delta != recomputed_delta:
                raise MethodAUnavailable("delta_assignment_mismatch", "record_validation")
            if recomputed_t is None:
                audit["records_outside_t"] += 1
                continue
            if recomputed_delta is None:
                audit["records_outside_delta"] += 1
                continue
            if npe <= 0.0:
                audit["records_nonpositive_response"] += 1
                continue
            allcuts = bool(record.get("allcuts", False))
            nommcuts = bool(record.get("nommcuts", False))
            if not allcuts and not nommcuts:
                audit["records_not_nommcuts_or_allcuts"] += 1
                continue
            source_label = str(record.get("source_label") or "")
            if source_label not in {"prompt", "rand", "dummy", "dummy_rand"}:
                raise MethodAUnavailable("response_source_label_invalid", "record_validation")
            source_entry = pion_source_provenance.get(source_label)
            if not isinstance(source_entry, dict):
                raise MethodAUnavailable("record_source_provenance_missing", "record_validation")
            weight = _finite(record.get("diagnostic_weight"), "diagnostic_weight")
            record_coefficient = _finite(record.get("coefficient"), "source_coefficient")
            provenance_coefficient = _finite(
                source_entry.get("coefficient"), "provenance_source_coefficient"
            )
            if not (
                math.isclose(weight, record_coefficient, rel_tol=1.0e-12, abs_tol=1.0e-12)
                and math.isclose(weight, provenance_coefficient, rel_tol=1.0e-12, abs_tol=1.0e-12)
            ):
                raise MethodAUnavailable("signed_source_coefficient_mismatch", "record_validation")
            expected_sign = {
                "prompt": 1.0,
                "rand": -1.0,
                "dummy": -1.0,
                "dummy_rand": 1.0,
            }[source_label]
            if weight == 0.0 or math.copysign(1.0, weight) != expected_sign:
                raise MethodAUnavailable("signed_source_algebra_mismatch", "record_validation")
            try:
                entry_index = int(record.get("entry_index", -1))
            except (TypeError, ValueError):
                raise MethodAUnavailable("entry_index_invalid", "record_validation")
            normalized = {
                "source_label": source_label,
                "entry_index": entry_index,
                "analysis_t": analysis_t,
                "t_index": recomputed_t,
                "delta": delta,
                "delta_index": recomputed_delta,
                "npe": npe,
                "x": x,
                "y": y,
                "allcuts": allcuts,
                "nommcuts": nommcuts,
                "weight": weight,
            }
            fingerprint_records.append(normalized)
            cell = cells[(recomputed_t, recomputed_delta)]
            if source_label == "prompt" and nommcuts:
                cell["prompt_nommcuts"].append(normalized)
                audit["prompt_positive_nommcuts_records"] += 1
            if source_label == "prompt" and allcuts:
                cell["prompt_allcuts"].append(normalized)
            if nommcuts:
                _add_signed(cell["signed"]["positive"], weight)
                _add_signed(
                    cell["signed"]["low" if npe <= 2.0 else "control"], weight
                )
                audit["signed_positive_nommcuts_records"] += 1

        if audit["prompt_positive_nommcuts_records"] <= 0:
            raise MethodAUnavailable("response_population_empty", "response_population")

        serialized_cells = []
        for key in sorted(cells):
            source = cells[key]
            prompt = _prompt_metrics(source["prompt_nommcuts"], resolved)
            allcuts = _prompt_metrics(source["prompt_allcuts"], resolved)
            signed = {
                name: _finish_signed(metrics)
                for name, metrics in source["signed"].items()
            }
            signed_ratio, signed_sigma = _signed_ratio(signed, resolved)
            cell = {
                "t_index": source["t_index"],
                "t_low": source["t_low"],
                "t_high": source["t_high"],
                "delta_index": source["delta_index"],
                "delta_low": source["delta_low"],
                "delta_high": source["delta_high"],
                "prompt_positive_count": prompt["positive"],
                "prompt_low_count": prompt["low"],
                "prompt_control_count": prompt["control"],
                "partition_closure_passed": prompt["partition_closure_passed"],
                "f_low": prompt["f_low"],
                "f_low_low": prompt["f_low_low"],
                "f_low_high": prompt["f_low_high"],
                "R_low_control": prompt["ratio"],
                "R_low_control_low": prompt["ratio_low"],
                "R_low_control_high": prompt["ratio_high"],
                "mean_npe": prompt["npe"]["mean"],
                "median_npe": prompt["npe"]["q50"],
                "npe_q10": prompt["npe"]["q10"],
                "npe_q25": prompt["npe"]["q25"],
                "npe_q50": prompt["npe"]["q50"],
                "npe_q75": prompt["npe"]["q75"],
                "npe_q90": prompt["npe"]["q90"],
                "minimum_npe": prompt["npe"]["minimum"],
                "maximum_npe": prompt["npe"]["maximum"],
                "hgcer_x_mean": prompt["x"]["mean"],
                "hgcer_x_rms": prompt["x"]["rms"],
                "hgcer_x_median": prompt["x"]["q50"],
                "hgcer_y_mean": prompt["y"]["mean"],
                "hgcer_y_rms": prompt["y"]["rms"],
                "hgcer_y_median": prompt["y"]["q50"],
                "signed_positive_yield": signed["positive"]["sumw"],
                "signed_low_yield": signed["low"]["sumw"],
                "signed_control_yield": signed["control"]["sumw"],
                "signed_positive_abs_support": signed["positive"]["sumabs"],
                "signed_low_abs_support": signed["low"]["sumabs"],
                "signed_control_abs_support": signed["control"]["sumabs"],
                "signed_positive_sumw2": signed["positive"]["sumw2"],
                "signed_low_sumw2": signed["low"]["sumw2"],
                "signed_control_sumw2": signed["control"]["sumw2"],
                "signed_positive_neff": signed["positive"]["neff"],
                "signed_low_neff": signed["low"]["neff"],
                "signed_control_neff": signed["control"]["neff"],
                "signed_R_low_control": signed_ratio,
                "signed_R_low_control_sigma": signed_sigma,
                "prompt_vs_signed_status": _prompt_vs_signed_status(
                    prompt, signed_ratio, signed_sigma, resolved
                ),
                "prompt_allcuts_positive_count": allcuts["positive"],
                "prompt_allcuts_low_count": allcuts["low"],
                "prompt_allcuts_control_count": allcuts["control"],
                "f_low_allcuts": allcuts["f_low"],
                "R_low_control_allcuts": allcuts["ratio"],
                "nommcuts_vs_allcuts_status": _nommcuts_vs_allcuts_status(
                    prompt, allcuts
                ),
                "support_class": prompt["support_class"],
                "method_A_status": "available" if prompt["available"] else "unavailable",
                "method_A_reason": None if prompt["available"] else "support_insufficient",
            }
            serialized_cells.append(cell)

        support_counts = {
            name: sum(1 for cell in serialized_cells if cell["support_class"] == name)
            for name in ("supported", "marginal", "unsupported")
        }
        response_mask = {
            "field": "P_hgcer_npeSum",
            "operator": ">",
            "value": 0.0,
            "selection": "prompt_noRF_nommcuts",
        }
        source_algebra = [
            {
                "source_label": label,
                "source_role": entry.get("source_role"),
                "coefficient": entry.get("coefficient"),
                "tree_name": entry.get("tree_name"),
                "rf_state": entry.get("rf_state"),
            }
            for label, entry in sorted(
                (((response.get("source_provenance") or {}).get("pion") or {}).items())
            )
        ]
        fingerprint_records.sort(
            key=lambda record: (
                record["source_label"], record["entry_index"], record["t_index"],
                record["delta_index"], record["npe"],
            )
        )
        event_population_fingerprint = _payload_hash(fingerprint_records)
        fingerprint_inputs = {
            "schema_version": METHOD_A_SCHEMA_VERSION,
            "method": METHOD_A_METHOD,
            "canonical_t_edges": t_edges,
            "delta_edges": delta_edges,
            "coordinate_fingerprint": coordinate_fingerprint,
            "phase_a_contract_fingerprint": phase.get("contract_fingerprint"),
            "physical_pion_control_mask_fingerprint": phase.get(
                "physical_pion_control_mask_fingerprint"
            ),
            "pion_pid_response_mask": response_mask,
            "pion_pid_response_mask_fingerprint": _payload_hash(response_mask),
            "source_algebra": source_algebra,
            "configuration": resolved,
            "event_population_fingerprint": event_population_fingerprint,
        }
        result = {
            "schema_version": METHOD_A_SCHEMA_VERSION,
            "method": METHOD_A_METHOD,
            "status": "available",
            "available": True,
            "reason": None,
            "diagnostic_stage": "complete",
            "non_authoritative": True,
            "production_objects_mutated": False,
            "refinement_applied": False,
            "rf_ct_required": False,
            "zerope_model_used": False,
            "t_edges": t_edges,
            "delta_edges": delta_edges,
            "coordinate_fingerprint": coordinate_fingerprint,
            "phase_a_contract_fingerprint": phase.get("contract_fingerprint"),
            "response_population_definition": "prompt_noRF_nommcuts_P_hgcer_npeSum_gt_0",
            "physical_control_definition": "P_hgcer_npeSum_gt_2",
            "low_response_definition": "0_lt_P_hgcer_npeSum_le_2",
            "pion_pid_response_mask_fingerprint": fingerprint_inputs[
                "pion_pid_response_mask_fingerprint"
            ],
            "configuration": resolved,
            "event_population_fingerprint": event_population_fingerprint,
            "fingerprint_inputs": fingerprint_inputs,
            "fingerprint": _payload_hash(fingerprint_inputs),
            "cells": serialized_cells,
            "summary": {
                **audit,
                "cell_count": len(serialized_cells),
                "support_counts": support_counts,
                "partition_closure_passed": all(
                    cell["partition_closure_passed"] for cell in serialized_cells
                ),
            },
        }
        result = _json_ready(result)
        json.dumps(result, allow_nan=False)
        return result
    except MethodAUnavailable as exc:
        return _unavailable(
            exc.reason, exc.stage, resolved, t_edges, delta_edges, exception=exc
        )
    except Exception as exc:
        return _unavailable(
            "unexpected_method_a_build_failure",
            "unexpected_exception",
            resolved,
            t_edges,
            delta_edges,
            exception=exc,
        )


def summarize_pion_hgcer_method_a(result):
    """Return a lightweight detached summary without duplicating cell payloads."""
    payload = result if isinstance(result, dict) else {}
    summary = {
        "schema_version": payload.get("schema_version", METHOD_A_SCHEMA_VERSION),
        "method": payload.get("method", METHOD_A_METHOD),
        "status": payload.get("status", "unavailable"),
        "available": bool(payload.get("available", False)),
        "reason": payload.get("reason"),
        "diagnostic_stage": payload.get("diagnostic_stage"),
        "fingerprint": payload.get("fingerprint"),
        "coordinate_fingerprint": payload.get("coordinate_fingerprint"),
        "phase_a_contract_fingerprint": payload.get("phase_a_contract_fingerprint"),
        "response_population_definition": payload.get("response_population_definition"),
        "physical_control_definition": payload.get("physical_control_definition"),
        "low_response_definition": payload.get("low_response_definition"),
        "t_edges": list(payload.get("t_edges") or ()),
        "delta_edges": list(payload.get("delta_edges") or ()),
        "cell_count": len(payload.get("cells") or ()),
        "support_counts": dict((payload.get("summary") or {}).get("support_counts") or {}),
        "prompt_positive_nommcuts_records": (payload.get("summary") or {}).get(
            "prompt_positive_nommcuts_records", 0
        ),
        "uncertainty_method": (payload.get("configuration") or {}).get(
            "uncertainty_method"
        ),
        "support_thresholds": dict(
            (payload.get("configuration") or {}).get("support") or {}
        ),
        "rf_ct_required": bool(payload.get("rf_ct_required", False)),
        "zerope_model_used": bool(payload.get("zerope_model_used", False)),
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
