#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2026-05-05 00:00:00 codex"
# ================================================================
#

import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

from ltsep import Misc
from ltsep import Root

lt = Root(os.path.realpath(__file__), "Plot_LTSep")

LTANAPATH = lt.LTANAPATH

sys.path.append("utility")
from utility import flatten_hist
from background_config import (
    MIN_PHI_BINS,
    MIN_T_BINS,
    PHI_BIN_MAX_DEG,
    PHI_BIN_MIN_DEG,
    PHI_BIN_MIN_EVENTS,
    T_BIN_ADJUST_MAX_ITERATIONS,
    T_BIN_ADJUST_TOLERANCE,
    T_BIN_EDGE_BIAS,
    T_BIN_MIN_EVENTS,
    get_proton_contamination_cleaning_config,
)

PRE_PARTICLE_SUBTRACTION_T_HIST_KEY = "_binning_H_t_DATA_pre_particle_subtraction"
PRE_PARTICLE_SUBTRACTION_PHI_HIST_KEY = "_binning_H_ph_q_DATA_pre_particle_subtraction"


def _clone_histlist(histlist):
    histlist_copy = []
    for hist in histlist:
        hist_copy = {}
        for key, val in hist.items():
            if hasattr(val, "Clone") and callable(getattr(val, "Clone")):
                hist_copy[key] = val.Clone()
            else:
                hist_copy[key] = val
        histlist_copy.append(hist_copy)
    return histlist_copy


def _get_binning_histograms(hist):
    """Prefer kaon distributions captured before slow-proton and pion subtraction."""
    t_hist = hist.get(PRE_PARTICLE_SUBTRACTION_T_HIST_KEY)
    phi_hist = hist.get(PRE_PARTICLE_SUBTRACTION_PHI_HIST_KEY)
    if t_hist is None:
        t_hist = hist["H_t_DATA"]
    if phi_hist is None:
        phi_hist = hist["H_ph_q_DATA"]
    return t_hist, phi_hist


def _collect_bin_samples(histlist, inpDict):
    tmin = float(inpDict["tmin"])
    tmax = float(inpDict["tmax"])

    h_t_right = np.array([])
    h_t_left = np.array([])
    h_t_center = np.array([])
    h_phi_right = np.array([])
    h_phi_left = np.array([])
    h_phi_center = np.array([])

    for hist in _clone_histlist(histlist):
        normfac_data = hist.get("normfac_data", inpDict["normfac_data"])
        t_hist, phi_hist = _get_binning_histograms(hist)
        t_hist.Scale(1.0 / normfac_data)
        phi_hist.Scale(1.0 / normfac_data)

        t_values = flatten_hist(t_hist)
        phi_deg = [(phi * (180.0 / math.pi)) for phi in flatten_hist(phi_hist)]

        if hist["phi_setting"] == "Right":
            h_t_right = np.append(h_t_right, t_values)
            h_phi_right = np.append(h_phi_right, phi_deg)
        elif hist["phi_setting"] == "Left":
            h_t_left = np.append(h_t_left, t_values)
            h_phi_left = np.append(h_phi_left, phi_deg)
        elif hist["phi_setting"] == "Center":
            h_t_center = np.append(h_t_center, t_values)
            h_phi_center = np.append(h_phi_center, phi_deg)

    h_t_combined = np.concatenate((h_t_right, h_t_left, h_t_center))
    if h_t_combined.size == 0:
        raise ValueError("No t-distribution data available for bin finding.")

    min_index = np.argmin(h_t_combined)
    max_index = np.argmax(h_t_combined)
    h_t_combined[min_index] = tmin
    h_t_combined[max_index] = tmax

    h_phi_combined = np.concatenate((h_phi_right, h_phi_left, h_phi_center))
    if h_phi_combined.size == 0:
        raise ValueError("No phi-distribution data available for bin finding.")

    return h_t_combined, h_phi_combined


def _find_phi_bins(phi_values, requested_num_phi_bins, quiet=False):
    if not quiet:
        print("\nFinding phi bins...")

    num_phi_bins = max(MIN_PHI_BINS, int(requested_num_phi_bins))
    bin_edges = _build_uniform_phi_bins(num_phi_bins)

    while True:
        counts, bins = np.histogram(phi_values, bin_edges)
        bad_bins = np.where(counts < PHI_BIN_MIN_EVENTS)[0]
        if len(bad_bins) == 0:
            break

        num_phi_bins -= 1
        if num_phi_bins < MIN_PHI_BINS:
            raise ValueError(
                "Only {} phi-bins achieve a minimum of {} events.".format(
                    num_phi_bins, PHI_BIN_MIN_EVENTS
                )
            )
        bin_edges = _build_uniform_phi_bins(num_phi_bins)

    if not quiet:
        print(
            "Using fixed equal-width phi bins from {:.1f} to {:.1f} degrees".format(
                float(PHI_BIN_MIN_DEG),
                float(PHI_BIN_MAX_DEG),
            )
        )
        for i, val in enumerate(counts):
            print(
                "Bin {} from {:.1f} to {:.1f} has {} events".format(
                    i + 1, bins[i], bins[i + 1], val
                )
            )
        print("phi_bins = ", bins)

    return bins, counts


def _build_uniform_phi_bins(num_phi_bins):
    return np.linspace(float(PHI_BIN_MIN_DEG), float(PHI_BIN_MAX_DEG), int(num_phi_bins) + 1)


def _calculate_edge_scaling(index, total_bins, edge_bias):
    normalized_index = index / (total_bins - 1)
    return (1 - abs(2 * normalized_index - 1)) ** edge_bias


def _adjust_t_bins(values, requested_num_t_bins, quiet=False):
    nbin = max(MIN_T_BINS, int(requested_num_t_bins)) + 1
    tmin = np.min(values)
    tmax = np.max(values)
    bin_edges = np.linspace(tmin, tmax, num=nbin)
    counts, _ = np.histogram(values, bins=bin_edges)

    iteration = 0
    while np.any(counts < T_BIN_MIN_EVENTS) and iteration < T_BIN_ADJUST_MAX_ITERATIONS:
        new_bin_edges = bin_edges.copy()
        for i in range(1, len(bin_edges) - 1):
            edge_scaling = _calculate_edge_scaling(i, len(bin_edges), T_BIN_EDGE_BIAS)
            if counts[i - 1] < T_BIN_MIN_EVENTS:
                new_bin_edges[i] += T_BIN_ADJUST_TOLERANCE * edge_scaling
            elif counts[i] < T_BIN_MIN_EVENTS:
                new_bin_edges[i] -= T_BIN_ADJUST_TOLERANCE * edge_scaling

        bin_edges = new_bin_edges
        counts, _ = np.histogram(values, bins=bin_edges)
        iteration += 1

        if not quiet:
            Misc.progressBar(iteration, T_BIN_ADJUST_MAX_ITERATIONS - 1, bar_length=25)

    if np.any(counts < T_BIN_MIN_EVENTS) and not quiet:
        print(
            "WARNING: Could not achieve minimum {} events in all bins after {} iterations.".format(
                T_BIN_MIN_EVENTS, T_BIN_ADJUST_MAX_ITERATIONS
            )
        )

    return bin_edges, counts


def _find_t_bins(t_values, requested_num_t_bins, quiet=False):
    if not quiet:
        print("\nFinding t bins...")

    try:
        bin_edges = _adjust_t_bins(t_values, requested_num_t_bins, quiet=quiet)
        bins, counts = bin_edges
        counts, bins = np.histogram(t_values, bins)
    except ValueError:
        raise ValueError("Unavoidable empty bins. Tighten t-range or adjust number of t-bins.")

    if np.size(counts) == 0:
        raise ValueError("t-binning failed: no valid bins available.")

    actual_num_t_bins = len(counts)
    if actual_num_t_bins < MIN_T_BINS:
        raise ValueError(
            "Only {} t-bin achieved a minimum of {} events.".format(
                actual_num_t_bins, T_BIN_MIN_EVENTS
            )
        )

    if not quiet:
        for i, val in enumerate(counts):
            print(
                "Bin {} from {:.3f} to {:.3f} has {} events".format(
                    i + 1, bins[i], bins[i + 1], val
                )
            )
        print("t_bins = ", bins)

    return bins, counts


def propose_bins(histlist, inpDict, num_t_bins=None, num_phi_bins=None, quiet=False):
    requested_t_bins = int(inpDict["NumtBins"] if num_t_bins is None else num_t_bins)
    requested_phi_bins = int(inpDict["NumPhiBins"] if num_phi_bins is None else num_phi_bins)

    t_values, phi_values = _collect_bin_samples(histlist, inpDict)
    phi_bins, phi_counts = _find_phi_bins(phi_values, requested_phi_bins, quiet=quiet)
    t_bins, t_counts = _find_t_bins(t_values, requested_t_bins, quiet=quiet)

    return {
        "t_bins": np.array(t_bins, dtype=float),
        "phi_bins": np.array(phi_bins, dtype=float),
        "requested_num_t_bins": requested_t_bins,
        "requested_num_phi_bins": requested_phi_bins,
        "actual_num_t_bins": len(t_bins) - 1,
        "actual_num_phi_bins": len(phi_bins) - 1,
        "t_counts": np.array(t_counts, dtype=int),
        "phi_counts": np.array(phi_counts, dtype=int),
    }


def apply_bin_proposal(inpDict, proposal):
    inpDict["NumtBins"] = int(proposal["actual_num_t_bins"])
    inpDict["NumPhiBins"] = int(proposal["actual_num_phi_bins"])
    inpDict["t_bins"] = np.array(proposal["t_bins"], dtype=float)
    inpDict["phi_bins"] = np.array(proposal["phi_bins"], dtype=float)
    return inpDict


def _canonical_t_binning_config(inp_dict):
    config = get_proton_contamination_cleaning_config(inp_dict=inp_dict)
    return dict(config.get("t_binning") or {})


def _canonical_epsset(inp_dict):
    value = str(inp_dict.get("EPSSET", "low")).strip().lower()
    return value[:-1] if value.endswith("e") else value


def _interval_paths(inp_dict):
    particle_type = inp_dict["ParticleType"]
    q2_token = inp_dict["Q2"].replace("p", "")
    w_token = inp_dict["W"].replace("p", "")
    phi_path = "{}/src/{}/phi_bin_interval_Q{}W{}".format(
        LTANAPATH, particle_type, q2_token, w_token
    )
    t_path = "{}/src/{}/t_bin_interval_Q{}W{}".format(
        LTANAPATH, particle_type, q2_token, w_token
    )
    return {
        "phi_path": phi_path,
        "t_path": t_path,
        "phi_metadata_path": "{}.metadata.json".format(phi_path),
        "t_metadata_path": "{}.metadata.json".format(t_path),
    }


def _json_ready(value):
    if isinstance(value, np.ndarray):
        return [_json_ready(item) for item in value.tolist()]
    if isinstance(value, (np.floating, float)):
        value = float(value)
        if not math.isfinite(value):
            raise ValueError("canonical-binning JSON cannot contain nonfinite values")
        return value
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    return value


def _canonical_binning_semantic_config(inp_dict, t_binning_config):
    """Return only stable choices that change the canonical t proposal."""
    return {
        "requested_num_t_bins": int(inp_dict["NumtBins"]),
        "tmin": float(inp_dict["tmin"]),
        "tmax": float(inp_dict["tmax"]),
        "minimum_t_bins": int(MIN_T_BINS),
        "minimum_t_events": int(T_BIN_MIN_EVENTS),
        "adjust_tolerance": float(T_BIN_ADJUST_TOLERANCE),
        "adjust_max_iterations": int(T_BIN_ADJUST_MAX_ITERATIONS),
        "edge_bias": float(T_BIN_EDGE_BIAS),
        "canonical_bin_support_metric": str(
            t_binning_config.get("canonical_bin_support_metric", "raw_event_count")
        ),
        "canonical_bin_support_threshold": float(
            t_binning_config.get("canonical_bin_support_threshold", T_BIN_MIN_EVENTS)
        ),
        "algorithm_identifier": "find_bins_adjust_t_bins",
        "algorithm_version": 1,
    }


def _canonical_binning_config_hash(semantic_config):
    encoded = json.dumps(
        _json_ready(semantic_config), sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _read_text_interval_edges(interval_path):
    with open(interval_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    if len(lines) < 2:
        raise ValueError("interval file has no edge line")
    tokens = [token.strip() for token in lines[1].split("\t") if token.strip()]
    if len(tokens) < 2:
        raise ValueError("interval file has fewer than two edges")
    return np.asarray([float(token) for token in tokens], dtype=float)


def _prepass_records(pre_subtraction_histograms):
    records = []
    for payload in pre_subtraction_histograms or []:
        for record in (payload or {}).get("records") or []:
            try:
                t_value = float(record["adj_t"])
                phi_value = float(record["phi_value"])
                coefficient = float(record["physical_coefficient"])
            except (KeyError, TypeError, ValueError):
                continue
            if not (math.isfinite(t_value) and math.isfinite(phi_value) and math.isfinite(coefficient)):
                continue
            records.append(
                {
                    "adj_t": t_value,
                    "phi_value": phi_value,
                    "physical_coefficient": coefficient,
                    "source_label": str(record.get("source_label", "")),
                    "entry_index": int(record.get("entry_index", -1)),
                }
            )
    return records


def _support_summaries(t_values, physical_weights, t_edges):
    raw_counts, _ = np.histogram(t_values, bins=t_edges)
    signed, _ = np.histogram(t_values, bins=t_edges, weights=physical_weights)
    absolute, _ = np.histogram(t_values, bins=t_edges, weights=np.abs(physical_weights))
    positive, _ = np.histogram(
        t_values, bins=t_edges, weights=np.maximum(physical_weights, 0.0)
    )
    negative, _ = np.histogram(
        t_values, bins=t_edges, weights=np.abs(np.minimum(physical_weights, 0.0))
    )
    return {
        "raw_event_count_by_t_bin": np.asarray(raw_counts, dtype=int),
        "signed_weighted_yield_by_t_bin": np.asarray(signed, dtype=float),
        "absolute_weighted_support_by_t_bin": np.asarray(absolute, dtype=float),
        "positive_weighted_support_by_t_bin": np.asarray(positive, dtype=float),
        "negative_weighted_support_by_t_bin": np.asarray(negative, dtype=float),
    }


def _find_t_bins_from_prepass(
    t_values,
    physical_weights,
    requested_num_t_bins,
    support_metric,
    support_threshold=T_BIN_MIN_EVENTS,
    quiet=False,
):
    allowed_metrics = {
        "raw_event_count": None,
        "absolute_weighted_support": np.abs(physical_weights),
        "positive_weighted_support": np.maximum(physical_weights, 0.0),
    }
    if support_metric == "signed_weighted_yield":
        raise ValueError("signed_weighted_yield cannot be a canonical support threshold metric")
    if support_metric not in allowed_metrics:
        raise ValueError("unsupported canonical t support metric '{}'".format(support_metric))
    if t_values.size == 0:
        raise ValueError("No selected pre-particle-subtraction t records available.")

    proposal_values = np.asarray(t_values, dtype=float).copy()
    proposal_values[np.argmin(proposal_values)] = float(np.min(t_values))
    proposal_values[np.argmax(proposal_values)] = float(np.max(t_values))
    support_weights = allowed_metrics[support_metric]
    n_edges = max(MIN_T_BINS, int(requested_num_t_bins)) + 1
    edges = np.linspace(float(np.min(proposal_values)), float(np.max(proposal_values)), n_edges)
    support, _ = np.histogram(proposal_values, bins=edges, weights=support_weights)
    iteration = 0
    support_threshold = float(support_threshold)
    if not math.isfinite(support_threshold) or support_threshold < 0.0:
        raise ValueError("canonical support threshold must be finite and nonnegative")
    while np.any(support < support_threshold) and iteration < T_BIN_ADJUST_MAX_ITERATIONS:
        adjusted = edges.copy()
        for index in range(1, len(edges) - 1):
            scaling = _calculate_edge_scaling(index, len(edges), T_BIN_EDGE_BIAS)
            if support[index - 1] < support_threshold:
                adjusted[index] += T_BIN_ADJUST_TOLERANCE * scaling
            elif support[index] < support_threshold:
                adjusted[index] -= T_BIN_ADJUST_TOLERANCE * scaling
        edges = adjusted
        support, _ = np.histogram(proposal_values, bins=edges, weights=support_weights)
        iteration += 1
        if not quiet:
            Misc.progressBar(iteration, T_BIN_ADJUST_MAX_ITERATIONS - 1, bar_length=25)
    return np.asarray(edges, dtype=float), np.asarray(support, dtype=float)


def _validate_authoritative_t_interval(inp_dict, consumer_epsilon):
    """Validate text and metadata together; return all rejection reasons."""
    paths = _interval_paths(inp_dict)
    t_binning = _canonical_t_binning_config(inp_dict)
    reasons = []
    metadata = None
    text_edges = None
    if not os.path.exists(paths["t_path"]):
        reasons.append("interval_file_missing")
    else:
        try:
            text_edges = _read_text_interval_edges(paths["t_path"])
        except Exception as exc:
            reasons.append("interval_file_unreadable:{}".format(exc))
    if not os.path.exists(paths["t_metadata_path"]):
        reasons.append("metadata_sidecar_missing")
    else:
        try:
            with open(paths["t_metadata_path"], "r", encoding="utf-8") as handle:
                metadata = json.load(handle)
        except Exception as exc:
            reasons.append("metadata_sidecar_unreadable:{}".format(exc))

    if metadata is not None:
        expected_config = _canonical_binning_semantic_config(inp_dict, t_binning)
        expected_hash = _canonical_binning_config_hash(expected_config)
        expected_eps = "low"
        checks = (
            (metadata.get("schema_version") == int(t_binning.get("metadata_schema_version", 1)), "metadata_schema_version_mismatch"),
            (metadata.get("particle_type") == inp_dict.get("ParticleType"), "particle_type_mismatch"),
            (metadata.get("Q2_token") == inp_dict.get("Q2"), "q2_token_mismatch"),
            (metadata.get("W_token") == inp_dict.get("W"), "w_token_mismatch"),
            (
                metadata.get("kinematic_key")
                == "Q{}W{}".format(inp_dict.get("Q2"), inp_dict.get("W")),
                "kinematic_key_mismatch",
            ),
            (metadata.get("source_stage") == "pre_particle_subtraction", "source_stage_mismatch"),
            (metadata.get("source_epsilon") == expected_eps, "source_epsilon_mismatch"),
            (bool(metadata.get("shared_with_high_epsilon", False)), "low_to_high_sharing_not_enabled"),
            (metadata.get("requested_num_t_bins") == int(inp_dict["NumtBins"]), "requested_bin_count_mismatch"),
            (metadata.get("actual_num_t_bins") == int(inp_dict["NumtBins"]), "actual_bin_count_mismatch"),
            (metadata.get("binning_config_hash") == expected_hash, "binning_config_hash_mismatch"),
            (metadata.get("algorithm_identifier") in set(t_binning.get("allowed_algorithm_identifiers") or ()), "algorithm_identifier_not_allowed"),
            (metadata.get("algorithm_version") in set(t_binning.get("allowed_algorithm_versions") or ()), "algorithm_version_not_allowed"),
        )
        for valid, reason in checks:
            if not valid:
                reasons.append(reason)
        try:
            metadata_edges = np.asarray(metadata.get("t_edges") or [], dtype=float)
        except (TypeError, ValueError):
            metadata_edges = np.asarray([], dtype=float)
            reasons.append("metadata_edges_invalid")
        tolerance = float(t_binning.get("edge_tolerance", 1.0e-9))
        if metadata_edges.size != int(inp_dict["NumtBins"]) + 1:
            reasons.append("metadata_edge_count_mismatch")
        if metadata_edges.size and metadata.get("actual_num_t_bins") != int(metadata_edges.size - 1):
            reasons.append("metadata_actual_bin_count_edge_count_mismatch")
        if metadata_edges.size and (not np.all(np.isfinite(metadata_edges)) or np.any(np.diff(metadata_edges) <= 0.0)):
            reasons.append("metadata_edges_not_strictly_increasing_finite")
        if metadata_edges.size and (
            abs(float(metadata_edges[0]) - float(inp_dict["tmin"])) > tolerance
            or abs(float(metadata_edges[-1]) - float(inp_dict["tmax"])) > tolerance
        ):
            reasons.append("metadata_edges_do_not_cover_configured_range")
        if text_edges is not None:
            if text_edges.size < 2:
                reasons.append("text_edge_count_mismatch")
            elif not np.all(np.isfinite(text_edges)) or np.any(np.diff(text_edges) <= 0.0):
                reasons.append("text_edges_not_strictly_increasing_finite")
            if text_edges.size != metadata_edges.size:
                reasons.append("text_metadata_edge_count_mismatch")
            elif not np.allclose(text_edges, metadata_edges, rtol=0.0, atol=tolerance):
                reasons.append("text_metadata_edges_mismatch")
    status = "validated_authoritative_interval_file" if not reasons else "metadata_mismatch"
    if reasons == ["metadata_sidecar_missing"]:
        status = "legacy_interval_rejected"
    return {
        "valid": not reasons,
        "validation_status": status,
        "validation_rejection_reasons": reasons,
        "interval_file": paths["t_path"],
        "metadata_file": paths["t_metadata_path"],
        "metadata": metadata,
        "t_edges": text_edges,
    }


def validate_canonical_t_edges(inp_dict, t_edges, context=""):
    canonical = dict(inp_dict.get("canonical_t_binning") or {})
    expected_source = canonical.get("t_edges")
    if expected_source is None:
        expected_source = inp_dict.get("t_bins")
    expected = np.asarray(expected_source if expected_source is not None else [], dtype=float)
    candidate = np.asarray(t_edges if t_edges is not None else [], dtype=float)
    tolerance = float(_canonical_t_binning_config(inp_dict).get("edge_tolerance", 1.0e-9))
    match = bool(expected.size == candidate.size and np.allclose(expected, candidate, rtol=0.0, atol=tolerance))
    max_difference = float(np.max(np.abs(expected - candidate))) if expected.size == candidate.size and expected.size else float("inf")
    diagnostics = {
        "context": str(context),
        "canonical_edge_match": match,
        "maximum_edge_difference": max_difference,
    }
    if not match and bool(get_proton_contamination_cleaning_config(inp_dict=inp_dict).get("strict_mode", False)):
        raise RuntimeError("canonical t edges differ in {}".format(context or "unknown context"))
    return diagnostics


def resolve_canonical_analysis_bins_pre_subtraction(
    pre_subtraction_histograms,
    inp_dict,
    *,
    allow_interval_file=True,
    write_interval_files=True,
    quiet=False,
):
    """Resolve immutable canonical analysis bins before proton cleaning.

    ``pre_subtraction_histograms`` is a sequence of lightweight prepass
    payloads.  It deliberately carries records and physical coefficients so
    support is independent of plotting-histogram normalization.
    """
    records = _prepass_records(pre_subtraction_histograms)
    if not records:
        raise ValueError("No selected pre-particle-subtraction records available for canonical bins.")
    t_values = np.asarray([record["adj_t"] for record in records], dtype=float)
    phi_values = np.asarray([record["phi_value"] for record in records], dtype=float)
    physical_weights = np.asarray(
        [record["physical_coefficient"] for record in records], dtype=float
    )
    t_binning = _canonical_t_binning_config(inp_dict)
    support_metric = str(t_binning.get("canonical_bin_support_metric", "raw_event_count"))
    allowed_metrics = set(t_binning.get("allowed_support_metrics") or ())
    if support_metric == "signed_weighted_yield":
        raise ValueError("signed_weighted_yield cannot be a canonical support threshold metric")
    if support_metric not in allowed_metrics:
        raise ValueError("canonical_bin_support_metric '{}' is not allowed".format(support_metric))

    requested_t = int(inp_dict["NumtBins"])
    requested_phi = int(inp_dict["NumPhiBins"])
    consumer_epsilon = _canonical_epsset(inp_dict)
    interval_validation = {
        "valid": False,
        "validation_status": "no_valid_source",
        "validation_rejection_reasons": ["authoritative_interval_disabled"],
    }
    if allow_interval_file and bool(t_binning.get("allow_authoritative_interval_file", True)):
        interval_validation = _validate_authoritative_t_interval(inp_dict, consumer_epsilon)

    if interval_validation.get("valid"):
        t_edges = np.asarray(interval_validation["t_edges"], dtype=float)
        source = "validated_authoritative_interval_file"
        validation_status = source
        metadata = dict(interval_validation.get("metadata") or {})
    elif consumer_epsilon == "high":
        # High epsilon is a consumer of low-epsilon canonical bins, never a
        # second independent producer for the same kinematic key.
        canonical = {
            "source": "no_valid_source",
            "validation_status": interval_validation.get("validation_status", "no_valid_source"),
            "validation_rejection_reasons": list(interval_validation.get("validation_rejection_reasons") or []),
            "consumer_epsilon": "high",
            "shared_from_epsilon": "low",
            "created_before_proton_cleaning": True,
        }
        inp_dict["canonical_t_binning"] = canonical
        raise RuntimeError("No compatible low-epsilon canonical t interval is available for high epsilon.")
    else:
        proposal_values = t_values.copy()
        proposal_values[np.argmin(proposal_values)] = float(inp_dict["tmin"])
        proposal_values[np.argmax(proposal_values)] = float(inp_dict["tmax"])
        t_edges, _ = _find_t_bins_from_prepass(
            proposal_values,
            physical_weights,
            requested_t,
            support_metric,
            t_binning.get("canonical_bin_support_threshold", T_BIN_MIN_EVENTS),
            quiet=quiet,
        )
        source = "computed_from_pre_particle_subtraction_histograms"
        validation_status = "computed_pre_particle_subtraction"
        metadata = {}

    phi_degrees = phi_values * (180.0 / math.pi)
    phi_edges, phi_counts = _find_phi_bins(phi_degrees, requested_phi, quiet=quiet)
    support = _support_summaries(t_values, physical_weights, t_edges)
    semantic_config = _canonical_binning_semantic_config(inp_dict, t_binning)
    paths = _interval_paths(inp_dict)
    source_stats = [
        _json_ready((payload or {}).get("source_stats") or {})
        for payload in pre_subtraction_histograms or []
    ]
    canonical = {
        "schema_version": int(t_binning.get("metadata_schema_version", 1)),
        "particle_type": inp_dict["ParticleType"],
        "Q2_token": inp_dict["Q2"],
        "W_token": inp_dict["W"],
        "kinematic_key": "Q{}W{}".format(inp_dict["Q2"], inp_dict["W"]),
        "source_epsilon": "low",
        "consumer_epsilon": consumer_epsilon,
        "shared_with_high_epsilon": True,
        "shared_from_epsilon": "low" if consumer_epsilon == "high" else None,
        "source_stage": "pre_particle_subtraction",
        "source": source,
        "resolution_status": validation_status,
        "validation_status": validation_status,
        "validation_rejection_reasons": list(interval_validation.get("validation_rejection_reasons") or []),
        "interval_file": paths["t_path"],
        "metadata_file": paths["t_metadata_path"],
        "requested_num_t_bins": requested_t,
        "actual_num_t_bins": len(t_edges) - 1,
        "tmin": float(inp_dict["tmin"]),
        "tmax": float(inp_dict["tmax"]),
        "t_edges": _json_ready(t_edges),
        "raw_event_count_by_t_bin": _json_ready(support["raw_event_count_by_t_bin"]),
        "signed_weighted_yield_by_t_bin": _json_ready(support["signed_weighted_yield_by_t_bin"]),
        "absolute_weighted_support_by_t_bin": _json_ready(support["absolute_weighted_support_by_t_bin"]),
        "positive_weighted_support_by_t_bin": _json_ready(support["positive_weighted_support_by_t_bin"]),
        "negative_weighted_support_by_t_bin": _json_ready(support["negative_weighted_support_by_t_bin"]),
        "support_metric": support_metric,
        "support_threshold": float(
            t_binning.get("canonical_bin_support_threshold", T_BIN_MIN_EVENTS)
        ),
        "binning_config": _json_ready(semantic_config),
        "binning_config_hash": _canonical_binning_config_hash(semantic_config),
        "algorithm_identifier": "find_bins_adjust_t_bins",
        "algorithm_version": 1,
        "input_provenance": {"source_stats": source_stats, "record_count": int(len(records))},
        "created_before_proton_cleaning": True,
        "creation_timestamp": datetime.now(timezone.utc).isoformat(),
        "edge_match": True,
    }
    # Preserve source-side provenance when a valid authoritative file wins,
    # while exposing current-run support in the active resolution payload.
    if metadata:
        canonical["authoritative_metadata"] = _json_ready(metadata)
    inp_dict["t_bins"] = np.asarray(t_edges, dtype=float)
    inp_dict["phi_bins"] = np.asarray(phi_edges, dtype=float)
    inp_dict["NumtBins"] = int(len(t_edges) - 1)
    inp_dict["NumPhiBins"] = int(len(phi_edges) - 1)
    inp_dict["canonical_t_binning"] = canonical
    artifacts = _interval_paths(inp_dict)
    if source == "computed_from_pre_particle_subtraction_histograms" and write_interval_files:
        artifacts = write_bin_interval_files(inp_dict, t_edges, phi_edges, canonical_metadata=canonical)
    canonical["artifacts"] = artifacts
    return {
        "t_bins": np.asarray(t_edges, dtype=float),
        "phi_bins": np.asarray(phi_edges, dtype=float),
        "t_counts": np.asarray(support["raw_event_count_by_t_bin"], dtype=int),
        "phi_counts": np.asarray(phi_counts, dtype=int),
        "requested_num_t_bins": requested_t,
        "actual_num_t_bins": len(t_edges) - 1,
        "source": source,
        "metadata": canonical,
        "support": support,
    }


def write_bin_interval_files(inpDict, t_bins, phi_bins, canonical_metadata=None):
    paths = _interval_paths(inpDict)
    phi_path = paths["phi_path"]
    t_path = paths["t_path"]

    phi_lines = ["\t{:.17g}".format(float(phi)) for phi in phi_bins]
    with open(phi_path, "w") as file:
        file.write(
            "{}\t{}\t{}\t{}\n".format(
                inpDict["Q2"].replace("p", "."),
                inpDict["W"].replace("p", "."),
                len(t_bins) - 1,
                len(phi_bins) - 1,
            )
        )
        file.writelines(phi_lines)

    t_lines = ["\t{:.17g}".format(float(t_val)) for t_val in t_bins]
    with open(t_path, "w") as file:
        file.write(
            "{}\t{}\t{}\t{}\n".format(
                inpDict["Q2"].replace("p", "."),
                inpDict["W"].replace("p", "."),
                len(t_bins) - 1,
                len(phi_bins) - 1,
            )
        )
        file.writelines(t_lines)

    artifacts = dict(paths)
    if canonical_metadata is not None:
        metadata = dict(canonical_metadata)
        metadata["interval_file"] = t_path
        metadata["metadata_file"] = paths["t_metadata_path"]
        metadata["t_edges"] = _json_ready(np.asarray(t_bins, dtype=float))
        metadata["actual_num_t_bins"] = len(t_bins) - 1
        metadata["created_before_proton_cleaning"] = True
        with open(paths["t_metadata_path"], "w", encoding="utf-8") as handle:
            json.dump(_json_ready(metadata), handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
        phi_metadata = {
            "schema_version": metadata.get("schema_version", 1),
            "particle_type": metadata.get("particle_type"),
            "Q2_token": metadata.get("Q2_token"),
            "W_token": metadata.get("W_token"),
            "source_epsilon": metadata.get("source_epsilon"),
            "source_stage": metadata.get("source_stage"),
            "interval_file": phi_path,
            "metadata_file": paths["phi_metadata_path"],
            "phi_edges": _json_ready(np.asarray(phi_bins, dtype=float)),
            "requested_num_phi_bins": len(phi_bins) - 1,
            "actual_num_phi_bins": len(phi_bins) - 1,
            "created_before_proton_cleaning": True,
            "creation_timestamp": metadata.get("creation_timestamp"),
        }
        with open(paths["phi_metadata_path"], "w", encoding="utf-8") as handle:
            json.dump(_json_ready(phi_metadata), handle, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return artifacts


def find_bins(histlist, inpDict):
    proposal = propose_bins(histlist, inpDict, quiet=False)
    apply_bin_proposal(inpDict, proposal)
    write_bin_interval_files(inpDict, proposal["t_bins"], proposal["phi_bins"])
    return proposal


def check_bins(histlist, inpDict):
    t_values, phi_values = _collect_bin_samples(histlist, inpDict)
    t_bins = histlist[0]["t_bins"]
    phi_bins = histlist[0]["phi_bins"]
    canonical_active = bool(inpDict.get("canonical_t_binning"))
    if canonical_active:
        edge_diagnostics = validate_canonical_t_edges(
            inpDict, t_bins, context="find_bins.check_bins"
        )
        inpDict["canonical_t_binning"]["check_bins_edge_validation"] = edge_diagnostics

    print("\nFinding phi bins...")
    phi_counts, phi_edges = np.histogram(phi_values, phi_bins)
    if len(phi_bins) - 1 != inpDict["NumPhiBins"]:
        print(
            "Number of phi-bins changed from {} to: {}".format(
                inpDict["NumPhiBins"], len(phi_counts)
            )
        )
        if not canonical_active:
            inpDict["NumPhiBins"] = len(phi_counts)
    for i, val in enumerate(phi_counts):
        print(
            "Bin {} from {:.1f} to {:.1f} has {} events".format(
                i + 1, phi_edges[i], phi_edges[i + 1], val
            )
        )
    print("phi_bins = ", phi_edges)

    print("\nFinding t bins...")
    t_counts, t_edges = np.histogram(t_values, t_bins)
    if len(t_bins) - 1 != inpDict["NumtBins"]:
        print(
            "Number of t-bins changed from {} to: {}".format(
                inpDict["NumtBins"], len(t_counts)
            )
        )
        if not canonical_active:
            inpDict["NumtBins"] = len(t_counts)
    for i, val in enumerate(t_counts):
        print(
            "Bin {} from {:.3f} to {:.3f} has {} events".format(
                i + 1, t_edges[i], t_edges[i + 1], val
            )
        )
    print("t_bins = ", t_edges)
    if canonical_active:
        inpDict["canonical_t_binning"]["check_bins_counts"] = {
            "t_counts": _json_ready(np.asarray(t_counts, dtype=int)),
            "phi_counts": _json_ready(np.asarray(phi_counts, dtype=int)),
        }
    return {"t_counts": t_counts, "phi_counts": phi_counts}
