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
import tempfile
import uuid
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


def _find_phi_bins(
    phi_values,
    requested_num_phi_bins,
    quiet=False,
    return_diagnostics=False,
    *,
    auto_rebin_phi=True,
    min_phi_bins=None,
    support_evaluator=None,
):
    """Resolve equal-width phi bins without losing the original request.

    The canonical sidecar must distinguish the configured request from the
    smaller grid accepted after the minimum-event reduction loop.
    """
    if not quiet:
        print("\nFinding phi bins...")

    original_requested = int(requested_num_phi_bins)
    if original_requested < 1:
        raise ValueError("requested_num_phi_bins must be positive")
    configured_minimum = MIN_PHI_BINS if min_phi_bins is None else int(min_phi_bins)
    # A configured minimum limits automatic *reduction*; it must never enlarge
    # an explicitly requested lower-resolution grid.
    effective_minimum = min(original_requested, max(1, configured_minimum))
    num_phi_bins = original_requested
    initial_num_phi_bins = int(num_phi_bins)
    reduction_iterations = 0
    bin_edges = _build_uniform_phi_bins(num_phi_bins)

    while True:
        counts, bins = np.histogram(phi_values, bin_edges)
        support = (
            dict(support_evaluator(bin_edges) or {})
            if support_evaluator is not None
            else {
                "passed": bool(np.all(counts >= PHI_BIN_MIN_EVENTS)),
                "final_phi_event_counts": [int(value) for value in counts],
                "support_policy": "all_phi_bins",
                "support_threshold": int(PHI_BIN_MIN_EVENTS),
            }
        )
        support.setdefault("final_phi_event_counts", [int(value) for value in counts])
        support.setdefault("support_policy", "all_phi_bins")
        support.setdefault("support_threshold", int(PHI_BIN_MIN_EVENTS))
        if bool(support.get("passed")):
            break
        if not bool(auto_rebin_phi):
            break
        num_phi_bins -= 1
        if num_phi_bins < effective_minimum:
            raise ValueError(
                "No phi grid from {} down to {} passes {} support.".format(
                    original_requested, effective_minimum, support["support_policy"]
                )
            )
        bin_edges = _build_uniform_phi_bins(num_phi_bins)
        reduction_iterations += 1

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

    diagnostics = {
        "requested_num_phi_bins": int(original_requested),
        "initial_num_phi_bins": int(initial_num_phi_bins),
        "actual_num_phi_bins": int(len(bins) - 1),
        "phi_bin_reduction_applied": bool(len(bins) - 1 < initial_num_phi_bins),
        "phi_bin_reduction_reason": (
            "support_requirement" if len(bins) - 1 < initial_num_phi_bins else None
        ),
        "phi_bin_reduction_iterations": int(reduction_iterations),
        "minimum_phi_events": int(support["support_threshold"]),
        "minimum_phi_bins": int(effective_minimum),
        "auto_rebin_phi": bool(auto_rebin_phi),
        "final_phi_event_counts": list(support["final_phi_event_counts"]),
        "support": _json_ready(support),
        "status": "accepted" if bool(support.get("passed")) else "continued_with_requested_phi_binning",
        "continued_with_requested_phi_binning": bool(
            not bool(support.get("passed")) and not bool(auto_rebin_phi)
        ),
    }
    if return_diagnostics:
        return bins, counts, diagnostics
    return bins, counts


def _build_uniform_phi_bins(num_phi_bins):
    return np.linspace(float(PHI_BIN_MIN_DEG), float(PHI_BIN_MAX_DEG), int(num_phi_bins) + 1)


def evaluate_t_phi_raw_support(records, t_edges, phi_edges, *, minimum_events=PHI_BIN_MIN_EVENTS, policy="all_cells"):
    """Evaluate canonical raw-record support without using signed yields.

    ``records`` intentionally uses the same selected pre-particle record
    representation as the canonical t resolver.  Detector-setting coverage is
    therefore naturally aggregateable while the caller may retain source-wise
    diagnostics separately.
    """
    if str(policy).strip().lower() != "all_cells":
        raise ValueError("Unsupported t/phi support policy '{}'".format(policy))
    t_values = []
    phi_degrees = []
    for record in records or ():
        try:
            t_value = float(record["adj_t"])
            phi_value = float(record["phi_value"]) * (180.0 / math.pi)
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(t_value) and math.isfinite(phi_value):
            t_values.append(t_value)
            phi_degrees.append(phi_value)
    counts, _, _ = np.histogram2d(
        np.asarray(t_values, dtype=float),
        np.asarray(phi_degrees, dtype=float),
        bins=(np.asarray(t_edges, dtype=float), np.asarray(phi_edges, dtype=float)),
    )
    counts = counts.astype(int)
    threshold = int(minimum_events)
    supported = counts >= threshold
    marginal = (counts > 0) & ~supported
    unsupported = counts == 0
    return {
        "support_policy": "all_cells",
        "support_threshold": threshold,
        "raw_event_count_matrix": counts.tolist(),
        "supported_cells": int(np.count_nonzero(supported)),
        "marginal_cells": int(np.count_nonzero(marginal)),
        "unsupported_cells": int(np.count_nonzero(unsupported)),
        "total_cells": int(counts.size),
        "support_fraction": float(np.count_nonzero(supported) / counts.size) if counts.size else 0.0,
        "passed": bool(counts.size and np.all(supported)),
        "final_phi_event_counts": [int(value) for value in np.sum(counts, axis=0)],
    }


def resolve_phi_bins_from_t_phi_support(
    records,
    t_edges,
    requested_num_phi_bins,
    *,
    auto_rebin_phi=True,
    min_phi_bins=None,
    minimum_events=PHI_BIN_MIN_EVENTS,
    policy="all_cells",
    quiet=False,
):
    """Use the established equal-width reducer with canonical t/phi support."""
    phi_values = [
        float(record["phi_value"]) * (180.0 / math.pi)
        for record in records or ()
        if isinstance(record, dict) and "phi_value" in record
    ]
    if not phi_values:
        raise ValueError("No selected phi records available for t/phi support")
    return _find_phi_bins(
        np.asarray(phi_values, dtype=float),
        requested_num_phi_bins,
        quiet=quiet,
        return_diagnostics=True,
        auto_rebin_phi=auto_rebin_phi,
        min_phi_bins=min_phi_bins,
        support_evaluator=lambda edges: evaluate_t_phi_raw_support(
            records,
            t_edges,
            edges,
            minimum_events=minimum_events,
            policy=policy,
        ),
    )


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


def _canonical_interval_pair_payload(t_metadata, phi_metadata):
    """Return the stable identity shared by canonical t and phi sidecars."""
    t_metadata = dict(t_metadata or {})
    phi_metadata = dict(phi_metadata or {})
    return {
        "kinematic_key": t_metadata.get("kinematic_key"),
        "particle_type": t_metadata.get("particle_type"),
        "source_epsilon": t_metadata.get("source_epsilon"),
        "requested_num_t_bins": t_metadata.get("requested_num_t_bins"),
        "actual_num_t_bins": t_metadata.get("actual_num_t_bins"),
        "requested_num_phi_bins": phi_metadata.get("requested_num_phi_bins"),
        "actual_num_phi_bins": phi_metadata.get("actual_num_phi_bins"),
        "t_edges": t_metadata.get("t_edges"),
        "phi_edges": phi_metadata.get("phi_edges"),
        "t_binning_config_hash": t_metadata.get("binning_config_hash"),
        "phi_binning_config_hash": phi_metadata.get("binning_config_hash"),
        "t_algorithm_identifier": t_metadata.get("algorithm_identifier"),
        "t_algorithm_version": t_metadata.get("algorithm_version"),
        "phi_algorithm_identifier": phi_metadata.get("algorithm_identifier"),
        "phi_algorithm_version": phi_metadata.get("algorithm_version"),
    }


def _canonical_interval_pair_hash(t_metadata, phi_metadata):
    return _canonical_binning_config_hash(
        _canonical_interval_pair_payload(t_metadata, phi_metadata)
    )


def _validate_canonical_interval_pair(t_validation, phi_validation):
    """Ensure individually valid sidecars are from the same resolution."""
    reasons = []
    t_metadata = dict((t_validation or {}).get("metadata") or {})
    phi_metadata = dict((phi_validation or {}).get("metadata") or {})
    t_pair_id = t_metadata.get("canonical_interval_pair_id")
    phi_pair_id = phi_metadata.get("canonical_interval_pair_id")
    t_pair_hash = t_metadata.get("canonical_interval_pair_hash")
    phi_pair_hash = phi_metadata.get("canonical_interval_pair_hash")
    if not t_pair_id:
        reasons.append("t_pair_id_missing")
    if not phi_pair_id:
        reasons.append("phi_pair_id_missing")
    if not t_pair_hash:
        reasons.append("t_pair_hash_missing")
    if not phi_pair_hash:
        reasons.append("phi_pair_hash_missing")
    if t_pair_id and phi_pair_id and str(t_pair_id) != str(phi_pair_id):
        reasons.append("pair_id_mismatch")
    if t_pair_hash and phi_pair_hash and str(t_pair_hash) != str(phi_pair_hash):
        reasons.append("pair_hash_mismatch")
    if t_metadata and phi_metadata:
        try:
            expected_hash = _canonical_interval_pair_hash(t_metadata, phi_metadata)
        except (TypeError, ValueError):
            expected_hash = None
            reasons.append("pair_hash_uncomputable")
        if expected_hash is not None:
            if t_pair_hash and str(t_pair_hash) != expected_hash:
                reasons.append("t_pair_hash_recomputed_mismatch")
            if phi_pair_hash and str(phi_pair_hash) != expected_hash:
                reasons.append("phi_pair_hash_recomputed_mismatch")
    return {
        "valid": not reasons,
        "validation_rejection_reasons": reasons,
        "canonical_interval_pair_id": t_pair_id if t_pair_id == phi_pair_id else None,
        "canonical_interval_pair_hash": t_pair_hash if t_pair_hash == phi_pair_hash else None,
    }


def _build_canonical_phi_metadata(canonical_metadata, phi_edges, paths):
    """Build the phi sidecar from the one canonical pair payload."""
    metadata = dict(canonical_metadata or {})
    return {
        "schema_version": metadata.get("schema_version", 1),
        "particle_type": metadata.get("particle_type"),
        "Q2_token": metadata.get("Q2_token"),
        "W_token": metadata.get("W_token"),
        "kinematic_key": metadata.get("kinematic_key"),
        "source_epsilon": metadata.get("source_epsilon"),
        "consumer_epsilon": metadata.get("consumer_epsilon"),
        "shared_with_high_epsilon": bool(metadata.get("shared_with_high_epsilon", False)),
        "source_stage": metadata.get("source_stage"),
        "interval_file": paths["phi_path"],
        "metadata_file": paths["phi_metadata_path"],
        "phi_edges": _json_ready(np.asarray(phi_edges, dtype=float)),
        "requested_num_phi_bins": int(metadata.get("requested_num_phi_bins")),
        "actual_num_phi_bins": int(len(phi_edges) - 1),
        "phi_bin_reduction_applied": bool(metadata.get("phi_bin_reduction_applied", False)),
        "phi_bin_reduction_reason": metadata.get("phi_bin_reduction_reason"),
        "phi_bin_reduction_iterations": int(metadata.get("phi_bin_reduction_iterations", 0)),
        "final_phi_event_counts": metadata.get("final_phi_event_counts", []),
        "t_phi_support_policy": metadata.get("t_phi_support_policy"),
        "t_phi_support_min_events": metadata.get("t_phi_support_min_events"),
        "shared_epsilon_support": metadata.get("shared_epsilon_support", {}),
        "shared_preflight_stage": metadata.get("shared_preflight_stage"),
        "binning_config": metadata.get("phi_binning_config"),
        "binning_config_hash": metadata.get("phi_binning_config_hash"),
        "algorithm_identifier": metadata.get("phi_algorithm_identifier", "find_bins_phi_minimum_events"),
        "algorithm_version": metadata.get("phi_algorithm_version", 1),
        "canonical_interval_pair_id": metadata.get("canonical_interval_pair_id"),
        "canonical_interval_pair_hash": metadata.get("canonical_interval_pair_hash"),
        "input_provenance": metadata.get("input_provenance", {}),
        "created_before_proton_cleaning": True,
        "creation_timestamp": metadata.get("creation_timestamp"),
    }


def _canonical_phi_binning_semantic_config(inp_dict, t_binning_config):
    """Stable choices that determine the shared canonical phi proposal."""
    requested = int(inp_dict["NumPhiBins"])
    configured_minimum = int(inp_dict.get("min_phi_bins", MIN_PHI_BINS))
    return {
        "requested_num_phi_bins": requested,
        "phi_min_deg": float(PHI_BIN_MIN_DEG),
        "phi_max_deg": float(PHI_BIN_MAX_DEG),
        # This is a reduction floor, not permission to increase an explicit
        # lower-resolution request such as Nphi=7.
        "minimum_phi_bins": int(min(requested, max(1, configured_minimum))),
        "minimum_phi_events": int(
            inp_dict.get("t_phi_support_min_events", PHI_BIN_MIN_EVENTS)
        ),
        "auto_rebin_phi": bool(inp_dict.get("auto_rebin_phi", True)),
        "t_phi_support_policy": str(inp_dict.get("t_phi_support_policy", "all_cells")),
        "algorithm_identifier": "find_bins_t_phi_raw_support",
        "algorithm_version": 2,
    }


def _read_text_interval_edges(interval_path):
    with open(interval_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    if len(lines) < 2:
        raise ValueError("interval file has no edge line")
    tokens = [token.strip() for token in lines[1].split("\t") if token.strip()]
    if len(tokens) < 2:
        raise ValueError("interval file has fewer than two edges")
    return np.asarray([float(token) for token in tokens], dtype=float)


def _read_interval_metadata(metadata_path, reasons):
    if not os.path.exists(metadata_path):
        reasons.append("metadata_sidecar_missing")
        return None
    try:
        with open(metadata_path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except Exception as exc:
        reasons.append("metadata_sidecar_unreadable:{}".format(exc))
        return None


def _validate_authoritative_phi_interval(inp_dict, consumer_epsilon):
    """Validate the phi member of the low-epsilon canonical interval pair."""
    paths = _interval_paths(inp_dict)
    t_binning = _canonical_t_binning_config(inp_dict)
    reasons = []
    metadata = None
    text_edges = None
    if not os.path.exists(paths["phi_path"]):
        reasons.append("interval_file_missing")
    else:
        try:
            text_edges = _read_text_interval_edges(paths["phi_path"])
        except Exception as exc:
            reasons.append("interval_file_unreadable:{}".format(exc))
    if bool(t_binning.get("require_phi_metadata_sidecar", True)):
        metadata = _read_interval_metadata(paths["phi_metadata_path"], reasons)
    elif os.path.exists(paths["phi_metadata_path"]):
        metadata = _read_interval_metadata(paths["phi_metadata_path"], reasons)
    else:
        reasons.append("metadata_sidecar_missing")

    if metadata is not None:
        expected_config = _canonical_phi_binning_semantic_config(inp_dict, t_binning)
        expected_hash = _canonical_binning_config_hash(expected_config)
        checks = (
            (metadata.get("schema_version") == int(t_binning.get("metadata_schema_version", 1)), "metadata_schema_version_mismatch"),
            (metadata.get("particle_type") == inp_dict.get("ParticleType"), "particle_type_mismatch"),
            (metadata.get("Q2_token") == inp_dict.get("Q2"), "q2_token_mismatch"),
            (metadata.get("W_token") == inp_dict.get("W"), "w_token_mismatch"),
            (metadata.get("kinematic_key") == "Q{}W{}".format(inp_dict.get("Q2"), inp_dict.get("W")), "kinematic_key_mismatch"),
            (metadata.get("source_stage") == "pre_particle_subtraction", "source_stage_mismatch"),
            (metadata.get("source_epsilon") == "low", "source_epsilon_mismatch"),
            (bool(metadata.get("shared_with_high_epsilon", False)), "low_to_high_sharing_not_enabled"),
            (metadata.get("requested_num_phi_bins") == int(inp_dict["NumPhiBins"]), "requested_bin_count_mismatch"),
            (metadata.get("binning_config_hash") == expected_hash, "binning_config_hash_mismatch"),
            (metadata.get("algorithm_identifier") in set(t_binning.get("allowed_phi_algorithm_identifiers") or ("find_bins_phi_minimum_events", "find_bins_t_phi_raw_support")), "algorithm_identifier_not_allowed"),
            (metadata.get("algorithm_version") in set(t_binning.get("allowed_phi_algorithm_versions") or (1, 2)), "algorithm_version_not_allowed"),
        )
        for valid, reason in checks:
            if not valid:
                reasons.append(reason)
        try:
            metadata_edges = np.asarray(metadata.get("phi_edges") or [], dtype=float)
        except (TypeError, ValueError):
            metadata_edges = np.asarray([], dtype=float)
            reasons.append("metadata_edges_invalid")
        tolerance = float(t_binning.get("edge_tolerance", 1.0e-9))
        if metadata_edges.size and metadata.get("actual_num_phi_bins") != int(metadata_edges.size - 1):
            reasons.append("metadata_actual_bin_count_edge_count_mismatch")
        if metadata_edges.size and (not np.all(np.isfinite(metadata_edges)) or np.any(np.diff(metadata_edges) <= 0.0)):
            reasons.append("metadata_edges_not_strictly_increasing_finite")
        if metadata_edges.size and (
            abs(float(metadata_edges[0]) - float(PHI_BIN_MIN_DEG)) > tolerance
            or abs(float(metadata_edges[-1]) - float(PHI_BIN_MAX_DEG)) > tolerance
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
        "interval_file": paths["phi_path"],
        "metadata_file": paths["phi_metadata_path"],
        "metadata": metadata,
        "phi_edges": text_edges,
    }


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


def validate_canonical_phi_edges(inp_dict, phi_edges, context=""):
    """Record and enforce the companion immutable phi edges."""
    canonical = dict(inp_dict.get("canonical_t_binning") or {})
    expected_source = canonical.get("phi_edges")
    if expected_source is None:
        expected_source = inp_dict.get("phi_bins")
    expected = np.asarray(expected_source if expected_source is not None else [], dtype=float)
    candidate = np.asarray(phi_edges if phi_edges is not None else [], dtype=float)
    tolerance = float(_canonical_t_binning_config(inp_dict).get("edge_tolerance", 1.0e-9))
    match = bool(expected.size == candidate.size and np.allclose(expected, candidate, rtol=0.0, atol=tolerance))
    max_difference = float(np.max(np.abs(expected - candidate))) if expected.size == candidate.size and expected.size else float("inf")
    diagnostics = {
        "context": str(context),
        "canonical_phi_edge_match": match,
        "maximum_phi_edge_difference": max_difference,
    }
    if not match and bool(get_proton_contamination_cleaning_config(inp_dict=inp_dict).get("strict_mode", False)):
        raise RuntimeError("canonical phi edges differ in {}".format(context or "unknown context"))
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
    phi_interval_validation = {
        "valid": False,
        "validation_status": "no_valid_source",
        "validation_rejection_reasons": ["authoritative_interval_disabled"],
    }
    pair_validation = {
        "valid": False,
        "validation_rejection_reasons": ["authoritative_interval_disabled"],
    }
    if allow_interval_file and bool(t_binning.get("allow_authoritative_interval_file", True)):
        interval_validation = _validate_authoritative_t_interval(inp_dict, consumer_epsilon)
        phi_interval_validation = _validate_authoritative_phi_interval(inp_dict, consumer_epsilon)
        pair_validation = _validate_canonical_interval_pair(
            interval_validation, phi_interval_validation
        )

    interval_pair_valid = bool(
        interval_validation.get("valid")
        and phi_interval_validation.get("valid")
        and pair_validation.get("valid")
    )
    pair_rejection_reasons = list(interval_validation.get("validation_rejection_reasons") or []) + [
        "phi_{}".format(reason)
        for reason in (phi_interval_validation.get("validation_rejection_reasons") or [])
    ] + ["pair_{}".format(reason) for reason in (pair_validation.get("validation_rejection_reasons") or [])]
    if interval_pair_valid:
        t_edges = np.asarray(interval_validation["t_edges"], dtype=float)
        phi_edges = np.asarray(phi_interval_validation["phi_edges"], dtype=float)
        source = "validated_authoritative_interval_file"
        validation_status = source
        metadata = dict(interval_validation.get("metadata") or {})
        phi_metadata = dict(phi_interval_validation.get("metadata") or {})
        phi_resolution = {
            "requested_num_phi_bins": int(phi_metadata.get("requested_num_phi_bins")),
            "actual_num_phi_bins": int(phi_metadata.get("actual_num_phi_bins")),
            "phi_bin_reduction_applied": bool(phi_metadata.get("phi_bin_reduction_applied", False)),
            "phi_bin_reduction_reason": phi_metadata.get("phi_bin_reduction_reason"),
            "phi_bin_reduction_iterations": int(phi_metadata.get("phi_bin_reduction_iterations", 0)),
            "final_phi_event_counts": list(phi_metadata.get("final_phi_event_counts") or []),
            "status": "authoritative_reuse",
        }
    elif consumer_epsilon == "high":
        # High epsilon is a consumer of low-epsilon canonical bins, never a
        # second independent producer for the same kinematic key.
        canonical = {
            "source": "no_valid_source",
            "validation_status": "no_valid_source",
            "validation_rejection_reasons": pair_rejection_reasons,
            "consumer_epsilon": "high",
            "shared_from_epsilon": "low",
            "created_before_proton_cleaning": True,
        }
        inp_dict["canonical_t_binning"] = canonical
        raise RuntimeError("No compatible low-epsilon canonical t/phi interval pair is available for high epsilon.")
    else:
        if bool(inp_dict.get("require_shared_canonical_preflight", False)) and not bool(
            inp_dict.get("allow_unpaired_canonical_binning", False)
        ):
            inp_dict["canonical_t_binning"] = {
                "source": "unpaired_canonical_preflight_missing",
                "resolution_status": "non_pairable_rejected",
                "validation_rejection_reasons": pair_rejection_reasons,
                "consumer_epsilon": "low",
                "shared_with_high_epsilon": False,
                "created_before_proton_cleaning": True,
            }
            raise RuntimeError(
                "No shared low/high canonical t/phi preflight pair is available. "
                "Set LT_ANALYSIS_ALLOW_UNPAIRED_CANONICAL_BINNING only for explicit development use."
            )
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
        phi_metadata = {}
        phi_degrees = phi_values * (180.0 / math.pi)
        phi_edges, _, phi_resolution = resolve_phi_bins_from_t_phi_support(
            records,
            t_edges,
            requested_phi,
            quiet=quiet,
            auto_rebin_phi=bool(inp_dict.get("auto_rebin_phi", True)),
            min_phi_bins=inp_dict.get("min_phi_bins", MIN_PHI_BINS),
            minimum_events=int(inp_dict.get("t_phi_support_min_events", PHI_BIN_MIN_EVENTS)),
            policy=str(inp_dict.get("t_phi_support_policy", "all_cells")),
        )

    phi_degrees = phi_values * (180.0 / math.pi)
    phi_counts, _ = np.histogram(phi_degrees, phi_edges)
    support = _support_summaries(t_values, physical_weights, t_edges)
    semantic_config = _canonical_binning_semantic_config(inp_dict, t_binning)
    phi_semantic_config = _canonical_phi_binning_semantic_config(inp_dict, t_binning)
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
        "validation_rejection_reasons": pair_rejection_reasons,
        "interval_file": paths["t_path"],
        "metadata_file": paths["t_metadata_path"],
        "requested_num_t_bins": requested_t,
        "actual_num_t_bins": len(t_edges) - 1,
        "tmin": float(inp_dict["tmin"]),
        "tmax": float(inp_dict["tmax"]),
        "t_edges": _json_ready(t_edges),
        "phi_interval_file": paths["phi_path"],
        "phi_metadata_file": paths["phi_metadata_path"],
        "requested_num_phi_bins": requested_phi,
        "actual_num_phi_bins": len(phi_edges) - 1,
        "phi_edges": _json_ready(phi_edges),
        "phi_bin_reduction_applied": bool(phi_resolution["phi_bin_reduction_applied"]),
        "phi_bin_reduction_reason": phi_resolution["phi_bin_reduction_reason"],
        "phi_bin_reduction_iterations": int(phi_resolution["phi_bin_reduction_iterations"]),
        "final_phi_event_counts": list(phi_resolution["final_phi_event_counts"]),
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
        "phi_binning_config": _json_ready(phi_semantic_config),
        "phi_binning_config_hash": _canonical_binning_config_hash(phi_semantic_config),
        "phi_algorithm_identifier": "find_bins_t_phi_raw_support",
        "phi_algorithm_version": 2,
        "algorithm_identifier": "find_bins_adjust_t_bins",
        "algorithm_version": 1,
        "input_provenance": {"source_stats": source_stats, "record_count": int(len(records))},
        "created_before_proton_cleaning": True,
        "creation_timestamp": datetime.now(timezone.utc).isoformat(),
        "edge_match": True,
        "canonical_interval_pair_validation": _json_ready(pair_validation),
    }
    # Preserve source-side provenance when a valid authoritative file wins,
    # while exposing current-run support in the active resolution payload.
    if metadata:
        canonical["authoritative_metadata"] = _json_ready(metadata)
    if phi_metadata:
        canonical["authoritative_phi_metadata"] = _json_ready(phi_metadata)
    if source == "validated_authoritative_interval_file":
        canonical["canonical_interval_pair_id"] = metadata.get("canonical_interval_pair_id")
        canonical["canonical_interval_pair_hash"] = metadata.get("canonical_interval_pair_hash")
    else:
        canonical["canonical_interval_pair_id"] = uuid.uuid4().hex
        pair_phi_metadata = _build_canonical_phi_metadata(canonical, phi_edges, paths)
        canonical["canonical_interval_pair_hash"] = _canonical_interval_pair_hash(
            canonical, pair_phi_metadata
        )
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


def resolve_shared_canonical_phi_preflight(
    prepass_by_epsilon,
    inp_dict,
    *,
    write_interval_files=True,
    quiet=False,
):
    """Resolve one low/high-epsilon canonical interval pair from raw support.

    This helper deliberately consumes only the lightweight records captured
    before particle subtraction.  It is suitable for the launcher-level common
    orchestration point: no fitted amplitudes, normalized yields, or signed
    weights contribute to the grid decision.  The t proposal remains
    low-epsilon-authoritative for compatibility; only the phi grid is required
    to pass raw support for *both* epsilon members.
    """
    if not isinstance(prepass_by_epsilon, dict):
        raise TypeError("prepass_by_epsilon must map low/high epsilon to prepass payloads")

    def _epsilon_payloads(epsilon):
        payloads = prepass_by_epsilon.get(epsilon)
        if isinstance(payloads, dict):
            payloads = list(payloads.values())
        return list(payloads or ())

    low_records = _prepass_records(_epsilon_payloads("low"))
    high_records = _prepass_records(_epsilon_payloads("high"))
    if not low_records or not high_records:
        missing = [
            epsilon
            for epsilon, records in (("low", low_records), ("high", high_records))
            if not records
        ]
        raise ValueError(
            "shared_canonical_preflight_missing_selected_records:{}".format(
                ",".join(missing)
            )
        )

    requested_t = int(inp_dict["NumtBins"])
    requested_phi = int(inp_dict["NumPhiBins"])
    t_binning = _canonical_t_binning_config(inp_dict)
    low_t_values = np.asarray([record["adj_t"] for record in low_records], dtype=float)
    low_weights = np.asarray(
        [record["physical_coefficient"] for record in low_records], dtype=float
    )
    proposal_values = low_t_values.copy()
    proposal_values[np.argmin(proposal_values)] = float(inp_dict["tmin"])
    proposal_values[np.argmax(proposal_values)] = float(inp_dict["tmax"])
    t_edges, _ = _find_t_bins_from_prepass(
        proposal_values,
        low_weights,
        requested_t,
        str(t_binning.get("canonical_bin_support_metric", "raw_event_count")),
        t_binning.get("canonical_bin_support_threshold", T_BIN_MIN_EVENTS),
        quiet=quiet,
    )
    policy = str(inp_dict.get("t_phi_support_policy", "all_cells"))
    threshold = int(inp_dict.get("t_phi_support_min_events", PHI_BIN_MIN_EVENTS))
    all_records = low_records + high_records

    def _shared_support(phi_edges):
        low_support = evaluate_t_phi_raw_support(
            low_records, t_edges, phi_edges, minimum_events=threshold, policy=policy
        )
        high_support = evaluate_t_phi_raw_support(
            high_records, t_edges, phi_edges, minimum_events=threshold, policy=policy
        )
        return {
            "passed": bool(low_support["passed"] and high_support["passed"]),
            "support_policy": policy,
            "support_threshold": threshold,
            "final_phi_event_counts": [
                int(value)
                for value in np.histogram(
                    [record["phi_value"] * (180.0 / math.pi) for record in all_records],
                    bins=phi_edges,
                )[0]
            ],
            "epsilon_support": {"low": low_support, "high": high_support},
        }

    phi_degrees = np.asarray(
        [record["phi_value"] * (180.0 / math.pi) for record in all_records], dtype=float
    )
    phi_edges, phi_counts, phi_resolution = _find_phi_bins(
        phi_degrees,
        requested_phi,
        quiet=quiet,
        return_diagnostics=True,
        auto_rebin_phi=bool(inp_dict.get("auto_rebin_phi", True)),
        min_phi_bins=inp_dict.get("min_phi_bins", MIN_PHI_BINS),
        support_evaluator=_shared_support,
    )
    support = _support_summaries(low_t_values, low_weights, t_edges)
    paths = _interval_paths(inp_dict)
    semantic_config = _canonical_binning_semantic_config(inp_dict, t_binning)
    phi_semantic_config = _canonical_phi_binning_semantic_config(inp_dict, t_binning)
    shared_support = dict(phi_resolution.get("support") or {})
    canonical = {
        "schema_version": int(t_binning.get("metadata_schema_version", 1)),
        "particle_type": inp_dict["ParticleType"],
        "Q2_token": inp_dict["Q2"],
        "W_token": inp_dict["W"],
        "kinematic_key": "Q{}W{}".format(inp_dict["Q2"], inp_dict["W"]),
        "source_epsilon": "low",
        "consumer_epsilon": "low",
        "shared_with_high_epsilon": True,
        # Keep the established authoritative-interval contract: this is still
        # pre-particle-subtraction data, with the shared orchestration recorded
        # separately below.
        "source_stage": "pre_particle_subtraction",
        "shared_preflight_stage": "low_high_common_orchestration",
        "source": "shared_low_high_raw_support_preflight",
        "resolution_status": "computed_shared_preflight",
        "validation_status": "computed_shared_preflight",
        "interval_file": paths["t_path"],
        "metadata_file": paths["t_metadata_path"],
        "requested_num_t_bins": requested_t,
        "actual_num_t_bins": int(len(t_edges) - 1),
        "tmin": float(inp_dict["tmin"]),
        "tmax": float(inp_dict["tmax"]),
        "t_edges": _json_ready(t_edges),
        "phi_interval_file": paths["phi_path"],
        "phi_metadata_file": paths["phi_metadata_path"],
        "requested_num_phi_bins": requested_phi,
        "actual_num_phi_bins": int(len(phi_edges) - 1),
        "phi_edges": _json_ready(phi_edges),
        "phi_bin_reduction_applied": bool(phi_resolution["phi_bin_reduction_applied"]),
        "phi_bin_reduction_reason": phi_resolution["phi_bin_reduction_reason"],
        "phi_bin_reduction_iterations": int(phi_resolution["phi_bin_reduction_iterations"]),
        "final_phi_event_counts": list(phi_resolution["final_phi_event_counts"]),
        "shared_epsilon_support": _json_ready(shared_support.get("epsilon_support") or {}),
        "t_phi_support_policy": policy,
        "t_phi_support_min_events": threshold,
        "raw_event_count_by_t_bin": _json_ready(support["raw_event_count_by_t_bin"]),
        "binning_config": _json_ready(semantic_config),
        "binning_config_hash": _canonical_binning_config_hash(semantic_config),
        "phi_binning_config": _json_ready(phi_semantic_config),
        "phi_binning_config_hash": _canonical_binning_config_hash(phi_semantic_config),
        "phi_algorithm_identifier": "find_bins_t_phi_raw_support",
        "phi_algorithm_version": 2,
        "algorithm_identifier": "find_bins_adjust_t_bins",
        "algorithm_version": 1,
        "input_provenance": {
            "low_record_count": int(len(low_records)),
            "high_record_count": int(len(high_records)),
            "selection": "raw_pre_particle_selected_records",
        },
        "created_before_proton_cleaning": True,
        "creation_timestamp": datetime.now(timezone.utc).isoformat(),
        "edge_match": True,
        "canonical_interval_pair_id": uuid.uuid4().hex,
    }
    phi_metadata = _build_canonical_phi_metadata(canonical, phi_edges, paths)
    canonical["canonical_interval_pair_hash"] = _canonical_interval_pair_hash(
        canonical, phi_metadata
    )
    inp_dict["t_bins"] = np.asarray(t_edges, dtype=float)
    inp_dict["phi_bins"] = np.asarray(phi_edges, dtype=float)
    inp_dict["canonical_t_binning"] = canonical
    artifacts = _interval_paths(inp_dict)
    if write_interval_files:
        artifacts = write_bin_interval_files(inp_dict, t_edges, phi_edges, canonical_metadata=canonical)
    canonical["artifacts"] = artifacts
    return {
        "t_bins": np.asarray(t_edges, dtype=float),
        "phi_bins": np.asarray(phi_edges, dtype=float),
        "phi_counts": np.asarray(phi_counts, dtype=int),
        "metadata": canonical,
        "support": shared_support,
    }


def write_bin_interval_files(inpDict, t_bins, phi_bins, canonical_metadata=None):
    paths = _interval_paths(inpDict)
    phi_path = paths["phi_path"]
    t_path = paths["t_path"]

    # A legacy bin finder must never replace canonical text while leaving an
    # old sidecar behind.  Canonical output is owned by the pre-subtraction
    # resolver, which supplies one matching metadata payload for both files.
    if canonical_metadata is None and (
        bool(inpDict.get("canonical_t_binning"))
        or os.path.exists(paths["t_metadata_path"])
        or os.path.exists(paths["phi_metadata_path"])
    ):
        raise RuntimeError(
            "canonical_interval_overwrite_refused: use resolve_canonical_analysis_bins_pre_subtraction"
        )

    header = "{}\t{}\t{}\t{}\n".format(
        inpDict["Q2"].replace("p", "."),
        inpDict["W"].replace("p", "."),
        len(t_bins) - 1,
        len(phi_bins) - 1,
    )
    phi_text = header + "".join("\t{:.17g}".format(float(phi)) for phi in phi_bins)
    t_text = header + "".join("\t{:.17g}".format(float(value)) for value in t_bins)
    artifacts = dict(paths)
    if canonical_metadata is None:
        with open(phi_path, "w", encoding="utf-8") as handle:
            handle.write(phi_text)
        with open(t_path, "w", encoding="utf-8") as handle:
            handle.write(t_text)
        return artifacts

    # Build all members before publishing any of them.  Individual replaces
    # cannot be a filesystem transaction, so the shared pair identity also
    # makes an interrupted mixed generation non-authoritative.
    metadata = dict(canonical_metadata)
    metadata.update({
        "interval_file": t_path,
        "metadata_file": paths["t_metadata_path"],
        "t_edges": _json_ready(np.asarray(t_bins, dtype=float)),
        "actual_num_t_bins": int(len(t_bins) - 1),
        "created_before_proton_cleaning": True,
        "canonical_interval_pair_id": metadata.get("canonical_interval_pair_id") or uuid.uuid4().hex,
    })
    phi_metadata = _build_canonical_phi_metadata(metadata, phi_bins, paths)
    metadata["canonical_interval_pair_hash"] = _canonical_interval_pair_hash(metadata, phi_metadata)
    phi_metadata = _build_canonical_phi_metadata(metadata, phi_bins, paths)
    if metadata["canonical_interval_pair_hash"] != _canonical_interval_pair_hash(metadata, phi_metadata):
        raise RuntimeError("canonical interval pair hash construction failed")

    staged = {}
    payloads = {
        "t_path": (t_path, t_text),
        "phi_path": (phi_path, phi_text),
        "t_metadata_path": (
            paths["t_metadata_path"],
            json.dumps(_json_ready(metadata), sort_keys=True, separators=(",", ":"), allow_nan=False),
        ),
        "phi_metadata_path": (
            paths["phi_metadata_path"],
            json.dumps(_json_ready(phi_metadata), sort_keys=True, separators=(",", ":"), allow_nan=False),
        ),
    }
    try:
        for key, (target, content) in payloads.items():
            descriptor, temporary_path = tempfile.mkstemp(
                prefix=".{}-".format(os.path.basename(target)), suffix=".tmp", dir=os.path.dirname(target)
            )
            with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
            staged[key] = temporary_path
        staged_t_metadata = _read_interval_metadata(staged["t_metadata_path"], [])
        staged_phi_metadata = _read_interval_metadata(staged["phi_metadata_path"], [])
        staged_pair = _validate_canonical_interval_pair(
            {"metadata": staged_t_metadata}, {"metadata": staged_phi_metadata}
        )
        if not staged_pair["valid"]:
            raise RuntimeError("canonical interval temporary pair is invalid: {}".format(
                ", ".join(staged_pair["validation_rejection_reasons"])
            ))
        if not np.allclose(_read_text_interval_edges(staged["t_path"]), np.asarray(t_bins, dtype=float), rtol=0.0, atol=0.0):
            raise RuntimeError("canonical temporary t interval does not round-trip")
        if not np.allclose(_read_text_interval_edges(staged["phi_path"]), np.asarray(phi_bins, dtype=float), rtol=0.0, atol=0.0):
            raise RuntimeError("canonical temporary phi interval does not round-trip")
        for key in ("t_path", "phi_path", "t_metadata_path", "phi_metadata_path"):
            os.replace(staged[key], payloads[key][0])
            staged.pop(key, None)
    finally:
        for temporary_path in staged.values():
            try:
                os.unlink(temporary_path)
            except OSError:
                pass
    return artifacts


def find_bins(histlist, inpDict):
    if bool(inpDict.get("canonical_t_binning")):
        raise RuntimeError(
            "canonical_interval_overwrite_refused: canonical bins are immutable after pre-subtraction resolution"
        )
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
        phi_edge_diagnostics = validate_canonical_phi_edges(
            inpDict, phi_bins, context="find_bins.check_bins"
        )
        inpDict["canonical_t_binning"]["check_bins_edge_validation"] = edge_diagnostics
        inpDict["canonical_t_binning"]["check_bins_phi_edge_validation"] = phi_edge_diagnostics

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
