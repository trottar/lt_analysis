#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2026-05-05 00:00:00 codex"
# ================================================================
#

import math
import os
import sys

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
    PHI_BIN_MIN_EVENTS,
    T_BIN_ADJUST_MAX_ITERATIONS,
    T_BIN_ADJUST_TOLERANCE,
    T_BIN_EDGE_BIAS,
    T_BIN_MIN_EVENTS,
)


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
        hist["H_t_DATA"].Scale(1.0 / normfac_data)
        hist["H_ph_q_DATA"].Scale(1.0 / normfac_data)

        t_values = flatten_hist(hist["H_t_DATA"])
        phi_deg = [(phi * (180.0 / math.pi)) for phi in flatten_hist(hist["H_ph_q_DATA"])]

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
    bin_edges = np.linspace(-180.0, 180.0, num_phi_bins + 1)

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
        bin_edges = np.linspace(-180.0, 180.0, num_phi_bins + 1)

    if not quiet:
        for i, val in enumerate(counts):
            print(
                "Bin {} from {:.1f} to {:.1f} has {} events".format(
                    i + 1, bins[i], bins[i + 1], val
                )
            )
        print("phi_bins = ", bins)

    return bins, counts


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


def write_bin_interval_files(inpDict, t_bins, phi_bins):
    particle_type = inpDict["ParticleType"]
    q2_token = inpDict["Q2"].replace("p", "")
    w_token = inpDict["W"].replace("p", "")

    phi_path = "{}/src/{}/phi_bin_interval_Q{}W{}".format(
        LTANAPATH, particle_type, q2_token, w_token
    )
    t_path = "{}/src/{}/t_bin_interval_Q{}W{}".format(
        LTANAPATH, particle_type, q2_token, w_token
    )

    phi_lines = ["\t{}".format(float(phi)) for phi in phi_bins]
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

    t_lines = ["\t{:.2f}".format(float(t_val)) for t_val in t_bins]
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

    return {"phi_path": phi_path, "t_path": t_path}


def find_bins(histlist, inpDict):
    proposal = propose_bins(histlist, inpDict, quiet=False)
    apply_bin_proposal(inpDict, proposal)
    write_bin_interval_files(inpDict, proposal["t_bins"], proposal["phi_bins"])
    return proposal


def check_bins(histlist, inpDict):
    t_values, phi_values = _collect_bin_samples(histlist, inpDict)
    t_bins = histlist[0]["t_bins"]
    phi_bins = histlist[0]["phi_bins"]

    print("\nFinding phi bins...")
    phi_counts, phi_edges = np.histogram(phi_values, phi_bins)
    if len(phi_bins) - 1 != inpDict["NumPhiBins"]:
        print(
            "Number of phi-bins changed from {} to: {}".format(
                inpDict["NumPhiBins"], len(phi_counts)
            )
        )
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
        inpDict["NumtBins"] = len(t_counts)
    for i, val in enumerate(t_counts):
        print(
            "Bin {} from {:.3f} to {:.3f} has {} events".format(
                i + 1, t_edges[i], t_edges[i + 1], val
            )
        )
    print("t_bins = ", t_edges)
