#! /usr/bin/python

"""
Shared analysis configuration for stage-2 background tuning and bin finding.

Edit the constants in this file when you want one place to control the
fallback BG scale, binning thresholds, and optimizer search behavior.
"""

from __future__ import annotations


BG_STAT_SCALE2 = 0.5

# Sigma cut used when ranking data/SIMC ratio agreement.
RATIO_SIGMA_THRESHOLD = 3.0

# Shared bin-finding controls migrated from find_bins.py.
PHI_BIN_MIN_EVENTS = 10
T_BIN_MIN_EVENTS = 500
MIN_PHI_BINS = 5
MIN_T_BINS = 2
PHI_BIN_MIN_DEG = -180.0
PHI_BIN_MAX_DEG = 180.0
T_BIN_ADJUST_TOLERANCE = 1e-3
T_BIN_ADJUST_MAX_ITERATIONS = 50000
T_BIN_EDGE_BIAS = 2.0

# Joint optimizer search space.
NUM_T_BIN_SEARCH_OFFSETS = (-2, -1, 0, 1, 2)
NUM_PHI_BIN_SEARCH_OFFSETS = (-2, -1, 0, 1, 2)

BG_STAT_SCALE2_COARSE_MIN = 0.0
BG_STAT_SCALE2_COARSE_MAX = 1.0
BG_STAT_SCALE2_COARSE_STEP = 0.1
BG_STAT_SCALE2_REFINE_WINDOW = 0.1
BG_STAT_SCALE2_REFINE_STEP = 0.025
BG_STAT_SCALE2_FINALIST_COUNT = 2

# Variables used for lightweight SIMC-vs-data kinematic scoring.
KINEMATIC_SCORE_VARS = ("Q2", "W", "theta_cm", "mm")


def _round_scale(value):
    return round(float(value), 6)


def get_bg_scale_setting_key(epsset, phi_setting):
    return "{}:{}".format(str(epsset).lower(), str(phi_setting))


def resolve_bg_stat_scale2(inpDict, phi_setting=None):
    if phi_setting is None:
        phi_setting = inpDict.get("phi_setting")
    key = get_bg_scale_setting_key(inpDict.get("EPSSET", ""), phi_setting)
    scale_map = inpDict.get("bg_stat_scale2_by_setting", {})
    if key in scale_map:
        return float(scale_map[key])
    if "bg_stat_scale2" in inpDict:
        return float(inpDict["bg_stat_scale2"])
    return float(BG_STAT_SCALE2)


def build_bin_count_candidates(requested_t_bins, requested_phi_bins):
    candidates = []
    seen = set()
    for t_offset in NUM_T_BIN_SEARCH_OFFSETS:
        for phi_offset in NUM_PHI_BIN_SEARCH_OFFSETS:
            t_bins = max(MIN_T_BINS, int(requested_t_bins) + int(t_offset))
            phi_bins = max(MIN_PHI_BINS, int(requested_phi_bins) + int(phi_offset))
            key = (t_bins, phi_bins)
            if key not in seen:
                seen.add(key)
                candidates.append(key)
    return candidates


def get_bg_scale_coarse_candidates():
    values = []
    candidate = BG_STAT_SCALE2_COARSE_MIN
    while candidate <= BG_STAT_SCALE2_COARSE_MAX + 1e-12:
        values.append(_round_scale(candidate))
        candidate += BG_STAT_SCALE2_COARSE_STEP
    if _round_scale(BG_STAT_SCALE2) not in values:
        values.append(_round_scale(BG_STAT_SCALE2))
    return sorted(set(values))


def get_bg_scale_refined_candidates(best_value):
    lo = max(BG_STAT_SCALE2_COARSE_MIN, float(best_value) - BG_STAT_SCALE2_REFINE_WINDOW)
    hi = min(BG_STAT_SCALE2_COARSE_MAX, float(best_value) + BG_STAT_SCALE2_REFINE_WINDOW)
    values = []
    candidate = lo
    while candidate <= hi + 1e-12:
        values.append(_round_scale(candidate))
        candidate += BG_STAT_SCALE2_REFINE_STEP
    values.append(_round_scale(best_value))
    return sorted(set(values))
