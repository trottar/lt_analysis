#! /usr/bin/python

import json
import math
import os
import subprocess
import sys
import csv

import numpy as np

from frozen_manifest import load_zeroth_iteration_input_bundle
from background_config import (
    BG_STAT_SCALE1,
    BG_STAT_SCALE2,
    BG_STAT_SCALE2_FINALIST_COUNT,
    KINEMATIC_SCORE_VARS,
    RATIO_SIGMA_THRESHOLD,
    build_bin_count_candidates,
    get_active_bg_profile,
    get_active_metric_weights,
    get_active_selection_mode,
    get_bg_scale1_coarse_candidates,
    get_bg_scale1_refined_candidates,
    get_bg_scale_coarse_candidates,
    get_bg_scale_refined_candidates,
    get_bg_scale_setting_key,
    get_common_epsilon_scale_behavior,
    get_forced_bg_scale1,
    get_forced_bg_scale2,
    get_resolved_bg_profile_settings,
    resolve_bg_stat_scale1,
    resolve_bg_stat_scale2,
    use_common_epsilon_scales,
)

sys.path.append("binning")
from find_bins import apply_bin_proposal, propose_bins, write_bin_interval_files


def _log(message):
    print("[BG OPT] {}".format(message))


def _format_metric_value(value, precision=4):
    try:
        val = float(value)
    except Exception:
        return "nan"
    if not math.isfinite(val):
        return "inf"
    return "{:.{prec}f}".format(val, prec=precision)


def _safe_float(value, default=float("inf")):
    try:
        val = float(value)
        if math.isfinite(val):
            return val
    except Exception:
        pass
    return default


def _get_selection_mode():
    return get_active_selection_mode()


def _normalize_metric_values(values, higher_is_better=False):
    numeric = np.asarray([_safe_float(value, default=float("nan")) for value in values], dtype=float)
    finite_mask = np.isfinite(numeric)
    normalized = np.zeros(len(numeric), dtype=float)
    if not finite_mask.any():
        return normalized

    finite_vals = numeric[finite_mask]
    lo = float(np.min(finite_vals))
    hi = float(np.max(finite_vals))
    if math.isclose(lo, hi, rel_tol=0.0, abs_tol=1.0e-12):
        return normalized

    if higher_is_better:
        normalized[finite_mask] = (hi - numeric[finite_mask]) / (hi - lo)
    else:
        normalized[finite_mask] = (numeric[finite_mask] - lo) / (hi - lo)
    return normalized


def _annotate_selection_scores(results):
    mode = _get_selection_mode()
    for result in results:
        result["selection_mode"] = mode
        result["selection_score"] = None

    if mode != "weighted":
        return

    valid_results = [result for result in results if result.get("valid")]
    if not valid_results:
        return

    metric_series = {
        "ratio_fail_count": _normalize_metric_values(
            [result.get("metrics", {}).get("ratio_fail_count") for result in valid_results],
            higher_is_better=False,
        ),
        "ratio_mean_dev": _normalize_metric_values(
            [result.get("metrics", {}).get("ratio_mean_dev") for result in valid_results],
            higher_is_better=False,
        ),
        "ratio_rms": _normalize_metric_values(
            [result.get("metrics", {}).get("ratio_rms") for result in valid_results],
            higher_is_better=False,
        ),
        "kinematic_score": _normalize_metric_values(
            [result.get("metrics", {}).get("kinematic_score") for result in valid_results],
            higher_is_better=False,
        ),
        "valid_ratio_bins": _normalize_metric_values(
            [result.get("metrics", {}).get("valid_ratio_bins") for result in valid_results],
            higher_is_better=True,
        ),
    }

    for idx, result in enumerate(valid_results):
        score = 0.0
        for metric_name, weight in get_active_metric_weights().items():
            if metric_name not in metric_series:
                continue
            score += float(weight) * float(metric_series[metric_name][idx])
        result["selection_score"] = float(score)


def _candidate_sort_key(result):
    return (
        int(result["metrics"]["ratio_fail_count"]),
        float(result["metrics"]["ratio_mean_dev"]),
        float(result["metrics"]["ratio_rms"]),
        float(result["metrics"]["kinematic_score"]),
    )


def _result_bin_count_key(result):
    proposal = result.get("proposal", {}) if isinstance(result, dict) else {}
    actual_t = proposal.get("actual_num_t_bins", result.get("actual_num_t_bins", result.get("requested_num_t_bins", 0)))
    actual_phi = proposal.get("actual_num_phi_bins", result.get("actual_num_phi_bins", result.get("requested_num_phi_bins", 0)))
    try:
        actual_t = int(actual_t)
    except Exception:
        actual_t = 0
    try:
        actual_phi = int(actual_phi)
    except Exception:
        actual_phi = 0
    return actual_t * actual_phi, actual_t, actual_phi


def _candidate_order_key(result):
    metrics = result.get("metrics", {})
    total_bins, actual_t, actual_phi = _result_bin_count_key(result)
    if _get_selection_mode() == "weighted":
        return (
            _safe_float(result.get("selection_score"), default=float("inf")),
            -int(total_bins),
            -int(actual_t),
            -int(actual_phi),
            int(metrics.get("ratio_fail_count", 10**9)),
            _safe_float(metrics.get("ratio_mean_dev"), default=float("inf")),
            _safe_float(metrics.get("ratio_rms"), default=float("inf")),
            _safe_float(metrics.get("kinematic_score"), default=float("inf")),
            -int(metrics.get("valid_ratio_bins", 0)),
        )
    return _candidate_sort_key(result)


def _is_better_result(candidate, incumbent):
    if candidate is None:
        return False
    if incumbent is None:
        return True
    if not candidate.get("valid"):
        return False
    if not incumbent.get("valid"):
        return True
    return _candidate_order_key(candidate) < _candidate_order_key(incumbent)


def _result_summary(result):
    if not result:
        return "none"
    metrics = result.get("metrics", {})
    score = result.get("selection_score")
    score_str = ""
    if score is not None and math.isfinite(_safe_float(score, default=float("inf"))):
        score_str = "score={} ".format(_format_metric_value(score))
    return (
        "{}BG_STAT_SCALE1={:.3f} BG_STAT_SCALE2={:.3f} fail={} mean_dev={} rms={} kin={}".format(
            score_str,
            float(result.get("bg_scale1", BG_STAT_SCALE1)),
            float(result.get("bg_scale2", BG_STAT_SCALE2)),
            int(metrics.get("ratio_fail_count", 10**9)),
            _format_metric_value(metrics.get("ratio_mean_dev", float("inf"))),
            _format_metric_value(metrics.get("ratio_rms", float("inf"))),
            _format_metric_value(metrics.get("kinematic_score", float("inf"))),
        )
    )


def _get_profile_metadata():
    resolved = get_resolved_bg_profile_settings()
    return {
        "active_profile": resolved.get("name"),
        "resolved_profile": resolved,
        "forced_bg_scale1": get_forced_bg_scale1(),
        "forced_bg_scale2": get_forced_bg_scale2(),
        "use_common_epsilon_scales": use_common_epsilon_scales(),
        "common_epsilon_scale_behavior": get_common_epsilon_scale_behavior(),
    }


def _copy_static_hist_context(base_hist, candidate_hist):
    candidate_hist.update({key: val for key, val in base_hist.items() if key not in candidate_hist})
    return candidate_hist


def _prepare_candidate_hist(base_hist, t_bins, phi_bins):
    candidate_hist = dict(base_hist)
    candidate_hist["t_bins"] = np.array(t_bins, dtype=float)
    candidate_hist["phi_bins"] = np.array(phi_bins, dtype=float)
    return candidate_hist


def _build_ratio_metrics(ratio_dict, phi_setting):
    ratio_entries = (
        ratio_dict.get("binned", {})
        .get(phi_setting, {})
        .get("ratio", {})
    )

    valid = []
    fail_count = 0
    for entry in ratio_entries.values():
        ratio = _safe_float(entry.get("ratio"))
        ratio_err = _safe_float(entry.get("ratio_err"))
        if ratio_err <= 0.0 or not math.isfinite(ratio):
            continue
        sigma = abs(ratio - 1.0) / ratio_err
        if sigma > RATIO_SIGMA_THRESHOLD:
            fail_count += 1
        valid.append((ratio, ratio_err))

    if not valid:
        return {
            "ratio_fail_count": int(10**9),
            "ratio_mean_dev": float("inf"),
            "ratio_rms": float("inf"),
            "valid_ratio_bins": 0,
        }

    ratios = np.array([val[0] for val in valid], dtype=float)
    errors = np.array([val[1] for val in valid], dtype=float)
    weights = 1.0 / np.square(errors)
    weighted_mean = np.sum(weights * ratios) / np.sum(weights)
    rms = math.sqrt(np.mean(np.square(ratios - 1.0)))

    return {
        "ratio_fail_count": int(fail_count),
        "ratio_mean_dev": float(abs(weighted_mean - 1.0)),
        "ratio_rms": float(rms),
        "valid_ratio_bins": int(len(valid)),
    }


def _hist_chi2_ndf(data_hist, simc_hist):
    if data_hist is None or simc_hist is None:
        return float("inf")
    if data_hist.GetNbinsX() != simc_hist.GetNbinsX():
        return float("inf")

    chi2 = 0.0
    n_used = 0
    for idx in range(1, data_hist.GetNbinsX() + 1):
        data_val = float(data_hist.GetBinContent(idx))
        simc_val = float(simc_hist.GetBinContent(idx))
        variance = abs(data_val) + abs(simc_val)
        if variance <= 0.0:
            continue
        chi2 += ((data_val - simc_val) ** 2) / variance
        n_used += 1

    if n_used <= 1:
        return float("inf")
    return chi2 / float(n_used - 1)


def _relative_moment_difference(data_hist, simc_hist, getter_name):
    data_val = float(getattr(data_hist, getter_name)())
    simc_val = float(getattr(simc_hist, getter_name)())
    denom = max(abs(data_val), abs(simc_val), 1.0e-12)
    return abs(data_val - simc_val) / denom


def _build_kinematic_metrics(hist):
    data_support = hist.get("_xsect_support_data")
    simc_support = hist.get("_xsect_support_simc")
    if not data_support or not simc_support:
        return {
            "kinematic_score": float("inf"),
            "kinematic_points": 0,
        }

    scores = []
    for var_name in KINEMATIC_SCORE_VARS:
        data_matrix = data_support.get(var_name, [])
        simc_matrix = simc_support.get(var_name, [])
        for data_row, simc_row in zip(data_matrix, simc_matrix):
            for data_hist, simc_hist in zip(data_row, simc_row):
                if data_hist is None or simc_hist is None:
                    continue
                chi2_ndf = _hist_chi2_ndf(data_hist, simc_hist)
                mean_rel = _relative_moment_difference(data_hist, simc_hist, "GetMean")
                rms_rel = _relative_moment_difference(data_hist, simc_hist, "GetRMS")
                if all(math.isfinite(val) for val in (chi2_ndf, mean_rel, rms_rel)):
                    scores.append(chi2_ndf + mean_rel + rms_rel)

    if not scores:
        return {
            "kinematic_score": float("inf"),
            "kinematic_points": 0,
        }

    return {
        "kinematic_score": float(np.mean(scores)),
        "kinematic_points": int(len(scores)),
    }


def _build_candidate_metrics(phi_setting, ratio_dict, hist):
    metrics = {}
    metrics.update(_build_ratio_metrics(ratio_dict, phi_setting))
    metrics.update(_build_kinematic_metrics(hist))
    return metrics


def _log_scale_result(stage_name, epsset, phi_setting, result):
    if not result.get("valid"):
        _log(
            "{} result {} {} BG_STAT_SCALE1={:.3f} BG_STAT_SCALE2={:.3f} -> invalid (valid_ratio_bins={} kin_points={})".format(
                stage_name,
                epsset,
                phi_setting,
                float(result.get("bg_scale1", BG_STAT_SCALE1)),
                float(result.get("bg_scale2", BG_STAT_SCALE2)),
                int(result.get("metrics", {}).get("valid_ratio_bins", 0)),
                int(result.get("metrics", {}).get("kinematic_points", 0)),
            )
        )
        return

    metrics = result["metrics"]
    score_str = ""
    if result.get("selection_score") is not None:
        score_str = " score={}".format(_format_metric_value(result.get("selection_score")))
    _log(
        "{} result {} {} BG_STAT_SCALE1={:.3f} BG_STAT_SCALE2={:.3f} -> fail={} mean_dev={} rms={} kin={} bins={} kin_pts={}{}".format(
            stage_name,
            epsset,
            phi_setting,
            float(result["bg_scale1"]),
            float(result["bg_scale2"]),
            int(metrics["ratio_fail_count"]),
            _format_metric_value(metrics["ratio_mean_dev"]),
            _format_metric_value(metrics["ratio_rms"]),
            _format_metric_value(metrics["kinematic_score"]),
            int(metrics["valid_ratio_bins"]),
            int(metrics["kinematic_points"]),
            score_str,
        )
    )


def _evaluate_scale_grid(
    base_hist,
    inpDict,
    phi_setting,
    t_bins,
    phi_bins,
    simc_yield_dict,
    simc_support,
    data_base_cache,
    stage_name,
    scale_axis,
    bg_scales,
    fixed_bg_scale1,
    fixed_bg_scale2,
):
    results = []
    total = len(bg_scales)
    running_best = None
    scale_label = scale_axis.upper()
    for idx, bg_scale in enumerate(bg_scales, start=1):
        bg_scale1 = float(bg_scale) if scale_axis == "bg_stat_scale1" else float(fixed_bg_scale1)
        bg_scale2 = float(bg_scale) if scale_axis == "bg_stat_scale2" else float(fixed_bg_scale2)
        _log(
            "{} scan {} {}: testing {}={:.3f} with BG_STAT_SCALE1={:.3f} BG_STAT_SCALE2={:.3f} ({}/{})".format(
                stage_name,
                inpDict["EPSSET"],
                phi_setting,
                scale_label,
                float(bg_scale),
                bg_scale1,
                bg_scale2,
                idx,
                total,
            )
        )
        result = _evaluate_phi_candidate(
            base_hist,
            inpDict,
            phi_setting,
            bg_scale1,
            bg_scale2,
            t_bins,
            phi_bins,
            simc_yield_dict,
            simc_support,
            data_base_cache,
        )
        result["scale_axis"] = scale_axis
        results.append(result)
        _annotate_selection_scores(results)
        _log_scale_result(stage_name, inpDict["EPSSET"], phi_setting, result)
        previous_best = running_best
        valid_results = [candidate for candidate in results if candidate.get("valid")]
        if valid_results:
            valid_results.sort(key=_candidate_order_key)
            running_best = valid_results[0]
            if previous_best is None or running_best is not previous_best:
                _log(
                    "{} leader update {} {} after {}/{}: {}".format(
                        stage_name,
                        inpDict["EPSSET"],
                        phi_setting,
                        idx,
                        total,
                        _result_summary(running_best),
                    )
                )
            else:
                _log(
                    "{} leader unchanged {} {} after {}/{}: {}".format(
                        stage_name,
                        inpDict["EPSSET"],
                        phi_setting,
                        idx,
                        total,
                        _result_summary(running_best),
                    )
                )
        elif running_best is not None:
            _log(
                "{} leader unchanged {} {} after {}/{}: {}".format(
                    stage_name,
                    inpDict["EPSSET"],
                    phi_setting,
                    idx,
                    total,
                    _result_summary(running_best),
                )
            )
        else:
            _log(
                "{} leader pending {} {} after {}/{}: no valid {} yet".format(
                    stage_name,
                    inpDict["EPSSET"],
                    phi_setting,
                    idx,
                    total,
                    scale_label,
                )
            )
    return results


def _evaluate_fixed_result(
    base_hist,
    inpDict,
    phi_setting,
    t_bins,
    phi_bins,
    simc_yield_dict,
    simc_support,
    data_base_cache,
    bg_scale1,
    bg_scale2,
    stage_name,
):
    result = _evaluate_phi_candidate(
        base_hist,
        inpDict,
        phi_setting,
        float(bg_scale1),
        float(bg_scale2),
        t_bins,
        phi_bins,
        simc_yield_dict,
        simc_support,
        data_base_cache,
    )
    result["selection_mode"] = _get_selection_mode()
    result["scan_skipped"] = True
    result["scan_skip_reason"] = stage_name
    return result


def _evaluate_phi_candidate(base_hist, inpDict, phi_setting, bg_scale1, bg_scale2, t_bins, phi_bins, simc_yield_dict, simc_support, data_base_cache):
    candidate_inp = dict(inpDict)
    candidate_inp["NumtBins"] = len(t_bins) - 1
    candidate_inp["NumPhiBins"] = len(phi_bins) - 1
    candidate_inp["yield_emit_plots"] = False
    candidate_inp["yield_show_progress"] = False
    candidate_inp["suppress_bg_opt_warnings"] = True
    candidate_inp["bg_opt_use_data_cache"] = True
    candidate_inp["bg_stat_scale1"] = float(bg_scale1)
    candidate_inp["bg_stat_scale2"] = float(bg_scale2)
    candidate_inp["bg_stat_scale1_by_setting"] = {
        get_bg_scale_setting_key(candidate_inp["EPSSET"], phi_setting): float(bg_scale1)
    }
    candidate_inp["bg_stat_scale2_by_setting"] = {
        get_bg_scale_setting_key(candidate_inp["EPSSET"], phi_setting): float(bg_scale2)
    }

    sys.path.append("binning")
    from calculate_yield import find_yield_data
    from calculate_ratio import find_ratio

    candidate_hist = _prepare_candidate_hist(base_hist, t_bins, phi_bins)
    candidate_hist["_xsect_support_simc"] = simc_support
    candidate_hist["_bg_opt_data_base_cache"] = data_base_cache

    yield_dict = {}
    yield_dict.update(find_yield_data([candidate_hist], candidate_inp))
    yield_dict.update(simc_yield_dict)

    ratio_dict = {}
    ratio_dict.update(find_ratio([candidate_hist], candidate_inp, yield_dict))

    metrics = _build_candidate_metrics(phi_setting, ratio_dict, candidate_hist)
    return {
        "phi_setting": phi_setting,
        "bg_scale1": float(bg_scale1),
        "bg_scale2": float(bg_scale2),
        "valid": all(math.isfinite(val) for val in (
            metrics["ratio_mean_dev"],
            metrics["ratio_rms"],
            metrics["kinematic_score"],
        )) and metrics["valid_ratio_bins"] > 0,
        "metrics": metrics,
        "hist": candidate_hist,
    }


def _build_simc_reference(base_hist, inpDict, t_bins, phi_bins):
    candidate_inp = dict(inpDict)
    candidate_inp["NumtBins"] = len(t_bins) - 1
    candidate_inp["NumPhiBins"] = len(phi_bins) - 1
    candidate_inp["yield_emit_plots"] = False
    candidate_inp["yield_show_progress"] = False
    candidate_inp["suppress_bg_opt_warnings"] = True

    sys.path.append("binning")
    from calculate_yield import find_yield_simc

    simc_hist = _prepare_candidate_hist(base_hist, t_bins, phi_bins)
    simc_yield_dict = {}
    simc_yield_dict.update(find_yield_simc([simc_hist], candidate_inp))
    return simc_yield_dict, simc_hist.get("_xsect_support_simc")


def _parse_test_model_chi2(param_path):
    chi2_vals = []
    if not os.path.exists(param_path):
        return None
    with open(param_path, "r") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 4:
                continue
            chi2 = _safe_float(parts[3], default=float("nan"))
            if math.isfinite(chi2) and chi2 > 0.0:
                chi2_vals.append(abs(chi2 - 1.0))
    if not chi2_vals:
        return None
    return float(np.mean(chi2_vals))


def _run_test_model_tiebreak(inpDict):
    ltanapath = inpDict.get("LTANAPATH")
    outpath = inpDict.get("OUTPATH")
    if not ltanapath or not outpath:
        return None

    q2 = inpDict["Q2"]
    w = inpDict["W"]
    tmin = str(inpDict["tmin"])
    tmax = str(inpDict["tmax"])
    particle_type = inpDict["ParticleType"]
    pol = str(inpDict["POL"])
    debug = "False"
    diag_date = "{}_bgopt".format(inpDict.get("formatted_date", "bgopt"))
    outfilename = "Testing_Q{}W{}_bgopt".format(q2, w)

    script_path = os.path.join(ltanapath, "test_model.sh")
    if not os.path.exists(script_path):
        return None

    try:
        subprocess.run(
            [
                "bash",
                script_path,
                q2,
                w,
                tmin,
                tmax,
                particle_type,
                pol,
                outfilename,
                diag_date,
                debug,
            ],
            check=True,
            cwd=ltanapath,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except Exception:
        return None

    q2_num = q2.replace("p", "")
    w_num = w.replace("p", "")
    if len(q2_num) < 2:
        q2_num = q2_num.ljust(2, "0")
    if len(w_num) < 3:
        w_num = w_num.ljust(3, "0")
    pol_str = "pl" if int(pol) == 1 else "mn"
    param_path = os.path.join(
        outpath,
        "testing_env",
        particle_type.lower(),
        "Q{}W{}".format(q2, w),
        diag_date,
        "new_par.{}_Q{}W{}.dat".format(pol_str, q2_num, w_num),
    )
    return _parse_test_model_chi2(param_path)


def _finalize_phi_results(results, inpDict):
    _annotate_selection_scores(results)
    valid_results = [result for result in results if result.get("valid")]
    if not valid_results:
        return {
            "phi_setting": results[0]["phi_setting"] if results else "unknown",
            "bg_scale1": float(BG_STAT_SCALE1),
            "bg_scale2": float(BG_STAT_SCALE2),
            "valid": False,
            "fallback": True,
            "metrics": {
                "ratio_fail_count": int(10**9),
                "ratio_mean_dev": float("inf"),
                "ratio_rms": float("inf"),
                "kinematic_score": float("inf"),
                "valid_ratio_bins": 0,
                "kinematic_points": 0,
            },
            "results": results,
        }

    valid_results.sort(key=_candidate_order_key)
    finalists = valid_results[: max(1, BG_STAT_SCALE2_FINALIST_COUNT)]
    if len(finalists) > 1:
        top_key = _candidate_order_key(finalists[0])
        tied = [result for result in finalists if _candidate_order_key(result) == top_key]
        if len(tied) > 1:
            _log(
                "Running test_model.sh tie-break for {} {} among BG_STAT_SCALE2={}".format(
                    inpDict["EPSSET"],
                    results[0]["phi_setting"],
                    ", ".join(
                        "({:.3f}, {:.3f})".format(
                            float(result["bg_scale1"]),
                            float(result["bg_scale2"]),
                        )
                        for result in tied
                    ),
                )
            )
            tiebreak_chi2 = _run_test_model_tiebreak(inpDict)
            if tiebreak_chi2 is not None:
                tied[0]["test_model_chi2"] = tiebreak_chi2
                finalists = tied

    best_result = dict(valid_results[0])
    if "metrics" in best_result and isinstance(best_result["metrics"], dict):
        best_result["metrics"] = dict(best_result["metrics"])
    best_result["results"] = results
    return best_result


def _optimize_phi_scale(base_hist, inpDict, t_bins, phi_bins):
    phi_setting = base_hist["phi_setting"]
    _log(
        "Scanning {} {} with NumtBins={} NumPhiBins={}".format(
            inpDict["EPSSET"],
            phi_setting,
            len(t_bins) - 1,
            len(phi_bins) - 1,
        )
    )
    sys.path.append("binning")
    from calculate_yield import prepare_bg_opt_data_base_cache

    simc_yield_dict, simc_support = _build_simc_reference(base_hist, inpDict, t_bins, phi_bins)
    data_base_cache = prepare_bg_opt_data_base_cache(base_hist, inpDict, t_bins, phi_bins)
    initial_scale1 = resolve_bg_stat_scale1(inpDict, phi_setting)
    initial_scale2 = resolve_bg_stat_scale2(inpDict, phi_setting)
    forced_scale1 = get_forced_bg_scale1()
    forced_scale2 = get_forced_bg_scale2()
    scan_skipped = {
        "bg_stat_scale1": forced_scale1 is not None,
        "bg_stat_scale2": forced_scale2 is not None,
    }

    if forced_scale1 is not None and forced_scale2 is not None:
        fixed_result = _evaluate_fixed_result(
            base_hist,
            inpDict,
            phi_setting,
            t_bins,
            phi_bins,
            simc_yield_dict,
            simc_support,
            data_base_cache,
            forced_scale1,
            forced_scale2,
            "profile_forced_both_scales",
        )
        fixed_result["coarse_fit1_best"] = fixed_result
        fixed_result["coarse_fit1_results"] = [fixed_result]
        fixed_result["coarse_fit2_best"] = fixed_result
        fixed_result["coarse_fit2_results"] = [fixed_result]
        fixed_result["refined_fit1_best"] = fixed_result
        fixed_result["refined_fit1_results"] = [fixed_result]
        fixed_result["refined_fit2_best"] = fixed_result
        fixed_result["refined_fit2_results"] = [fixed_result]
        fixed_result["forced_bg_scale1"] = forced_scale1
        fixed_result["forced_bg_scale2"] = forced_scale2
        fixed_result["scan_skipped"] = scan_skipped
        return fixed_result

    if forced_scale1 is not None:
        coarse_candidates1 = [float(forced_scale1)]
        _log(
            "Skipping BG_STAT_SCALE1 scan for {} {} due to active profile; fixed at {:.3f}".format(
                inpDict["EPSSET"],
                phi_setting,
                float(forced_scale1),
            )
        )
    else:
        coarse_candidates1 = list(get_bg_scale1_coarse_candidates())
        _log(
            "Coarse BG_STAT_SCALE1 grid for {} {}: {}".format(
                inpDict["EPSSET"],
                phi_setting,
                ", ".join("{:.3f}".format(float(scale)) for scale in coarse_candidates1),
            )
        )
    coarse_fit1_results = _evaluate_scale_grid(
        base_hist,
        inpDict,
        phi_setting,
        t_bins,
        phi_bins,
        simc_yield_dict,
        simc_support,
        data_base_cache,
        "Coarse",
        "bg_stat_scale1",
        coarse_candidates1,
        forced_scale1 if forced_scale1 is not None else initial_scale1,
        forced_scale2 if forced_scale2 is not None else initial_scale2,
    )

    best_coarse_fit1 = _finalize_phi_results(coarse_fit1_results, inpDict)
    if not best_coarse_fit1.get("valid"):
        _log("No valid coarse BG_STAT_SCALE1 found for {} {}".format(inpDict["EPSSET"], phi_setting))
        best_coarse_fit1["forced_bg_scale1"] = forced_scale1
        best_coarse_fit1["forced_bg_scale2"] = forced_scale2
        best_coarse_fit1["scan_skipped"] = scan_skipped
        return best_coarse_fit1

    _log(
        "Best coarse {} {} BG_STAT_SCALE1={:.3f} with BG_STAT_SCALE2={:.3f} fail={} mean_dev={:.4f} rms={:.4f}".format(
            inpDict["EPSSET"],
            phi_setting,
            float(best_coarse_fit1["bg_scale1"]),
            float(best_coarse_fit1["bg_scale2"]),
            int(best_coarse_fit1["metrics"]["ratio_fail_count"]),
            float(best_coarse_fit1["metrics"]["ratio_mean_dev"]),
            float(best_coarse_fit1["metrics"]["ratio_rms"]),
        )
    )

    if forced_scale2 is not None:
        coarse_candidates2 = [float(forced_scale2)]
        _log(
            "Skipping BG_STAT_SCALE2 scan for {} {} due to active profile; fixed at {:.3f}".format(
                inpDict["EPSSET"],
                phi_setting,
                float(forced_scale2),
            )
        )
    else:
        coarse_candidates2 = list(get_bg_scale_coarse_candidates())
        _log(
            "Coarse BG_STAT_SCALE2 grid for {} {}: {}".format(
                inpDict["EPSSET"],
                phi_setting,
                ", ".join("{:.3f}".format(float(scale)) for scale in coarse_candidates2),
            )
        )
    coarse_fit2_results = _evaluate_scale_grid(
        base_hist,
        inpDict,
        phi_setting,
        t_bins,
        phi_bins,
        simc_yield_dict,
        simc_support,
        data_base_cache,
        "Coarse",
        "bg_stat_scale2",
        coarse_candidates2,
        best_coarse_fit1["bg_scale1"],
        forced_scale2 if forced_scale2 is not None else initial_scale2,
    )

    best_coarse_fit2 = _finalize_phi_results(coarse_fit2_results, inpDict)
    if not best_coarse_fit2.get("valid"):
        best_coarse_fit2 = best_coarse_fit1

    _log(
        "Best coarse {} {} BG_STAT_SCALE2={:.3f} with BG_STAT_SCALE1={:.3f} fail={} mean_dev={:.4f} rms={:.4f}".format(
            inpDict["EPSSET"],
            phi_setting,
            float(best_coarse_fit2["bg_scale2"]),
            float(best_coarse_fit2["bg_scale1"]),
            int(best_coarse_fit2["metrics"]["ratio_fail_count"]),
            float(best_coarse_fit2["metrics"]["ratio_mean_dev"]),
            float(best_coarse_fit2["metrics"]["ratio_rms"]),
        )
    )

    if forced_scale1 is not None:
        refined_fit1_results = list(coarse_fit1_results)
        best_refined_fit1 = dict(best_coarse_fit1)
    else:
        refined_candidates1 = list(get_bg_scale1_refined_candidates(best_coarse_fit1["bg_scale1"]))
        _log(
            "Refined BG_STAT_SCALE1 grid for {} {}: {}".format(
                inpDict["EPSSET"],
                phi_setting,
                ", ".join("{:.3f}".format(float(scale)) for scale in refined_candidates1),
            )
        )
        refined_fit1_results = _evaluate_scale_grid(
            base_hist,
            inpDict,
            phi_setting,
            t_bins,
            phi_bins,
            simc_yield_dict,
            simc_support,
            data_base_cache,
            "Refined",
            "bg_stat_scale1",
            refined_candidates1,
            best_coarse_fit1["bg_scale1"],
            best_coarse_fit2["bg_scale2"],
        )
        best_refined_fit1 = _finalize_phi_results(refined_fit1_results, inpDict)
        if not best_refined_fit1.get("valid"):
            best_refined_fit1 = best_coarse_fit2

    if forced_scale2 is not None:
        refined_fit2_results = list(coarse_fit2_results)
        best_refined_fit2 = dict(best_coarse_fit2)
    else:
        refined_candidates2 = list(get_bg_scale_refined_candidates(best_coarse_fit2["bg_scale2"]))
        _log(
            "Refined BG_STAT_SCALE2 grid for {} {}: {}".format(
                inpDict["EPSSET"],
                phi_setting,
                ", ".join("{:.3f}".format(float(scale)) for scale in refined_candidates2),
            )
        )
        refined_fit2_results = _evaluate_scale_grid(
            base_hist,
            inpDict,
            phi_setting,
            t_bins,
            phi_bins,
            simc_yield_dict,
            simc_support,
            data_base_cache,
            "Refined",
            "bg_stat_scale2",
            refined_candidates2,
            best_refined_fit1["bg_scale1"],
            best_coarse_fit2["bg_scale2"],
        )
        best_refined_fit2 = _finalize_phi_results(refined_fit2_results, inpDict)

    selected_result = None
    for candidate in (best_refined_fit2, best_refined_fit1, best_coarse_fit2, best_coarse_fit1):
        if candidate.get("valid"):
            selected_result = dict(candidate)
            break
    if selected_result is None:
        selected_result = dict(best_refined_fit2)

    if "metrics" in selected_result and isinstance(selected_result["metrics"], dict):
        selected_result["metrics"] = dict(selected_result["metrics"])
    selected_result["coarse_fit1_best"] = best_coarse_fit1
    selected_result["coarse_fit1_results"] = coarse_fit1_results
    selected_result["coarse_fit2_best"] = best_coarse_fit2
    selected_result["coarse_fit2_results"] = coarse_fit2_results
    selected_result["refined_fit1_best"] = best_refined_fit1
    selected_result["refined_fit1_results"] = refined_fit1_results
    selected_result["refined_fit2_best"] = best_refined_fit2
    selected_result["refined_fit2_results"] = refined_fit2_results
    selected_result["forced_bg_scale1"] = forced_scale1
    selected_result["forced_bg_scale2"] = forced_scale2
    selected_result["scan_skipped"] = scan_skipped

    if selected_result.get("valid"):
        _log(
            "Selected {} {} BG_STAT_SCALE1={:.3f} BG_STAT_SCALE2={:.3f} fail={} mean_dev={:.4f} rms={:.4f} kin={:.4f}".format(
                inpDict["EPSSET"],
                phi_setting,
                float(selected_result["bg_scale1"]),
                float(selected_result["bg_scale2"]),
                int(selected_result["metrics"]["ratio_fail_count"]),
                float(selected_result["metrics"]["ratio_mean_dev"]),
                float(selected_result["metrics"]["ratio_rms"]),
                float(selected_result["metrics"]["kinematic_score"]),
            )
        )
    return selected_result


def _aggregate_bin_candidate_result(phi_results):
    metrics = {
        "ratio_fail_count": 0,
        "ratio_mean_dev": 0.0,
        "ratio_rms": 0.0,
        "kinematic_score": 0.0,
        "valid_ratio_bins": 0,
        "kinematic_points": 0,
    }

    if not phi_results:
        metrics.update({
            "ratio_fail_count": int(10**9),
            "ratio_mean_dev": float("inf"),
            "ratio_rms": float("inf"),
            "kinematic_score": float("inf"),
        })
        return metrics

    valid_count = 0
    for result in phi_results:
        if not result.get("valid"):
            metrics["ratio_fail_count"] = int(10**9)
            metrics["ratio_mean_dev"] = float("inf")
            metrics["ratio_rms"] = float("inf")
            metrics["kinematic_score"] = float("inf")
            return metrics
        valid_count += 1
        metrics["ratio_fail_count"] += int(result["metrics"]["ratio_fail_count"])
        metrics["ratio_mean_dev"] += float(result["metrics"]["ratio_mean_dev"])
        metrics["ratio_rms"] += float(result["metrics"]["ratio_rms"])
        metrics["kinematic_score"] += float(result["metrics"]["kinematic_score"])
        metrics["valid_ratio_bins"] += int(result["metrics"]["valid_ratio_bins"])
        metrics["kinematic_points"] += int(result["metrics"]["kinematic_points"])

    metrics["ratio_mean_dev"] /= float(valid_count)
    metrics["ratio_rms"] /= float(valid_count)
    metrics["kinematic_score"] /= float(valid_count)
    return metrics


def _serialize_result(result, seen=None):
    if seen is None:
        seen = set()

    result_id = id(result)
    if isinstance(result, (dict, list, tuple)) and result_id in seen:
        return "<cycle>"

    if isinstance(result, dict):
        seen.add(result_id)
        serializable = {}
        for key, val in result.items():
            if key == "hist":
                continue
            serializable[key] = _serialize_result(val, seen)
        seen.discard(result_id)
        return serializable
    if isinstance(result, list):
        seen.add(result_id)
        serializable = [_serialize_result(entry, seen) for entry in result]
        seen.discard(result_id)
        return serializable
    if isinstance(result, tuple):
        seen.add(result_id)
        serializable = [_serialize_result(entry, seen) for entry in result]
        seen.discard(result_id)
        return serializable
    if isinstance(result, np.ndarray):
        return result.tolist()
    if isinstance(result, np.generic):
        return result.item()
    return result


def _build_summary_report(inpDict, mode, proposal, phi_results):
    profile_meta = _get_profile_metadata()
    lines = [
        "BG/Tuning Summary ({})".format(mode),
        "Active profile: {}".format(profile_meta["active_profile"]),
        "Selection mode: {}".format(_get_selection_mode()),
        "Use common epsilon scales: {}".format(profile_meta["use_common_epsilon_scales"]),
        "Common epsilon scale behavior: {}".format(profile_meta["common_epsilon_scale_behavior"]),
        "Shared bins: NumtBins={} NumPhiBins={}".format(
            len(proposal["t_bins"]) - 1,
            len(proposal["phi_bins"]) - 1,
        ),
    ]
    for result in phi_results:
        status = "fallback" if result.get("fallback") else "selected"
        lines.append(
            "{} {}: BG_STAT_SCALE1={:.3f} BG_STAT_SCALE2={:.3f} score={}".format(
                result.get("phi_setting", "unknown"),
                status,
                float(result.get("bg_scale1", BG_STAT_SCALE1)),
                float(result.get("bg_scale2", BG_STAT_SCALE2)),
                _format_metric_value(result.get("selection_score")),
            )
        )
    return "\n".join(lines)


def _metric_columns(metrics):
    metric_dict = metrics or {}
    return {
        "ratio_fail_count": metric_dict.get("ratio_fail_count"),
        "ratio_mean_dev": metric_dict.get("ratio_mean_dev"),
        "ratio_rms": metric_dict.get("ratio_rms"),
        "kinematic_score": metric_dict.get("kinematic_score"),
        "valid_ratio_bins": metric_dict.get("valid_ratio_bins"),
        "kinematic_points": metric_dict.get("kinematic_points"),
    }


def _base_csv_row(
    inpDict,
    summary,
    row_kind,
    proposal=None,
    phi_setting="",
    stage="",
    bg_scale1=None,
    bg_scale2=None,
    scale_axis="",
    valid=None,
):
    proposal = proposal or {}
    row = {
        "mode": summary.get("mode", ""),
        "active_profile": summary.get("active_profile"),
        "selection_mode": summary.get("selection_mode", _get_selection_mode()),
        "epsset": inpDict.get("EPSSET", ""),
        "particle": inpDict.get("ParticleType", ""),
        "q2": inpDict.get("Q2", ""),
        "w": inpDict.get("W", ""),
        "outfilename": inpDict.get("OutFilename", ""),
        "row_kind": row_kind,
        "stage": stage,
        "scale_axis": scale_axis,
        "phi_setting": phi_setting,
        "bg_stat_scale1": bg_scale1,
        "bg_stat_scale2": bg_scale2,
        "valid": valid,
        "fallback": summary.get("fallback", False),
        "forced_bg_scale1": summary.get("forced_bg_scale1"),
        "forced_bg_scale2": summary.get("forced_bg_scale2"),
        "use_common_epsilon_scales": summary.get("use_common_epsilon_scales"),
        "common_epsilon_scale_behavior": summary.get("common_epsilon_scale_behavior"),
        "resolved_profile_json": json.dumps(summary.get("resolved_profile", {}), sort_keys=True),
        "requested_num_t_bins": proposal.get("requested_num_t_bins"),
        "requested_num_phi_bins": proposal.get("requested_num_phi_bins"),
        "actual_num_t_bins": proposal.get("actual_num_t_bins"),
        "actual_num_phi_bins": proposal.get("actual_num_phi_bins"),
    }
    return row


def _proposal_matches(selected_proposal, candidate_proposal):
    if selected_proposal is candidate_proposal:
        return True
    if not selected_proposal or not candidate_proposal:
        return False
    try:
        return (
            np.array_equal(
                np.asarray(selected_proposal.get("t_bins", []), dtype=float),
                np.asarray(candidate_proposal.get("t_bins", []), dtype=float),
            )
            and np.array_equal(
                np.asarray(selected_proposal.get("phi_bins", []), dtype=float),
                np.asarray(candidate_proposal.get("phi_bins", []), dtype=float),
            )
        )
    except Exception:
        return False


def _scales_match(result, selected_result):
    if not result or not selected_result:
        return False
    try:
        return (
            math.isclose(float(result.get("bg_scale1", float("nan"))), float(selected_result.get("bg_scale1", float("nan"))), rel_tol=0.0, abs_tol=1.0e-9)
            and math.isclose(float(result.get("bg_scale2", float("nan"))), float(selected_result.get("bg_scale2", float("nan"))), rel_tol=0.0, abs_tol=1.0e-9)
        )
    except Exception:
        return False


def _append_scale_rows(rows, inpDict, summary, proposal, selected_result, results, stage_name, scale_axis):
    for result in results or []:
        row = _base_csv_row(
            inpDict,
            summary,
            row_kind="scale_candidate",
            proposal=proposal,
            phi_setting=result.get("phi_setting", selected_result.get("phi_setting", "")),
            stage=stage_name,
            bg_scale1=result.get("bg_scale1"),
            bg_scale2=result.get("bg_scale2"),
            scale_axis=scale_axis,
            valid=result.get("valid"),
        )
        row.update(_metric_columns(result.get("metrics")))
        row["selection_score"] = result.get("selection_score")
        row["selected_for_stage"] = result.get("valid") and _scales_match(result, selected_result)
        row["test_model_chi2"] = result.get("test_model_chi2")
        rows.append(row)


def build_optimization_csv_rows(summary, inpDict):
    rows = []
    mode = summary.get("mode")

    if mode == "low":
        for candidate in summary.get("candidate_results", []):
            proposal = candidate.get("proposal", {})
            candidate_row = _base_csv_row(
                inpDict,
                summary,
                row_kind="bin_candidate",
                proposal=proposal,
                stage="aggregate",
                valid=candidate.get("valid"),
            )
            candidate_row.update(_metric_columns(candidate.get("metrics")))
            candidate_row["selection_score"] = candidate.get("selection_score")
            candidate_row["selected_bin_candidate"] = bool(
                not summary.get("fallback")
                and _proposal_matches(summary.get("proposal"), proposal)
            )
            candidate_row["error"] = candidate.get("error")
            rows.append(candidate_row)

            for phi_result in candidate.get("phi_results", []):
                _append_scale_rows(rows, inpDict, summary, proposal, phi_result.get("coarse_fit1_best", {}), phi_result.get("coarse_fit1_results", []), "coarse", "bg_stat_scale1")
                _append_scale_rows(rows, inpDict, summary, proposal, phi_result.get("coarse_fit2_best", {}), phi_result.get("coarse_fit2_results", []), "coarse", "bg_stat_scale2")
                _append_scale_rows(rows, inpDict, summary, proposal, phi_result.get("refined_fit1_best", {}), phi_result.get("refined_fit1_results", []), "refined", "bg_stat_scale1")
                _append_scale_rows(rows, inpDict, summary, proposal, phi_result.get("refined_fit2_best", {}), phi_result.get("refined_fit2_results", []), "refined", "bg_stat_scale2")

                selected_row = _base_csv_row(
                    inpDict,
                    summary,
                    row_kind="phi_selection",
                    proposal=proposal,
                    phi_setting=phi_result.get("phi_setting", ""),
                    stage="selected",
                    bg_scale1=phi_result.get("bg_scale1"),
                    bg_scale2=phi_result.get("bg_scale2"),
                    valid=phi_result.get("valid"),
                )
                selected_row.update(_metric_columns(phi_result.get("metrics")))
                selected_row["selection_score"] = phi_result.get("selection_score")
                selected_row["selected_for_stage"] = bool(
                    not summary.get("fallback")
                    and _proposal_matches(summary.get("proposal"), proposal)
                )
                selected_row["test_model_chi2"] = phi_result.get("test_model_chi2")
                rows.append(selected_row)

    elif mode == "high":
        proposal = summary.get("proposal", {})
        aggregate_row = _base_csv_row(
            inpDict,
            summary,
            row_kind="bin_candidate",
            proposal=proposal,
            stage="aggregate",
            valid=not summary.get("fallback", False),
        )
        aggregate_row.update(_metric_columns(summary.get("metrics")))
        aggregate_row["selection_score"] = summary.get("selection_score")
        aggregate_row["selected_bin_candidate"] = True
        rows.append(aggregate_row)

        for phi_result in summary.get("candidate_results", []):
            _append_scale_rows(rows, inpDict, summary, proposal, phi_result.get("coarse_fit1_best", {}), phi_result.get("coarse_fit1_results", []), "coarse", "bg_stat_scale1")
            _append_scale_rows(rows, inpDict, summary, proposal, phi_result.get("coarse_fit2_best", {}), phi_result.get("coarse_fit2_results", []), "coarse", "bg_stat_scale2")
            _append_scale_rows(rows, inpDict, summary, proposal, phi_result.get("refined_fit1_best", {}), phi_result.get("refined_fit1_results", []), "refined", "bg_stat_scale1")
            _append_scale_rows(rows, inpDict, summary, proposal, phi_result.get("refined_fit2_best", {}), phi_result.get("refined_fit2_results", []), "refined", "bg_stat_scale2")

            selected_row = _base_csv_row(
                inpDict,
                summary,
                row_kind="phi_selection",
                proposal=proposal,
                phi_setting=phi_result.get("phi_setting", ""),
                stage="selected",
                bg_scale1=phi_result.get("bg_scale1"),
                bg_scale2=phi_result.get("bg_scale2"),
                valid=phi_result.get("valid"),
            )
            selected_row.update(_metric_columns(phi_result.get("metrics")))
            selected_row["selection_score"] = phi_result.get("selection_score")
            selected_row["selected_for_stage"] = True
            selected_row["test_model_chi2"] = phi_result.get("test_model_chi2")
            rows.append(selected_row)

    return rows


def _rerun_final_rand_sub(histlist, inpDict, selected_bg_scale1s, selected_bg_scale2s):
    sys.path.append("cuts")
    from rand_sub import rand_sub

    final_inp = dict(inpDict)
    final_inp["bg_stat_scale1_by_setting"] = dict(selected_bg_scale1s)
    final_inp["bg_stat_scale2_by_setting"] = dict(selected_bg_scale2s)

    refreshed_histlist = []
    for base_hist in histlist:
        phi_setting = base_hist["phi_setting"]
        key = get_bg_scale_setting_key(inpDict["EPSSET"], phi_setting)
        scale1 = float(selected_bg_scale1s.get(key, BG_STAT_SCALE1))
        scale2 = float(selected_bg_scale2s.get(key, BG_STAT_SCALE2))
        _log(
            "Final rand_sub rerun for {} {} with BG_STAT_SCALE1={:.3f} BG_STAT_SCALE2={:.3f}".format(
                inpDict["EPSSET"],
                phi_setting,
                scale1,
                scale2,
            )
        )
        refreshed_hist = rand_sub(phi_setting, final_inp, shift_mode="raw", emit_plots=True)
        refreshed_histlist.append(_copy_static_hist_context(base_hist, refreshed_hist))
    return refreshed_histlist


def _load_common_scale_maps_from_low(inpDict, histlist):
    bundle = load_zeroth_iteration_input_bundle(
        inpDict["OUTPATH"],
        inpDict["ParticleType"],
        inpDict["Q2"],
        inpDict["W"],
        active_profile=inpDict.get("bg_active_profile"),
    ) or {}
    low_outfilename = bundle.get("low", {}).get("inp_dict", {}).get("OutFilename")
    if not low_outfilename:
        low_outfilename = str(inpDict.get("OutFilename", "")).replace("highe", "lowe")
    low_summary_path = os.path.join(
        inpDict["OUTPATH"],
        "{}_{}_bg_opt_low.json".format(inpDict["ParticleType"], low_outfilename),
    )
    if not os.path.exists(low_summary_path):
        raise FileNotFoundError(
            "Common-epsilon profile requires the low-epsilon BG summary, but it was not found: {}".format(
                low_summary_path
            )
        )
    with open(low_summary_path, "r") as handle:
        low_summary = json.load(handle)
    selected_bg_scale1s = deepcopy(low_summary.get("selected_bg_scale1s", {}))
    selected_bg_scale2s = deepcopy(low_summary.get("selected_bg_scale2s", {}))
    phi_results = []
    proposal = {
        "t_bins": np.array(histlist[0]["t_bins"], dtype=float),
        "phi_bins": np.array(histlist[0]["phi_bins"], dtype=float),
    }
    for base_hist in histlist:
        phi_setting = base_hist["phi_setting"]
        key = get_bg_scale_setting_key(inpDict["EPSSET"], phi_setting)
        low_key = get_bg_scale_setting_key("low", phi_setting)
        common_scale1 = float(selected_bg_scale1s.get(low_key, selected_bg_scale1s.get(key, BG_STAT_SCALE1)))
        common_scale2 = float(selected_bg_scale2s.get(low_key, selected_bg_scale2s.get(key, BG_STAT_SCALE2)))
        selected_bg_scale1s[key] = common_scale1
        selected_bg_scale2s[key] = common_scale2
        simc_yield_dict, simc_support = _build_simc_reference(base_hist, inpDict, proposal["t_bins"], proposal["phi_bins"])
        sys.path.append("binning")
        from calculate_yield import prepare_bg_opt_data_base_cache

        data_base_cache = prepare_bg_opt_data_base_cache(base_hist, inpDict, proposal["t_bins"], proposal["phi_bins"])
        result = _evaluate_fixed_result(
            base_hist,
            inpDict,
            phi_setting,
            proposal["t_bins"],
            proposal["phi_bins"],
            simc_yield_dict,
            simc_support,
            data_base_cache,
            common_scale1,
            common_scale2,
            "common_epsilon_low_freeze",
        )
        result["common_scale_source"] = "low_epsilon_freeze"
        phi_results.append(result)
    return selected_bg_scale1s, selected_bg_scale2s, phi_results


def optimize_low_epsilon_configuration(histlist, inpDict):
    _log("Starting low-e optimization prepass")
    profile_meta = _get_profile_metadata()
    inpDict["bg_active_profile"] = profile_meta["active_profile"]
    inpDict["bg_resolved_profile"] = profile_meta["resolved_profile"]
    inpDict["bg_use_common_epsilon_scales"] = profile_meta["use_common_epsilon_scales"]
    inpDict["bg_common_epsilon_scale_behavior"] = profile_meta["common_epsilon_scale_behavior"]
    requested_t_bins = int(inpDict["NumtBins"])
    requested_phi_bins = int(inpDict["NumPhiBins"])
    candidate_results = []

    for candidate_t_bins, candidate_phi_bins in build_bin_count_candidates(
        requested_t_bins, requested_phi_bins
    ):
        _log("Trying shared bin candidate NumtBins={} NumPhiBins={}".format(candidate_t_bins, candidate_phi_bins))
        try:
            proposal = propose_bins(
                histlist,
                inpDict,
                num_t_bins=candidate_t_bins,
                num_phi_bins=candidate_phi_bins,
                quiet=True,
            )
        except Exception as exc:
            candidate_results.append(
                {
                    "valid": False,
                    "requested_num_t_bins": candidate_t_bins,
                    "requested_num_phi_bins": candidate_phi_bins,
                    "error": str(exc),
                }
            )
            continue

        phi_results = [
            _optimize_phi_scale(base_hist, inpDict, proposal["t_bins"], proposal["phi_bins"])
            for base_hist in histlist
        ]
        aggregate_metrics = _aggregate_bin_candidate_result(phi_results)
        _log(
            "Candidate NumtBins={} NumPhiBins={} -> fail={} mean_dev={:.4f} rms={:.4f} kin={:.4f}".format(
                proposal["actual_num_t_bins"],
                proposal["actual_num_phi_bins"],
                int(aggregate_metrics["ratio_fail_count"]),
                float(aggregate_metrics["ratio_mean_dev"]),
                float(aggregate_metrics["ratio_rms"]),
                float(aggregate_metrics["kinematic_score"]),
            )
        )
        candidate_results.append(
            {
                "valid": math.isfinite(float(aggregate_metrics["ratio_mean_dev"])),
                "requested_num_t_bins": candidate_t_bins,
                "requested_num_phi_bins": candidate_phi_bins,
                "proposal": proposal,
                "phi_results": phi_results,
                "metrics": aggregate_metrics,
            }
        )

    valid_candidates = [result for result in candidate_results if result.get("valid")]
    if not valid_candidates:
        fallback_proposal = propose_bins(histlist, inpDict, quiet=True)
        apply_bin_proposal(inpDict, fallback_proposal)
        write_bin_interval_files(inpDict, fallback_proposal["t_bins"], fallback_proposal["phi_bins"])
        fallback_map = {
            get_bg_scale_setting_key(inpDict["EPSSET"], hist["phi_setting"]): float(BG_STAT_SCALE2)
            for hist in histlist
        }
        fallback_map1 = {
            get_bg_scale_setting_key(inpDict["EPSSET"], hist["phi_setting"]): float(BG_STAT_SCALE1)
            for hist in histlist
        }
        inpDict["bg_stat_scale1_by_setting"] = fallback_map1
        inpDict["bg_stat_scale2_by_setting"] = fallback_map
        _log("Falling back to configured binning/scales for low epsilon")
        summary = {
            "mode": "low",
            "selection_mode": _get_selection_mode(),
            "active_profile": profile_meta["active_profile"],
            "resolved_profile": profile_meta["resolved_profile"],
            "forced_bg_scale1": profile_meta["forced_bg_scale1"],
            "forced_bg_scale2": profile_meta["forced_bg_scale2"],
            "use_common_epsilon_scales": profile_meta["use_common_epsilon_scales"],
            "common_epsilon_scale_behavior": profile_meta["common_epsilon_scale_behavior"],
            "fallback": True,
            "proposal": fallback_proposal,
            "selected_bg_scale1s": fallback_map1,
            "candidate_results": candidate_results,
            "selected_bg_scale2s": fallback_map,
            "selected_bg_scales": fallback_map,
        }
        inpDict["bg_optimization_summary"] = _serialize_result(summary)
        inpDict["bg_optimization_report"] = _build_summary_report(
            inpDict,
            "low fallback",
            fallback_proposal,
            [
                {
                    "phi_setting": hist["phi_setting"],
                    "bg_scale1": BG_STAT_SCALE1,
                    "bg_scale2": BG_STAT_SCALE2,
                    "fallback": True,
                }
                for hist in histlist
            ],
        )
        return histlist, summary

    _annotate_selection_scores(candidate_results)
    valid_candidates.sort(key=_candidate_order_key)
    best_candidate = valid_candidates[0]
    _log(
        "Selected shared low-e bins NumtBins={} NumPhiBins={} score={}".format(
            best_candidate["proposal"]["actual_num_t_bins"],
            best_candidate["proposal"]["actual_num_phi_bins"],
            _format_metric_value(best_candidate.get("selection_score")),
        )
    )
    apply_bin_proposal(inpDict, best_candidate["proposal"])
    write_bin_interval_files(inpDict, best_candidate["proposal"]["t_bins"], best_candidate["proposal"]["phi_bins"])

    selected_histlist = []
    selected_bg_scale1s = {}
    selected_bg_scale2s = {}
    for result in best_candidate["phi_results"]:
        key = get_bg_scale_setting_key(inpDict["EPSSET"], result["phi_setting"])
        selected_bg_scale1s[key] = float(result["bg_scale1"])
        selected_bg_scale2s[key] = float(result["bg_scale2"])
        _log(
            "Selected low-e scales for {}: BG_STAT_SCALE1={:.3f} BG_STAT_SCALE2={:.3f}".format(
                result["phi_setting"],
                float(result["bg_scale1"]),
                float(result["bg_scale2"]),
            )
        )
    selected_histlist = _rerun_final_rand_sub(histlist, inpDict, selected_bg_scale1s, selected_bg_scale2s)

    inpDict["bg_stat_scale1_by_setting"] = selected_bg_scale1s
    inpDict["bg_stat_scale2_by_setting"] = selected_bg_scale2s
    summary = {
        "mode": "low",
        "selection_mode": _get_selection_mode(),
        "active_profile": profile_meta["active_profile"],
        "resolved_profile": profile_meta["resolved_profile"],
        "forced_bg_scale1": profile_meta["forced_bg_scale1"],
        "forced_bg_scale2": profile_meta["forced_bg_scale2"],
        "use_common_epsilon_scales": profile_meta["use_common_epsilon_scales"],
        "common_epsilon_scale_behavior": profile_meta["common_epsilon_scale_behavior"],
        "fallback": False,
        "proposal": best_candidate["proposal"],
        "selected_bg_scale1s": selected_bg_scale1s,
        "selected_bg_scale2s": selected_bg_scale2s,
        "selected_bg_scales": selected_bg_scale2s,
        "metrics": best_candidate["metrics"],
        "candidate_results": candidate_results,
        "selection_score": best_candidate.get("selection_score"),
    }
    inpDict["bg_optimization_summary"] = _serialize_result(summary)
    inpDict["bg_optimization_report"] = _build_summary_report(
        inpDict,
        "low",
        best_candidate["proposal"],
        best_candidate["phi_results"],
    )
    return selected_histlist, summary


def optimize_high_epsilon_configuration(histlist, inpDict):
    _log("Starting high-e optimization prepass")
    profile_meta = _get_profile_metadata()
    inpDict["bg_active_profile"] = profile_meta["active_profile"]
    inpDict["bg_resolved_profile"] = profile_meta["resolved_profile"]
    inpDict["bg_use_common_epsilon_scales"] = profile_meta["use_common_epsilon_scales"]
    inpDict["bg_common_epsilon_scale_behavior"] = profile_meta["common_epsilon_scale_behavior"]
    proposal = {
        "t_bins": np.array(histlist[0]["t_bins"], dtype=float),
        "phi_bins": np.array(histlist[0]["phi_bins"], dtype=float),
        "requested_num_t_bins": len(histlist[0]["t_bins"]) - 1,
        "requested_num_phi_bins": len(histlist[0]["phi_bins"]) - 1,
        "actual_num_t_bins": len(histlist[0]["t_bins"]) - 1,
        "actual_num_phi_bins": len(histlist[0]["phi_bins"]) - 1,
    }
    _log(
        "Using low-e shared bins for high-e: NumtBins={} NumPhiBins={}".format(
            len(proposal["t_bins"]) - 1,
            len(proposal["phi_bins"]) - 1,
        )
    )

    if use_common_epsilon_scales():
        _log(
            "Active profile requests '{}' behavior; reusing the frozen low-e scale map on the high-e pass".format(
                profile_meta["common_epsilon_scale_behavior"]
            )
        )
        selected_bg_scale1s, selected_bg_scale2s, phi_results = _load_common_scale_maps_from_low(inpDict, histlist)
    else:
        phi_results = [
            _optimize_phi_scale(base_hist, inpDict, proposal["t_bins"], proposal["phi_bins"])
            for base_hist in histlist
        ]

        selected_bg_scale1s = {}
        selected_bg_scale2s = {}
        for result in phi_results:
            key = get_bg_scale_setting_key(inpDict["EPSSET"], result["phi_setting"])
            selected_bg_scale1s[key] = float(result.get("bg_scale1", resolve_bg_stat_scale1(inpDict, result["phi_setting"])))
            selected_bg_scale2s[key] = float(result.get("bg_scale2", resolve_bg_stat_scale2(inpDict, result["phi_setting"])))
            _log(
                "Selected high-e scales for {}: BG_STAT_SCALE1={:.3f} BG_STAT_SCALE2={:.3f}".format(
                    result["phi_setting"],
                    float(selected_bg_scale1s[key]),
                    float(selected_bg_scale2s[key]),
                )
            )
    selected_histlist = _rerun_final_rand_sub(histlist, inpDict, selected_bg_scale1s, selected_bg_scale2s)

    inpDict["bg_stat_scale1_by_setting"] = selected_bg_scale1s
    inpDict["bg_stat_scale2_by_setting"] = selected_bg_scale2s
    summary = {
        "mode": "high",
        "selection_mode": _get_selection_mode(),
        "active_profile": profile_meta["active_profile"],
        "resolved_profile": profile_meta["resolved_profile"],
        "forced_bg_scale1": profile_meta["forced_bg_scale1"],
        "forced_bg_scale2": profile_meta["forced_bg_scale2"],
        "use_common_epsilon_scales": profile_meta["use_common_epsilon_scales"],
        "common_epsilon_scale_behavior": profile_meta["common_epsilon_scale_behavior"],
        "fallback": False,
        "proposal": proposal,
        "selected_bg_scale1s": selected_bg_scale1s,
        "selected_bg_scale2s": selected_bg_scale2s,
        "selected_bg_scales": selected_bg_scale2s,
        "metrics": _aggregate_bin_candidate_result(phi_results),
        "candidate_results": phi_results,
        "selection_score": None,
    }
    if use_common_epsilon_scales():
        summary["common_scale_source"] = "low_epsilon_freeze"
    inpDict["bg_optimization_summary"] = _serialize_result(summary)
    inpDict["bg_optimization_report"] = _build_summary_report(inpDict, "high", proposal, phi_results)
    return selected_histlist, summary


def write_optimization_summary(summary, summary_path):
    serializable = _serialize_result(summary)
    with open(summary_path, "w") as handle:
        json.dump(serializable, handle, indent=2)
    return summary_path


def write_optimization_csv(summary, inpDict, csv_path):
    rows = build_optimization_csv_rows(summary, inpDict)
    fieldnames = [
        "mode",
        "active_profile",
        "selection_mode",
        "epsset",
        "particle",
        "q2",
        "w",
        "outfilename",
        "row_kind",
        "stage",
        "scale_axis",
        "phi_setting",
        "bg_stat_scale1",
        "bg_stat_scale2",
        "valid",
        "fallback",
        "forced_bg_scale1",
        "forced_bg_scale2",
        "use_common_epsilon_scales",
        "common_epsilon_scale_behavior",
        "resolved_profile_json",
        "requested_num_t_bins",
        "requested_num_phi_bins",
        "actual_num_t_bins",
        "actual_num_phi_bins",
        "selected_bin_candidate",
        "selected_for_stage",
        "ratio_fail_count",
        "ratio_mean_dev",
        "ratio_rms",
        "kinematic_score",
        "valid_ratio_bins",
        "kinematic_points",
        "selection_score",
        "test_model_chi2",
        "error",
    ]
    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})
    return csv_path
