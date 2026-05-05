#! /usr/bin/python

import json
import math
import os
import subprocess
import sys

import numpy as np

from background_config import (
    BG_STAT_SCALE2,
    BG_STAT_SCALE2_FINALIST_COUNT,
    KINEMATIC_SCORE_VARS,
    RATIO_SIGMA_THRESHOLD,
    build_bin_count_candidates,
    get_bg_scale_coarse_candidates,
    get_bg_scale_refined_candidates,
    get_bg_scale_setting_key,
    resolve_bg_stat_scale2,
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


def _candidate_sort_key(result):
    return (
        int(result["metrics"]["ratio_fail_count"]),
        float(result["metrics"]["ratio_mean_dev"]),
        float(result["metrics"]["ratio_rms"]),
        float(result["metrics"]["kinematic_score"]),
    )


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
            "{} result {} {} BG_STAT_SCALE2={:.3f} -> invalid (valid_ratio_bins={} kin_points={})".format(
                stage_name,
                epsset,
                phi_setting,
                float(result.get("bg_scale", BG_STAT_SCALE2)),
                int(result.get("metrics", {}).get("valid_ratio_bins", 0)),
                int(result.get("metrics", {}).get("kinematic_points", 0)),
            )
        )
        return

    metrics = result["metrics"]
    _log(
        "{} result {} {} BG_STAT_SCALE2={:.3f} -> fail={} mean_dev={} rms={} kin={} bins={} kin_pts={}".format(
            stage_name,
            epsset,
            phi_setting,
            float(result["bg_scale"]),
            int(metrics["ratio_fail_count"]),
            _format_metric_value(metrics["ratio_mean_dev"]),
            _format_metric_value(metrics["ratio_rms"]),
            _format_metric_value(metrics["kinematic_score"]),
            int(metrics["valid_ratio_bins"]),
            int(metrics["kinematic_points"]),
        )
    )


def _evaluate_scale_grid(base_hist, inpDict, phi_setting, t_bins, phi_bins, simc_yield_dict, simc_support, stage_name, bg_scales):
    results = []
    total = len(bg_scales)
    for idx, bg_scale in enumerate(bg_scales, start=1):
        _log(
            "{} scan {} {}: testing BG_STAT_SCALE2={:.3f} ({}/{})".format(
                stage_name,
                inpDict["EPSSET"],
                phi_setting,
                float(bg_scale),
                idx,
                total,
            )
        )
        result = _evaluate_phi_candidate(
            base_hist,
            inpDict,
            phi_setting,
            bg_scale,
            t_bins,
            phi_bins,
            simc_yield_dict,
            simc_support,
        )
        _log_scale_result(stage_name, inpDict["EPSSET"], phi_setting, result)
        results.append(result)
    return results


def _evaluate_phi_candidate(base_hist, inpDict, phi_setting, bg_scale, t_bins, phi_bins, simc_yield_dict, simc_support):
    candidate_inp = dict(inpDict)
    candidate_inp["NumtBins"] = len(t_bins) - 1
    candidate_inp["NumPhiBins"] = len(phi_bins) - 1
    candidate_inp["yield_emit_plots"] = False
    candidate_inp["suppress_bg_opt_warnings"] = True
    candidate_inp["bg_stat_scale2"] = float(bg_scale)
    candidate_inp["bg_stat_scale2_by_setting"] = {
        get_bg_scale_setting_key(candidate_inp["EPSSET"], phi_setting): float(bg_scale)
    }

    sys.path.append("binning")
    from calculate_yield import find_yield_data
    from calculate_ratio import find_ratio

    candidate_hist = _prepare_candidate_hist(base_hist, t_bins, phi_bins)
    candidate_hist["_xsect_support_simc"] = simc_support

    yield_dict = {}
    yield_dict.update(find_yield_data([candidate_hist], candidate_inp))
    yield_dict.update(simc_yield_dict)

    ratio_dict = {}
    ratio_dict.update(find_ratio([candidate_hist], candidate_inp, yield_dict))

    metrics = _build_candidate_metrics(phi_setting, ratio_dict, candidate_hist)
    return {
        "phi_setting": phi_setting,
        "bg_scale": float(bg_scale),
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
    valid_results = [result for result in results if result.get("valid")]
    if not valid_results:
        return {
            "phi_setting": results[0]["phi_setting"] if results else "unknown",
            "bg_scale": float(BG_STAT_SCALE2),
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

    valid_results.sort(key=_candidate_sort_key)
    finalists = valid_results[: max(1, BG_STAT_SCALE2_FINALIST_COUNT)]
    if len(finalists) > 1:
        top_key = _candidate_sort_key(finalists[0])
        tied = [result for result in finalists if _candidate_sort_key(result) == top_key]
        if len(tied) > 1:
            _log(
                "Running test_model.sh tie-break for {} {} among BG_STAT_SCALE2={}".format(
                    inpDict["EPSSET"],
                    results[0]["phi_setting"],
                    ", ".join("{:.3f}".format(float(result["bg_scale"])) for result in tied),
                )
            )
            tiebreak_chi2 = _run_test_model_tiebreak(inpDict)
            if tiebreak_chi2 is not None:
                tied[0]["test_model_chi2"] = tiebreak_chi2
                finalists = tied

    best_result = valid_results[0]
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
    simc_yield_dict, simc_support = _build_simc_reference(base_hist, inpDict, t_bins, phi_bins)
    coarse_candidates = list(get_bg_scale_coarse_candidates())
    _log(
        "Coarse BG_STAT_SCALE2 grid for {} {}: {}".format(
            inpDict["EPSSET"],
            phi_setting,
            ", ".join("{:.3f}".format(float(scale)) for scale in coarse_candidates),
        )
    )
    coarse_results = _evaluate_scale_grid(
        base_hist,
        inpDict,
        phi_setting,
        t_bins,
        phi_bins,
        simc_yield_dict,
        simc_support,
        "Coarse",
        coarse_candidates,
    )

    best_coarse = _finalize_phi_results(coarse_results, inpDict)
    if not best_coarse.get("valid"):
        _log("No valid coarse BG scale found for {} {}".format(inpDict["EPSSET"], phi_setting))
        return best_coarse

    _log(
        "Best coarse {} {} BG_STAT_SCALE2={:.3f} fail={} mean_dev={:.4f} rms={:.4f}".format(
            inpDict["EPSSET"],
            phi_setting,
            float(best_coarse["bg_scale"]),
            int(best_coarse["metrics"]["ratio_fail_count"]),
            float(best_coarse["metrics"]["ratio_mean_dev"]),
            float(best_coarse["metrics"]["ratio_rms"]),
        )
    )

    refined_candidates = list(get_bg_scale_refined_candidates(best_coarse["bg_scale"]))
    _log(
        "Refined BG_STAT_SCALE2 grid for {} {}: {}".format(
            inpDict["EPSSET"],
            phi_setting,
            ", ".join("{:.3f}".format(float(scale)) for scale in refined_candidates),
        )
    )
    refined_results = _evaluate_scale_grid(
        base_hist,
        inpDict,
        phi_setting,
        t_bins,
        phi_bins,
        simc_yield_dict,
        simc_support,
        "Refined",
        refined_candidates,
    )
    best_refined = _finalize_phi_results(refined_results, inpDict)
    best_refined["coarse_results"] = coarse_results
    if best_refined.get("valid"):
        _log(
            "Selected {} {} BG_STAT_SCALE2={:.3f} fail={} mean_dev={:.4f} rms={:.4f} kin={:.4f}".format(
                inpDict["EPSSET"],
                phi_setting,
                float(best_refined["bg_scale"]),
                int(best_refined["metrics"]["ratio_fail_count"]),
                float(best_refined["metrics"]["ratio_mean_dev"]),
                float(best_refined["metrics"]["ratio_rms"]),
                float(best_refined["metrics"]["kinematic_score"]),
            )
        )
    return best_refined


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


def _serialize_result(result):
    if isinstance(result, dict):
        serializable = {}
        for key, val in result.items():
            if key == "hist":
                continue
            serializable[key] = _serialize_result(val)
        return serializable
    if isinstance(result, list):
        return [_serialize_result(entry) for entry in result]
    if isinstance(result, tuple):
        return [_serialize_result(entry) for entry in result]
    if isinstance(result, np.ndarray):
        return result.tolist()
    if isinstance(result, np.generic):
        return result.item()
    return result


def _build_summary_report(inpDict, mode, proposal, phi_results):
    lines = [
        "BG/Tuning Summary ({})".format(mode),
        "Shared bins: NumtBins={} NumPhiBins={}".format(
            len(proposal["t_bins"]) - 1,
            len(proposal["phi_bins"]) - 1,
        ),
    ]
    for result in phi_results:
        scale = result.get("bg_scale", BG_STAT_SCALE2)
        status = "fallback" if result.get("fallback") else "selected"
        lines.append(
            "{} {}: BG_STAT_SCALE2={:.3f}".format(
                result.get("phi_setting", "unknown"),
                status,
                float(scale),
            )
        )
    return "\n".join(lines)


def _rerun_final_rand_sub(histlist, inpDict, selected_bg_scales):
    sys.path.append("cuts")
    from rand_sub import rand_sub

    final_inp = dict(inpDict)
    final_inp["bg_stat_scale2_by_setting"] = dict(selected_bg_scales)

    refreshed_histlist = []
    for base_hist in histlist:
        phi_setting = base_hist["phi_setting"]
        scale = float(selected_bg_scales.get(get_bg_scale_setting_key(inpDict["EPSSET"], phi_setting), BG_STAT_SCALE2))
        _log("Final rand_sub rerun for {} {} with BG_STAT_SCALE2={:.3f}".format(inpDict["EPSSET"], phi_setting, scale))
        refreshed_hist = rand_sub(phi_setting, final_inp, shift_mode="raw", emit_plots=True)
        refreshed_histlist.append(_copy_static_hist_context(base_hist, refreshed_hist))
    return refreshed_histlist


def optimize_low_epsilon_configuration(histlist, inpDict):
    _log("Starting low-e optimization prepass")
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
        inpDict["bg_stat_scale2_by_setting"] = fallback_map
        _log("Falling back to configured binning/scales for low epsilon")
        summary = {
            "mode": "low",
            "fallback": True,
            "proposal": fallback_proposal,
            "candidate_results": candidate_results,
            "selected_bg_scales": fallback_map,
        }
        inpDict["bg_optimization_summary"] = _serialize_result(summary)
        inpDict["bg_optimization_report"] = _build_summary_report(
            inpDict,
            "low fallback",
            fallback_proposal,
            [{"phi_setting": hist["phi_setting"], "bg_scale": BG_STAT_SCALE2, "fallback": True} for hist in histlist],
        )
        return histlist, summary

    valid_candidates.sort(
        key=lambda result: (
            int(result["metrics"]["ratio_fail_count"]),
            float(result["metrics"]["ratio_mean_dev"]),
            float(result["metrics"]["ratio_rms"]),
            float(result["metrics"]["kinematic_score"]),
        )
    )
    best_candidate = valid_candidates[0]
    _log(
        "Selected shared low-e bins NumtBins={} NumPhiBins={}".format(
            best_candidate["proposal"]["actual_num_t_bins"],
            best_candidate["proposal"]["actual_num_phi_bins"],
        )
    )
    apply_bin_proposal(inpDict, best_candidate["proposal"])
    write_bin_interval_files(inpDict, best_candidate["proposal"]["t_bins"], best_candidate["proposal"]["phi_bins"])

    selected_histlist = []
    selected_bg_scales = {}
    for result in best_candidate["phi_results"]:
        key = get_bg_scale_setting_key(inpDict["EPSSET"], result["phi_setting"])
        selected_bg_scales[key] = float(result["bg_scale"])
        _log(
            "Selected low-e scale for {}: BG_STAT_SCALE2={:.3f}".format(
                result["phi_setting"],
                float(result["bg_scale"]),
            )
        )
    selected_histlist = _rerun_final_rand_sub(histlist, inpDict, selected_bg_scales)

    inpDict["bg_stat_scale2_by_setting"] = selected_bg_scales
    summary = {
        "mode": "low",
        "fallback": False,
        "proposal": best_candidate["proposal"],
        "selected_bg_scales": selected_bg_scales,
        "metrics": best_candidate["metrics"],
        "candidate_results": candidate_results,
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
    proposal = {
        "t_bins": np.array(histlist[0]["t_bins"], dtype=float),
        "phi_bins": np.array(histlist[0]["phi_bins"], dtype=float),
    }

    phi_results = [
        _optimize_phi_scale(base_hist, inpDict, proposal["t_bins"], proposal["phi_bins"])
        for base_hist in histlist
    ]

    selected_bg_scales = {}
    for result in phi_results:
        key = get_bg_scale_setting_key(inpDict["EPSSET"], result["phi_setting"])
        selected_bg_scales[key] = float(result.get("bg_scale", resolve_bg_stat_scale2(inpDict, result["phi_setting"])))
        _log(
            "Selected high-e scale for {}: BG_STAT_SCALE2={:.3f}".format(
                result["phi_setting"],
                float(selected_bg_scales[key]),
            )
        )
    selected_histlist = _rerun_final_rand_sub(histlist, inpDict, selected_bg_scales)

    inpDict["bg_stat_scale2_by_setting"] = selected_bg_scales
    summary = {
        "mode": "high",
        "fallback": False,
        "proposal": proposal,
        "selected_bg_scales": selected_bg_scales,
        "metrics": _aggregate_bin_candidate_result(phi_results),
        "candidate_results": phi_results,
    }
    inpDict["bg_optimization_summary"] = _serialize_result(summary)
    inpDict["bg_optimization_report"] = _build_summary_report(inpDict, "high", proposal, phi_results)
    return selected_histlist, summary


def write_optimization_summary(summary, summary_path):
    serializable = _serialize_result(summary)
    with open(summary_path, "w") as handle:
        json.dump(serializable, handle, indent=2)
    return summary_path
