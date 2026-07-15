#! /usr/bin/python

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Any


DEFAULT_ROOT_HISTOGRAM_NAMES = (
    "H_proton_cleaning_global_pid",
)


def _load_json(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _is_finite_number(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except Exception:
        return False


def _numbers_close(left: Any, right: Any, abs_tol: float, rel_tol: float) -> bool:
    if not (_is_finite_number(left) and _is_finite_number(right)):
        return left == right
    return math.isclose(float(left), float(right), abs_tol=abs_tol, rel_tol=rel_tol)


def _append_issue(issues: list[dict[str, Any]], scope: str, message: str, details: dict[str, Any] | None = None) -> None:
    entry = {
        "scope": str(scope),
        "message": str(message),
    }
    if details:
        entry["details"] = details
    issues.append(entry)


def _compare_value(
    issues: list[dict[str, Any]],
    scope: str,
    label: str,
    left: Any,
    right: Any,
    abs_tol: float,
    rel_tol: float,
) -> None:
    if _is_finite_number(left) and _is_finite_number(right):
        if not _numbers_close(left, right, abs_tol, rel_tol):
            _append_issue(
                issues,
                scope,
                "{} mismatch".format(label),
                {"python": float(left), "macro": float(right)},
            )
        return
    if left != right:
        _append_issue(
            issues,
            scope,
            "{} mismatch".format(label),
            {"python": left, "macro": right},
        )


def _compare_shape_collection(
    issues: list[dict[str, Any]],
    scope: str,
    python_shapes: list[dict[str, Any]],
    macro_shapes: list[dict[str, Any]],
    keys: tuple[str, ...],
    abs_tol: float,
    rel_tol: float,
) -> None:
    if len(python_shapes) != len(macro_shapes):
        _append_issue(
            issues,
            scope,
            "shape-count mismatch",
            {"python_count": len(python_shapes), "macro_count": len(macro_shapes)},
        )
        return
    for index, (left, right) in enumerate(zip(python_shapes, macro_shapes)):
        sub_scope = "{}[{}]".format(scope, index)
        for key in keys:
            _compare_value(
                issues,
                sub_scope,
                key,
                (left or {}).get(key),
                (right or {}).get(key),
                abs_tol,
                rel_tol,
            )


def _compare_nested_slice_fits(
    issues: list[dict[str, Any]],
    python_payload: dict[str, Any],
    macro_payload: dict[str, Any],
    abs_tol: float,
    rel_tol: float,
) -> None:
    python_delta = python_payload.get("delta_slice_fits") or []
    macro_delta = macro_payload.get("delta_slice_fits") or []
    if len(python_delta) != len(macro_delta):
        _append_issue(
            issues,
            "delta_slice_fits",
            "delta-bin count mismatch",
            {"python_count": len(python_delta), "macro_count": len(macro_delta)},
        )
        return
    keys = (
        "fit_attempted",
        "valid",
        "fit_status_code",
        "kaon_amplitude",
        "kaon_amplitude_error",
        "proton_amplitude",
        "proton_amplitude_error",
        "other_amplitude",
        "other_amplitude_error",
        "kaon_yield",
        "proton_yield",
        "other_yield",
        "data_yield",
        "model_yield",
        "model_data_ratio",
        "chi2_data",
        "chi2_ndf",
        "goodness_ndf",
        "rejection_reason",
    )
    for delta_index, (left_slices, right_slices) in enumerate(zip(python_delta, macro_delta)):
        if len(left_slices or []) != len(right_slices or []):
            _append_issue(
                issues,
                "delta_slice_fits[{}]".format(delta_index),
                "aero-slice count mismatch",
                {"python_count": len(left_slices or []), "macro_count": len(right_slices or [])},
            )
            continue
        _compare_shape_collection(
            issues,
            "delta_slice_fits[{}]".format(delta_index),
            list(left_slices or []),
            list(right_slices or []),
            keys,
            abs_tol,
            rel_tol,
        )


def _default_histogram_names(payload: dict[str, Any]) -> list[str]:
    histogram_names = list(DEFAULT_ROOT_HISTOGRAM_NAMES)
    aero_edges = payload.get("aero_edges") or (0.0, 3.0, 6.0, 10.0, 15.0, 25.0)
    delta_edges = payload.get("delta_edges") or tuple(float(i) for i in range(11))
    for aero_index in range(max(len(aero_edges) - 1, 0)):
        histogram_names.append("H_proton_cleaning_global_time_slice_{}".format(aero_index))
    for delta_index in range(max(len(delta_edges) - 1, 0)):
        histogram_names.append("H_proton_cleaning_pid_delta_{}".format(delta_index))
        for aero_index in range(max(len(aero_edges) - 1, 0)):
            histogram_names.append(
                "H_proton_cleaning_time_delta_{}_aero_{}".format(delta_index, aero_index)
            )
    return histogram_names


def _compare_root_histograms(
    issues: list[dict[str, Any]],
    python_root: str,
    macro_root: str,
    histogram_names: list[str],
    abs_tol: float,
    rel_tol: float,
) -> None:
    import ROOT  # noqa: PLC0415

    left_file = ROOT.TFile.Open(str(python_root), "READ")
    right_file = ROOT.TFile.Open(str(macro_root), "READ")
    if not left_file or left_file.IsZombie():
        _append_issue(issues, "root", "failed to open python ROOT file", {"path": python_root})
        return
    if not right_file or right_file.IsZombie():
        _append_issue(issues, "root", "failed to open macro ROOT file", {"path": macro_root})
        return
    try:
        for hist_name in histogram_names:
            left_hist = left_file.Get(hist_name)
            right_hist = right_file.Get(hist_name)
            if left_hist is None or right_hist is None:
                _append_issue(
                    issues,
                    "root/{}".format(hist_name),
                    "missing histogram",
                    {
                        "python_exists": bool(left_hist is not None),
                        "macro_exists": bool(right_hist is not None),
                    },
                )
                continue
            if str(left_hist.ClassName()) != str(right_hist.ClassName()):
                _append_issue(
                    issues,
                    "root/{}".format(hist_name),
                    "histogram class mismatch",
                    {"python": str(left_hist.ClassName()), "macro": str(right_hist.ClassName())},
                )
            left_nx = int(left_hist.GetNbinsX())
            right_nx = int(right_hist.GetNbinsX())
            left_ny = int(left_hist.GetNbinsY()) if hasattr(left_hist, "GetNbinsY") else 0
            right_ny = int(right_hist.GetNbinsY()) if hasattr(right_hist, "GetNbinsY") else 0
            if left_nx != right_nx or left_ny != right_ny:
                _append_issue(
                    issues,
                    "root/{}".format(hist_name),
                    "histogram binning mismatch",
                    {
                        "python_nx": left_nx,
                        "macro_nx": right_nx,
                        "python_ny": left_ny,
                        "macro_ny": right_ny,
                    },
                )
                continue
            _compare_value(
                issues,
                "root/{}".format(hist_name),
                "entries",
                left_hist.GetEntries(),
                right_hist.GetEntries(),
                abs_tol,
                rel_tol,
            )
            _compare_value(
                issues,
                "root/{}".format(hist_name),
                "effective_entries",
                left_hist.GetEffectiveEntries(),
                right_hist.GetEffectiveEntries(),
                abs_tol,
                rel_tol,
            )
            for x_bin in range(1, left_nx + 1):
                _compare_value(
                    issues,
                    "root/{}".format(hist_name),
                    "x_low_edge({})".format(x_bin),
                    left_hist.GetXaxis().GetBinLowEdge(x_bin),
                    right_hist.GetXaxis().GetBinLowEdge(x_bin),
                    abs_tol,
                    rel_tol,
                )
                _compare_value(
                    issues,
                    "root/{}".format(hist_name),
                    "x_up_edge({})".format(x_bin),
                    left_hist.GetXaxis().GetBinUpEdge(x_bin),
                    right_hist.GetXaxis().GetBinUpEdge(x_bin),
                    abs_tol,
                    rel_tol,
                )
            if left_ny > 1:
                for y_bin in range(1, left_ny + 1):
                    _compare_value(
                        issues,
                        "root/{}".format(hist_name),
                        "y_low_edge({})".format(y_bin),
                        left_hist.GetYaxis().GetBinLowEdge(y_bin),
                        right_hist.GetYaxis().GetBinLowEdge(y_bin),
                        abs_tol,
                        rel_tol,
                    )
                    _compare_value(
                        issues,
                        "root/{}".format(hist_name),
                        "y_up_edge({})".format(y_bin),
                        left_hist.GetYaxis().GetBinUpEdge(y_bin),
                        right_hist.GetYaxis().GetBinUpEdge(y_bin),
                        abs_tol,
                        rel_tol,
                    )
                for x_bin in range(0, left_nx + 2):
                    for y_bin in range(0, left_ny + 2):
                        _compare_value(
                            issues,
                            "root/{}".format(hist_name),
                            "bin({}, {})".format(x_bin, y_bin),
                            left_hist.GetBinContent(x_bin, y_bin),
                            right_hist.GetBinContent(x_bin, y_bin),
                            abs_tol,
                            rel_tol,
                        )
                        _compare_value(
                            issues,
                            "root/{}".format(hist_name),
                            "binerr({}, {})".format(x_bin, y_bin),
                            left_hist.GetBinError(x_bin, y_bin),
                            right_hist.GetBinError(x_bin, y_bin),
                            abs_tol,
                            rel_tol,
                        )
                continue
            for x_bin in range(0, left_nx + 2):
                _compare_value(
                    issues,
                    "root/{}".format(hist_name),
                    "bin({})".format(x_bin),
                    left_hist.GetBinContent(x_bin),
                    right_hist.GetBinContent(x_bin),
                    abs_tol,
                    rel_tol,
                )
                _compare_value(
                    issues,
                    "root/{}".format(hist_name),
                    "binerr({})".format(x_bin),
                    left_hist.GetBinError(x_bin),
                    right_hist.GetBinError(x_bin),
                    abs_tol,
                    rel_tol,
                )
    finally:
        left_file.Close()
        right_file.Close()


def _load_weight_map(path: str) -> dict[str, float]:
    payload = _load_json(path)
    if isinstance(payload, dict):
        return {str(key): float(value) for key, value in payload.items()}
    if isinstance(payload, list):
        result = {}
        for row in payload:
            if not isinstance(row, dict):
                continue
            signature = row.get("signature")
            value = row.get("proton_weight")
            if signature is None or value is None:
                continue
            result[str(signature)] = float(value)
        return result
    raise TypeError("Unsupported weight payload format for '{}'".format(path))


def _compare_weight_maps(
    issues: list[dict[str, Any]],
    python_weights: dict[str, float],
    macro_weights: dict[str, float],
    abs_tol: float,
    rel_tol: float,
) -> None:
    python_keys = set(python_weights.keys())
    macro_keys = set(macro_weights.keys())
    if python_keys != macro_keys:
        _append_issue(
            issues,
            "event_weights",
            "signature set mismatch",
            {
                "python_only": sorted(list(python_keys - macro_keys))[:25],
                "macro_only": sorted(list(macro_keys - python_keys))[:25],
            },
        )
    for signature in sorted(python_keys & macro_keys):
        _compare_value(
            issues,
            "event_weights",
            signature,
            python_weights[signature],
            macro_weights[signature],
            abs_tol,
            rel_tol,
        )


def _summarize_weight_maps(
    python_weights: dict[str, float] | None,
    macro_weights: dict[str, float] | None,
) -> dict[str, Any]:
    if not python_weights or not macro_weights:
        return {
            "python_weight_count": len(python_weights or {}),
            "macro_weight_count": len(macro_weights or {}),
            "shared_weight_count": 0,
            "missing_from_python_count": 0,
            "missing_from_macro_count": 0,
            "max_abs_weight_difference": None,
        }
    python_keys = set(python_weights.keys())
    macro_keys = set(macro_weights.keys())
    shared_keys = sorted(python_keys & macro_keys)
    max_abs_difference = None
    max_abs_difference_key = None
    for key in shared_keys:
        difference = abs(float(python_weights[key]) - float(macro_weights[key]))
        if max_abs_difference is None or difference > max_abs_difference:
            max_abs_difference = difference
            max_abs_difference_key = key
    return {
        "python_weight_count": len(python_weights),
        "macro_weight_count": len(macro_weights),
        "shared_weight_count": len(shared_keys),
        "missing_from_python_count": len(macro_keys - python_keys),
        "missing_from_macro_count": len(python_keys - macro_keys),
        "max_abs_weight_difference": max_abs_difference,
        "max_abs_weight_difference_key": max_abs_difference_key,
    }


def _write_text_report(path: str, report: dict[str, Any]) -> None:
    lines = [
        "Proton-cleaning parity validation",
        "passed={}".format(report["passed"]),
        "issue_count={}".format(len(report["issues"])),
        "legacy_c_parity_issue_count={}".format(
            int((report.get("legacy_c_parity") or {}).get("issue_count", 0))
        ),
        "python_tof_extension_issue_count={}".format(
            int((report.get("python_tof_extension") or {}).get("issue_count", 0))
        ),
    ]
    extension = report.get("python_tof_extension") or {}
    observations = extension.get("observations") or {}
    if observations:
        lines.append("python_tof_extension_observations={}".format(json.dumps(observations, sort_keys=True)))
    for issue in report["issues"]:
        line = "[{}] {}".format(issue.get("scope", "unknown"), issue.get("message", ""))
        details = issue.get("details")
        if details:
            line += " :: {}".format(json.dumps(details, sort_keys=True))
        lines.append(line)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate prompt-only proton-cleaning parity between the exact Python path and check_slow_protons.C outputs."
    )
    parser.add_argument("--python-json", required=True, help="Serialized Python proton-cleaning result JSON.")
    parser.add_argument("--macro-json", required=True, help="Serialized macro proton-cleaning result JSON.")
    parser.add_argument("--python-root", default="", help="Optional Python ROOT diagnostics file.")
    parser.add_argument("--macro-root", default="", help="Optional macro ROOT diagnostics file.")
    parser.add_argument("--python-weights", default="", help="Optional Python event-weight JSON.")
    parser.add_argument("--macro-weights", default="", help="Optional macro event-weight JSON.")
    parser.add_argument("--report-json", required=True, help="Path to write the JSON report.")
    parser.add_argument("--report-text", required=True, help="Path to write the text report.")
    parser.add_argument("--abs-tol", type=float, default=1e-6)
    parser.add_argument("--rel-tol", type=float, default=1e-6)
    args = parser.parse_args()

    legacy_issues: list[dict[str, Any]] = []
    extension_issues: list[dict[str, Any]] = []
    python_payload = _load_json(args.python_json)
    macro_payload = _load_json(args.macro_json)

    _compare_shape_collection(
        legacy_issues,
        "global_shapes",
        list(python_payload.get("global_shapes") or []),
        list(macro_payload.get("global_shapes") or []),
        (
            "fit_attempted",
            "valid",
            "fit_status_code",
            "kaon_amplitude",
            "kaon_amplitude_error",
            "kaon_mean",
            "kaon_sigma",
            "proton_amplitude",
            "proton_amplitude_error",
            "proton_mean",
            "proton_sigma",
            "other_amplitude",
            "separation",
            "kaon_significance",
            "proton_significance",
            "chi2_data",
            "chi2_ndf",
            "goodness_ndf",
            "rejection_reason",
        ),
        args.abs_tol,
        args.rel_tol,
    )

    if args.python_root and args.macro_root:
        histogram_names = _default_histogram_names(python_payload)
        _compare_root_histograms(
            legacy_issues,
            args.python_root,
            args.macro_root,
            histogram_names,
            args.abs_tol,
            args.rel_tol,
        )

    python_weights = None
    macro_weights = None
    if args.python_weights and args.macro_weights:
        python_weights = _load_weight_map(args.python_weights)
        macro_weights = _load_weight_map(args.macro_weights)

    extension_observations = {
        "note": (
            "TOF-corrected support labels, selected timing branch, delta-slice fits, "
            "and event weights are reported here and are not strict legacy C parity checks."
        ),
        "python_accepted": python_payload.get("accepted"),
        "macro_accepted": macro_payload.get("accepted"),
        "python_selected_timing_branch": python_payload.get("selected_timing_branch"),
        "macro_selected_timing_branch": macro_payload.get("selected_timing_branch"),
        "python_support_by_delta": python_payload.get("support_by_delta"),
        "macro_support_by_delta": macro_payload.get("support_by_delta"),
        "weight_map_summary": _summarize_weight_maps(python_weights, macro_weights),
    }

    issues = legacy_issues + extension_issues

    report = {
        "passed": len(legacy_issues) == 0 and len(extension_issues) == 0,
        "python_json": os.path.abspath(args.python_json),
        "macro_json": os.path.abspath(args.macro_json),
        "python_root": os.path.abspath(args.python_root) if args.python_root else "",
        "macro_root": os.path.abspath(args.macro_root) if args.macro_root else "",
        "python_weights": os.path.abspath(args.python_weights) if args.python_weights else "",
        "macro_weights": os.path.abspath(args.macro_weights) if args.macro_weights else "",
        "abs_tol": float(args.abs_tol),
        "rel_tol": float(args.rel_tol),
        "legacy_c_parity": {
            "passed": len(legacy_issues) == 0,
            "issue_count": len(legacy_issues),
            "issues": legacy_issues,
        },
        "python_tof_extension": {
            "passed": len(extension_issues) == 0,
            "issue_count": len(extension_issues),
            "issues": extension_issues,
            "observations": extension_observations,
        },
        "issues": issues,
    }
    with open(args.report_json, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    _write_text_report(args.report_text, report)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
