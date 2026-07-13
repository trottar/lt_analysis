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
    "H_MM_before_proton_cleaning",
    "H_MM_estimated_proton",
    "H_MM_after_proton_cleaning",
    "H_proton_fraction_vs_MM",
    "H_proton_weight_vs_delta",
    "H_proton_weight_vs_delta_aero",
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
            if hasattr(left_hist, "GetNbinsY") and int(left_hist.GetNbinsY()) > 1:
                if (
                    int(left_hist.GetNbinsX()) != int(right_hist.GetNbinsX())
                    or int(left_hist.GetNbinsY()) != int(right_hist.GetNbinsY())
                ):
                    _append_issue(
                        issues,
                        "root/{}".format(hist_name),
                        "histogram binning mismatch",
                    )
                    continue
                for x_bin in range(1, int(left_hist.GetNbinsX()) + 1):
                    for y_bin in range(1, int(left_hist.GetNbinsY()) + 1):
                        _compare_value(
                            issues,
                            "root/{}".format(hist_name),
                            "bin({}, {})".format(x_bin, y_bin),
                            left_hist.GetBinContent(x_bin, y_bin),
                            right_hist.GetBinContent(x_bin, y_bin),
                            abs_tol,
                            rel_tol,
                        )
                continue
            if int(left_hist.GetNbinsX()) != int(right_hist.GetNbinsX()):
                _append_issue(
                    issues,
                    "root/{}".format(hist_name),
                    "histogram binning mismatch",
                )
                continue
            for x_bin in range(1, int(left_hist.GetNbinsX()) + 1):
                _compare_value(
                    issues,
                    "root/{}".format(hist_name),
                    "bin({})".format(x_bin),
                    left_hist.GetBinContent(x_bin),
                    right_hist.GetBinContent(x_bin),
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


def _write_text_report(path: str, report: dict[str, Any]) -> None:
    lines = [
        "Proton-cleaning parity validation",
        "passed={}".format(report["passed"]),
        "issue_count={}".format(len(report["issues"])),
    ]
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

    issues: list[dict[str, Any]] = []
    python_payload = _load_json(args.python_json)
    macro_payload = _load_json(args.macro_json)

    _compare_value(issues, "summary", "accepted", python_payload.get("accepted"), macro_payload.get("accepted"), args.abs_tol, args.rel_tol)
    _compare_value(issues, "summary", "selected_timing_branch", python_payload.get("selected_timing_branch"), macro_payload.get("selected_timing_branch"), args.abs_tol, args.rel_tol)
    _compare_value(issues, "summary", "support_by_delta", python_payload.get("support_by_delta"), macro_payload.get("support_by_delta"), args.abs_tol, args.rel_tol)

    _compare_shape_collection(
        issues,
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
    _compare_nested_slice_fits(issues, python_payload, macro_payload, args.abs_tol, args.rel_tol)

    if args.python_root and args.macro_root:
        histogram_names = _default_histogram_names(python_payload)
        _compare_root_histograms(
            issues,
            args.python_root,
            args.macro_root,
            histogram_names,
            args.abs_tol,
            args.rel_tol,
        )

    if args.python_weights and args.macro_weights:
        _compare_weight_maps(
            issues,
            _load_weight_map(args.python_weights),
            _load_weight_map(args.macro_weights),
            args.abs_tol,
            args.rel_tol,
        )

    report = {
        "passed": len(issues) == 0,
        "python_json": os.path.abspath(args.python_json),
        "macro_json": os.path.abspath(args.macro_json),
        "python_root": os.path.abspath(args.python_root) if args.python_root else "",
        "macro_root": os.path.abspath(args.macro_root) if args.macro_root else "",
        "python_weights": os.path.abspath(args.python_weights) if args.python_weights else "",
        "macro_weights": os.path.abspath(args.macro_weights) if args.macro_weights else "",
        "abs_tol": float(args.abs_tol),
        "rel_tol": float(args.rel_tol),
        "issues": issues,
    }
    with open(args.report_json, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    _write_text_report(args.report_text, report)
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
