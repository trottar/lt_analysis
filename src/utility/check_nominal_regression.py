#! /usr/bin/python

from __future__ import annotations

import json
import math
import os

from background_config import REGRESSION_ABS_TOL, REGRESSION_REL_TOL


def _load_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


def _write_json(path, payload):
    with open(path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _write_text(path, text):
    with open(path, "w") as handle:
        handle.write(text)


def _try_float(value):
    if value in (None, ""):
        return None
    try:
        return float(value)
    except Exception:
        return None


def _values_close(candidate, reference, abs_tol, rel_tol):
    return abs(float(candidate) - float(reference)) <= float(abs_tol) + float(rel_tol) * abs(float(reference))


def _flatten_structure(prefix, value, output):
    if isinstance(value, dict):
        for key in sorted(value.keys()):
            next_prefix = "{}.{}".format(prefix, key) if prefix else str(key)
            _flatten_structure(next_prefix, value[key], output)
        return
    if isinstance(value, list):
        for index, item in enumerate(value):
            next_prefix = "{}[{}]".format(prefix, index) if prefix else "[{}]".format(index)
            _flatten_structure(next_prefix, item, output)
        return
    output[prefix] = value


def _extract_scalar(rows, candidate_keys):
    if not rows:
        return None
    for row in rows:
        for key in candidate_keys:
            if key in row and row[key] not in ("", None):
                numeric = _try_float(row[key])
                return numeric if numeric is not None else row[key]
    return None


def _collect_regression_fields(summary_payload):
    fields = {}
    _flatten_structure("t_bin_edges", summary_payload.get("t_bin_edges", []), fields)
    _flatten_structure("phi_bin_edges", summary_payload.get("phi_bin_edges", []), fields)
    _flatten_structure("fit1_scales", summary_payload.get("fit1_scales", {}), fields)
    _flatten_structure("fit2_scales", summary_payload.get("fit2_scales", {}), fields)

    low_corr = summary_payload.get("low_epsilon_corrections", {})
    high_corr = summary_payload.get("high_epsilon_corrections", {})
    for prefix, payload in (("low", low_corr), ("high", high_corr)):
        fields["{}.fit1_fractional_correction".format(prefix)] = payload.get("fit1_fractional_correction")
        fields["{}.fit2_fractional_correction".format(prefix)] = payload.get("fit2_fractional_correction")
        fields["{}.total_empirical_fraction".format(prefix)] = payload.get("total_empirical_fraction")
        fields["{}.final_lambda_window".format(prefix)] = (
            payload.get("combined_stage_yields", {}) or {}
        ).get("final_lambda_window")

    low_final = _try_float(fields.get("low.final_lambda_window")) or 0.0
    high_final = _try_float(fields.get("high.final_lambda_window")) or 0.0
    fields["combined.final_lambda_window"] = low_final + high_final

    xsect_outputs = summary_payload.get("xsect_outputs", {})
    fields["xsect.unseparated_cross_section"] = _extract_scalar(
        xsect_outputs.get("unseparated_csv") or [],
        ["xsec", "unsep_xsec", "cross_section"],
    )
    fields["xsect.sigma_L"] = _extract_scalar(
        xsect_outputs.get("separated_csv") or [],
        ["sigma_L", "sigL", "L"],
    )
    fields["xsect.sigma_T"] = _extract_scalar(
        xsect_outputs.get("separated_csv") or [],
        ["sigma_T", "sigT", "T"],
    )
    fields["xsect.sigma_LT"] = _extract_scalar(
        xsect_outputs.get("separated_csv") or [],
        ["sigma_LT", "sigLT", "LT"],
    )
    fields["xsect.sigma_TT"] = _extract_scalar(
        xsect_outputs.get("separated_csv") or [],
        ["sigma_TT", "sigTT", "TT"],
    )

    ratio_summary = summary_payload.get("ratio_summary") or summary_payload.get("data_simc_ratio_summary")
    if ratio_summary is not None:
        _flatten_structure("ratio_summary", ratio_summary, fields)
    return fields


def build_nominal_regression_report(reference_summary, candidate_summary, abs_tol=REGRESSION_ABS_TOL, rel_tol=REGRESSION_REL_TOL):
    reference_fields = _collect_regression_fields(reference_summary)
    candidate_fields = _collect_regression_fields(candidate_summary)
    field_names = sorted(set(reference_fields.keys()) | set(candidate_fields.keys()))
    comparisons = []
    passed = True

    for field_name in field_names:
        reference_value = reference_fields.get(field_name)
        candidate_value = candidate_fields.get(field_name)
        if reference_value is None and candidate_value is None:
            comparisons.append(
                {
                    "field": field_name,
                    "status": "skipped",
                    "reason": "both_missing",
                    "reference": reference_value,
                    "candidate": candidate_value,
                }
            )
            continue
        if reference_value is None or candidate_value is None:
            passed = False
            comparisons.append(
                {
                    "field": field_name,
                    "status": "different",
                    "reason": "missing_value",
                    "reference": reference_value,
                    "candidate": candidate_value,
                }
            )
            continue

        reference_numeric = _try_float(reference_value)
        candidate_numeric = _try_float(candidate_value)
        if reference_numeric is not None and candidate_numeric is not None:
            is_close = _values_close(candidate_numeric, reference_numeric, abs_tol, rel_tol)
            if not is_close:
                passed = False
            comparisons.append(
                {
                    "field": field_name,
                    "status": "matched" if is_close else "different",
                    "reference": reference_numeric,
                    "candidate": candidate_numeric,
                    "abs_diff": abs(candidate_numeric - reference_numeric),
                    "rel_diff": _safe_rel_diff(candidate_numeric, reference_numeric),
                }
            )
            continue

        exact_match = candidate_value == reference_value
        if not exact_match:
            passed = False
        comparisons.append(
            {
                "field": field_name,
                "status": "matched" if exact_match else "different",
                "reference": reference_value,
                "candidate": candidate_value,
            }
        )

    return {
        "passed": passed,
        "abs_tol": abs_tol,
        "rel_tol": rel_tol,
        "comparison_count": len(comparisons),
        "difference_count": sum(1 for item in comparisons if item.get("status") == "different"),
        "comparisons": comparisons,
    }


def _safe_rel_diff(candidate, reference):
    reference = float(reference)
    if abs(reference) <= 1.0e-12:
        return None if abs(float(candidate)) <= 1.0e-12 else float("inf")
    return abs(float(candidate) - reference) / abs(reference)


def _markdown_report(report, reference_path, candidate_path):
    lines = [
        "# Nominal Regression Check",
        "",
        "| Field | Value |",
        "| --- | --- |",
        "| Passed | {} |".format(report.get("passed")),
        "| Reference | {} |".format(reference_path),
        "| Candidate | {} |".format(candidate_path),
        "| Abs tol | {} |".format(report.get("abs_tol")),
        "| Rel tol | {} |".format(report.get("rel_tol")),
        "| Differences | {} |".format(report.get("difference_count")),
        "",
    ]
    differing = [item for item in report.get("comparisons", []) if item.get("status") == "different"]
    skipped = [item for item in report.get("comparisons", []) if item.get("status") == "skipped"]
    if differing:
        lines.extend(
            [
                "## Differences",
                "",
                "| Field | Reference | Candidate | Abs diff | Rel diff |",
                "| --- | ---: | ---: | ---: | ---: |",
            ]
        )
        for item in differing:
            lines.append(
                "| {field} | {reference} | {candidate} | {abs_diff} | {rel_diff} |".format(
                    field=item.get("field"),
                    reference=item.get("reference"),
                    candidate=item.get("candidate"),
                    abs_diff=item.get("abs_diff", ""),
                    rel_diff=item.get("rel_diff", ""),
                )
            )
        lines.append("")
    if skipped:
        lines.extend(
            [
                "## Skipped Fields",
                "",
                "| Field | Reason |",
                "| --- | --- |",
            ]
        )
        for item in skipped:
            lines.append("| {} | {} |".format(item.get("field"), item.get("reason")))
        lines.append("")
    return "\n".join(lines)


def write_nominal_regression_report(report, reference_path, candidate_path, output_prefix=None):
    if output_prefix is None:
        candidate_root, _ = os.path.splitext(candidate_path)
        output_prefix = candidate_root + "_nominal_regression_check"
    json_path = output_prefix + ".json"
    md_path = output_prefix + ".md"
    payload = {
        "reference_summary": reference_path,
        "candidate_summary": candidate_path,
        **report,
    }
    _write_json(json_path, payload)
    _write_text(md_path, _markdown_report(report, reference_path, candidate_path))
    return [json_path, md_path]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Compare a candidate nominal final-summary JSON against a trusted reference.")
    parser.add_argument("reference_summary")
    parser.add_argument("candidate_summary")
    parser.add_argument("--output-prefix")
    parser.add_argument("--abs-tol", type=float, default=REGRESSION_ABS_TOL)
    parser.add_argument("--rel-tol", type=float, default=REGRESSION_REL_TOL)
    args = parser.parse_args()

    reference_summary = _load_json(args.reference_summary)
    candidate_summary = _load_json(args.candidate_summary)
    report = build_nominal_regression_report(
        reference_summary,
        candidate_summary,
        abs_tol=args.abs_tol,
        rel_tol=args.rel_tol,
    )
    write_nominal_regression_report(
        report,
        args.reference_summary,
        args.candidate_summary,
        output_prefix=args.output_prefix,
    )
    if not report.get("passed"):
        raise SystemExit(1)
