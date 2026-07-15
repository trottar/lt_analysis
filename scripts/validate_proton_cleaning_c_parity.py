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
    result = {}
    for row in _load_weight_payloads(path):
        signature = row.get("signature")
        value = row.get("proton_weight")
        if signature is None or value is None:
            continue
        result[str(signature)] = float(value)
    if result:
        return result
    raise TypeError("Unsupported weight payload format for '{}'".format(path))


def _load_weight_payloads(path: str) -> list[dict[str, Any]]:
    payload = _load_json(path)
    if isinstance(payload, dict):
        if all(not isinstance(value, dict) for value in payload.values()):
            return [
                {"signature": str(key), "proton_weight": value}
                for key, value in payload.items()
            ]
        rows = []
        for key, value in payload.items():
            if not isinstance(value, dict):
                continue
            row = dict(value)
            row.setdefault("signature", str(key))
            rows.append(row)
        return rows
    if isinstance(payload, list):
        return [dict(row) for row in payload if isinstance(row, dict)]
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


def _float_or_none(value: Any) -> float | None:
    try:
        numeric = float(value)
    except Exception:
        return None
    return numeric if math.isfinite(numeric) else None


def _diagnostics(payload: dict[str, Any]) -> dict[str, Any]:
    return payload.get("diagnostics") or {}


def _modulo_distance(left: float, right: float, period: float) -> float:
    if period <= 0.0 or not math.isfinite(period):
        return abs(left - right)
    return abs(((left - right + 0.5 * period) % period) - 0.5 * period)


def _validate_tof_summary_rows(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
) -> None:
    diag = _diagnostics(payload)
    rows = list(diag.get("tof_summary_by_delta") or payload.get("tof_summary_by_delta") or [])
    if not rows:
        _append_issue(issues, "python_tof_extension/tof_summary", "missing TOF summary rows")
        return
    default_cfg = diag.get("tof_summary_validation") or {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        delta_index = int(row.get("delta_index", -1))
        scope = "python_tof_extension/tof_summary[{}]".format(delta_index)
        cfg = row.get("tof_summary_validation") or default_cfg
        min_events = int(
            (cfg or {}).get(
                "minimum_prompt_events",
                (cfg or {}).get("minimum_prompt_tof_events", 30),
            )
            or 30
        )
        min_valid_events = int((cfg or {}).get("minimum_valid_tof_events", min_events) or min_events)
        min_fraction = float((cfg or {}).get("minimum_valid_tof_fraction", 0.90) or 0.90)
        counter_keys = (
            "prompt_events_seen",
            "prompt_events_with_selected_timing",
            "prompt_events_inside_timing_range",
            "prompt_events_inside_aero_range",
            "prompt_events_inside_timing_and_aero_domain",
            "prompt_events_with_valid_tof",
            "prompt_events_used",
        )
        counters = [int(row.get(key, 0) or 0) for key in counter_keys]
        if any(left < right for left, right in zip(counters, counters[1:])):
            _append_issue(
                issues,
                scope,
                "TOF counter ordering is inconsistent",
                dict(zip(counter_keys, counters)),
            )
        valid_tof_fraction = _float_or_none(row.get("valid_tof_fraction"))
        if valid_tof_fraction is None:
            _append_issue(issues, scope, "missing valid_tof_fraction")
        denominator = int(row.get("prompt_events_inside_timing_and_aero_domain", 0) or 0)
        usable_tof = int(row.get("usable_tof_events", row.get("prompt_events_used", 0)) or 0)
        expected_fraction = float(usable_tof / denominator) if denominator > 0 else 0.0
        if valid_tof_fraction is not None and abs(valid_tof_fraction - expected_fraction) > 1.0e-12:
            _append_issue(
                issues,
                scope,
                "valid_tof_fraction uses wrong denominator",
                {
                    "valid_tof_fraction": valid_tof_fraction,
                    "expected": expected_fraction,
                    "usable_tof_events": usable_tof,
                    "prompt_events_inside_timing_and_aero_domain": denominator,
                },
            )
        if usable_tof > denominator:
            _append_issue(
                issues,
                scope,
                "usable TOF events exceed active timing+aerogel denominator",
                {"usable_tof_events": usable_tof, "denominator": denominator},
            )
        required_fields = (
            "mean_delta_t_pk_ns",
            "mean_P_gtr_p",
            "mean_shms_path_length_cm",
        )
        finite_required = all(
            (_float_or_none(row.get(field_name)) is not None and _float_or_none(row.get(field_name)) > 0.0)
            for field_name in required_fields
        )
        expected_valid = bool(
            denominator >= min_events
            and usable_tof >= min_valid_events
            and valid_tof_fraction is not None
            and valid_tof_fraction >= min_fraction
            and finite_required
        )
        if bool(row.get("valid", False)) != expected_valid:
            _append_issue(
                issues,
                scope,
                "TOF summary valid flag does not match serialized thresholds",
                {
                    "valid": row.get("valid"),
                    "expected_valid": expected_valid,
                    "prompt_events_inside_timing_and_aero_domain": denominator,
                    "minimum_prompt_events": min_events,
                    "usable_tof_events": usable_tof,
                    "minimum_valid_tof_events": min_valid_events,
                    "valid_tof_fraction": valid_tof_fraction,
                    "minimum_valid_tof_fraction": min_fraction,
                },
            )
        if not bool(row.get("valid", False)) and not row.get("rejection_reasons"):
            _append_issue(issues, scope, "invalid TOF summary has no rejection_reasons")


def _validate_low_aero_offset_metadata(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
) -> None:
    diag = _diagnostics(payload)
    cfg = diag.get("low_aero_offset_config") or {}
    reference_npe = _float_or_none(diag.get("aerogel_reference_npe"))
    if reference_npe is None:
        _append_issue(issues, "python_tof_extension/low_aero", "missing aerogel_reference_npe")
    full_range = diag.get("full_aerogel_range")
    if not isinstance(full_range, list) or len(full_range) != 2:
        _append_issue(issues, "python_tof_extension/low_aero", "missing full_aerogel_range")
    else:
        if _float_or_none(full_range[0]) != 0.0 or _float_or_none(full_range[1]) != 25.0:
            _append_issue(issues, "python_tof_extension/low_aero", "wrong full diagnostic aerogel range", {"full_aerogel_range": full_range})
    full_rows = list(diag.get("full_aero_tof_summary_by_delta") or [])
    selected_rows = list(
        diag.get("selected_tof_summary_by_delta")
        or diag.get("tof_summary_by_delta")
        or payload.get("tof_summary_by_delta")
        or []
    )
    offsets = list(payload.get("delta_timing_offset_fits") or diag.get("delta_timing_offset_fits") or [])
    allowed_selected_sources = {"low_aero_0_5_fit", "low_aero_0_6_fit", "stable_global_center_fallback"}
    for index, row in enumerate(selected_rows):
        if not isinstance(row, dict):
            continue
        scope = "python_tof_extension/low_aero[{}]".format(index)
        selected_source = str(row.get("selected_timing_center_source") or "unavailable")
        expected_source = (
            str((offsets[index] or {}).get("selected_timing_center_source") or "unavailable")
            if index < len(offsets) and isinstance(offsets[index], dict)
            else "unavailable"
        )
        if selected_source != expected_source:
            _append_issue(
                issues,
                scope,
                "selected TOF summary does not follow the resolved timing-center source",
                {"summary_source": selected_source, "resolved_source": expected_source},
            )
        expected_role = (
            "offset_fit_input"
            if selected_source in {"low_aero_0_5_fit", "low_aero_0_6_fit"}
            else "stable_fallback_diagnostic"
            if selected_source == "stable_global_center_fallback"
            else "unavailable"
        )
        if str(row.get("tof_summary_role") or "unavailable") != expected_role:
            _append_issue(
                issues,
                scope,
                "selected TOF summary has the wrong diagnostic role",
                {"role": row.get("tof_summary_role"), "expected_role": expected_role},
            )
        mode = str(row.get("offset_fit_aero_mode") or "unavailable")
        aero_min = _float_or_none(row.get("offset_fit_aero_min"))
        aero_max = _float_or_none(row.get("offset_fit_aero_max"))
        if selected_source == "low_aero_0_5_fit":
            if mode != "low_aero_0_5":
                _append_issue(issues, scope, "0-5 source has the wrong TOF summary mode", {"mode": mode})
            if aero_min != 0.0 or aero_max is None or aero_max > 5.0:
                _append_issue(issues, scope, "primary offset mode is not restricted to 0-5 NPE", {"aero_min": aero_min, "aero_max": aero_max})
        if selected_source == "low_aero_0_6_fit":
            if mode != "low_aero_0_6_fallback":
                _append_issue(issues, scope, "0-6 source has the wrong TOF summary mode", {"mode": mode})
            if aero_min != 0.0 or aero_max is None or aero_max > 6.0:
                _append_issue(issues, scope, "fallback offset mode is not restricted to 0-6 NPE", {"aero_min": aero_min, "aero_max": aero_max})
        if index < len(full_rows) and isinstance(full_rows[index], dict):
            full_max = _float_or_none(full_rows[index].get("offset_fit_aero_max"))
            if full_max is not None and full_max < 25.0:
                _append_issue(issues, scope, "full diagnostic row does not cover 0-25 NPE", {"full_aero_max": full_max})
    for index, row in enumerate(offsets):
        if not isinstance(row, dict):
            continue
        source = str(row.get("selected_timing_center_source") or row.get("timing_center_source") or "unavailable")
        if source not in allowed_selected_sources:
            _append_issue(
                issues,
                "python_tof_extension/low_aero_offset_fit[{}]".format(index),
                "offset fit selected unsupported timing-center source",
                {"selected_timing_center_source": source},
            )
            continue
        primary_attempt = row.get("primary_offset_attempt") or {}
        fallback_attempt = row.get("fallback_offset_attempt") or {}
        primary_valid = bool(row.get("primary_offset_valid", primary_attempt.get("valid", False)))
        fallback_valid = bool(row.get("fallback_offset_valid", fallback_attempt.get("valid", False)))
        if source == "low_aero_0_5_fit" and not primary_valid:
            _append_issue(issues, "python_tof_extension/low_aero_offset_fit[{}]".format(index), "0-5 timing-center source selected without a valid primary fit")
        if source == "low_aero_0_6_fit":
            if primary_valid:
                _append_issue(issues, "python_tof_extension/low_aero_offset_fit[{}]".format(index), "0-6 fallback selected even though primary fit is valid")
            if not fallback_valid:
                _append_issue(issues, "python_tof_extension/low_aero_offset_fit[{}]".format(index), "0-6 timing-center source selected without a valid fallback fit")
        if source == "stable_global_center_fallback":
            if primary_valid or fallback_valid:
                _append_issue(
                    issues,
                    "python_tof_extension/low_aero_offset_fit[{}]".format(index),
                    "stable center fallback selected even though an offset fit is valid",
                    {"primary_valid": primary_valid, "fallback_valid": fallback_valid},
                )
            if bool(row.get("offset_refinement_valid", False)) or bool(row.get("offset_refinement_applied", False)):
                _append_issue(
                    issues,
                    "python_tof_extension/low_aero_offset_fit[{}]".format(index),
                    "stable center fallback marked as applied/valid offset refinement",
                )
            if not bool(row.get("timing_center_model_valid", False)):
                _append_issue(
                    issues,
                    "python_tof_extension/low_aero_offset_fit[{}]".format(index),
                    "stable center fallback did not produce a valid timing-center model",
                )
    if not isinstance(cfg, dict) or not cfg:
        _append_issue(issues, "python_tof_extension/low_aero", "missing serialized low_aero_offset_config")


def _validate_offset_rows(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
) -> None:
    diag = _diagnostics(payload)
    rows = list(payload.get("delta_timing_offset_fits") or diag.get("delta_timing_offset_fits") or [])
    if not rows:
        _append_issue(issues, "python_tof_extension/delta_offsets", "missing delta offset fit rows")
        return
    for row in rows:
        if not isinstance(row, dict):
            continue
        delta_index = int(row.get("delta_index", -1))
        scope = "python_tof_extension/delta_offsets[{}]".format(delta_index)
        center_source = str(row.get("selected_timing_center_source") or row.get("timing_center_source") or "")
        if center_source == "stable_global_center_fallback":
            if bool(row.get("valid", False)):
                _append_issue(issues, scope, "stable center fallback row is marked as a valid offset fit")
            if bool(row.get("offset_refinement_valid", False)) or bool(row.get("offset_refinement_applied", False)):
                _append_issue(issues, scope, "stable center fallback row is marked as applied/valid offset refinement")
            if not bool(row.get("timing_center_model_valid", False)):
                _append_issue(issues, scope, "stable center fallback row lacks a valid timing-center model")
            if not row.get("reference_kaon_mean") and _float_or_none(row.get("reference_kaon_mean")) is None:
                _append_issue(issues, scope, "stable center fallback row missing reference kaon mean")
            if not row.get("reference_proton_mean") and _float_or_none(row.get("reference_proton_mean")) is None:
                _append_issue(issues, scope, "stable center fallback row missing reference proton mean")
            continue
        cfg = row.get("tof_offset_validation") or {}
        for key in (
            "maximum_offset_error_ns",
            "maximum_chi2_ndf",
            "minimum_component_significance",
            "minimum_smaller_component_fraction",
        ):
            if key not in cfg:
                _append_issue(issues, scope, "missing serialized offset threshold", {"key": key})
        if row.get("smaller_component_fraction_definition") != "min(K_amp,p_amp)/(K_amp+p_amp)":
            _append_issue(
                issues,
                scope,
                "wrong smaller-component fraction definition",
                {"definition": row.get("smaller_component_fraction_definition")},
            )
        valid = bool(row.get("valid", False))
        if valid:
            if int(row.get("fit_status_code", -999) or -999) != 0:
                _append_issue(issues, scope, "valid offset fit does not have ROOT status 0")
            offset_error = _float_or_none(row.get("delta_offset_error"))
            chi2_ndf = _float_or_none(row.get("chi2_ndf"))
            max_error = _float_or_none(cfg.get("maximum_offset_error_ns"))
            max_chi2 = _float_or_none(cfg.get("maximum_chi2_ndf"))
            if offset_error is None or (max_error is not None and offset_error > max_error):
                _append_issue(issues, scope, "valid offset error violates serialized threshold")
            if chi2_ndf is None or (max_chi2 is not None and chi2_ndf > max_chi2):
                _append_issue(issues, scope, "valid offset chi2/ndf violates serialized threshold")
            if row.get("rejection_reasons"):
                _append_issue(issues, scope, "valid offset row has rejection_reasons", {"reasons": row.get("rejection_reasons")})
        elif not row.get("rejection_reasons"):
            _append_issue(issues, scope, "invalid offset row has no rejection_reasons")


def _validate_timing_constraints(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
    center_tolerance: float = 1.0e-6,
) -> None:
    diag = _diagnostics(payload)
    probe_kind = str(diag.get("selected_probe_kind") or (payload.get("settings") or {}).get("selected_timing_probe_kind") or "ct")
    beam_spacing = _float_or_none(((payload.get("settings") or {}).get("global_fit") or {}).get("beam_bunch_spacing_ns"))
    if beam_spacing is None:
        beam_spacing = 2.0
    delta_slices = list(payload.get("delta_slice_fits") or diag.get("delta_slice_fits") or [])
    global_shapes = list(payload.get("global_shapes") or [])
    for delta_index, slice_rows in enumerate(delta_slices):
        for aero_index, row in enumerate(slice_rows or []):
            if not isinstance(row, dict) or not bool(row.get("valid", False)):
                continue
            scope = "python_tof_extension/delta_slice_fits[{}][{}]".format(delta_index, aero_index)
            k_mean = _float_or_none(row.get("predicted_kaon_mean"))
            p_mean = _float_or_none(row.get("predicted_proton_mean"))
            p_mean_raw = _float_or_none(row.get("predicted_proton_mean_raw", row.get("predicted_proton_mean")))
            delta_t = _float_or_none(row.get("mean_delta_t_pk_ns"))
            if k_mean is None or p_mean is None or delta_t is None:
                _append_issue(issues, scope, "missing timing-constraint mean fields")
                continue
            if probe_kind == "rf":
                distance = _modulo_distance(abs(p_mean - k_mean), delta_t, beam_spacing)
            else:
                distance = abs(abs(p_mean - k_mean) - delta_t)
            if distance > center_tolerance:
                _append_issue(
                    issues,
                    scope,
                    "timing-constraint center inconsistency",
                    {
                        "probe_kind": probe_kind,
                        "distance": distance,
                        "tolerance": center_tolerance,
                        "k_mean": k_mean,
                        "p_mean": p_mean,
                        "mean_delta_t_pk_ns": delta_t,
                    },
                )
            for sigma_key, reference_key in (
                ("kaon_sigma", "reference_global_kaon_sigma"),
                ("proton_sigma", "reference_global_proton_sigma"),
            ):
                sigma = _float_or_none(row.get(sigma_key))
                reference_sigma = _float_or_none(row.get(reference_key))
                if reference_sigma is None and 0 <= aero_index < len(global_shapes):
                    reference_sigma = _float_or_none((global_shapes[aero_index] or {}).get(sigma_key))
                if reference_sigma is None:
                    continue
                if sigma is None or abs(sigma - reference_sigma) > center_tolerance:
                    _append_issue(
                        issues,
                        scope,
                        "fixed-width consistency failure",
                        {"sigma_key": sigma_key, "sigma": sigma, "reference": reference_sigma},
                    )


def _validate_invalid_offset_propagation(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
) -> None:
    diag = _diagnostics(payload)
    offsets = list(payload.get("delta_timing_offset_fits") or diag.get("delta_timing_offset_fits") or [])
    delta_slices = list(payload.get("delta_slice_fits") or [])
    support = list(payload.get("support_by_delta") or [])
    for row in offsets:
        if not isinstance(row, dict) or bool(row.get("valid", False)):
            continue
        delta_index = int(row.get("delta_index", -1))
        if delta_index < 0:
            continue
        scope = "python_tof_extension/invalid_offset_propagation[{}]".format(delta_index)
        if bool(row.get("timing_center_model_valid", False)):
            source = str(row.get("selected_timing_center_source") or row.get("timing_center_source") or "")
            if source != "stable_global_center_fallback":
                _append_issue(
                    issues,
                    scope,
                    "invalid offset has timing-center model from unexpected source",
                    {"selected_timing_center_source": source},
                )
            continue
        if delta_index < len(support) and str(support[delta_index]) != "unsupported":
            _append_issue(issues, scope, "invalid offset did not force unsupported support label", {"support": support[delta_index]})
        for aero_index, slice_row in enumerate(delta_slices[delta_index] if delta_index < len(delta_slices) else []):
            if isinstance(slice_row, dict) and bool(slice_row.get("valid", False)):
                _append_issue(issues, scope, "invalid offset produced valid slice fit", {"aero_index": aero_index})


def _validate_cell_fit_attempts(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
) -> None:
    diag = _diagnostics(payload)
    debug_rows = list(diag.get("delta_support_debug_rows") or [])
    for delta_index, row in enumerate(debug_rows):
        if not isinstance(row, dict):
            continue
        if not bool(row.get("timing_center_model_valid", False)):
            continue
        candidate_cells = 0
        attempted_cells = 0
        for aero_index, slice_row in enumerate(row.get("slice_rows") or []):
            if not isinstance(slice_row, dict):
                continue
            if not bool(slice_row.get("global_shape_valid", False)):
                continue
            if not bool(slice_row.get("timing_center_model_valid", False)):
                continue
            support_entries = int(slice_row.get("support_entries", 0) or 0)
            minimum_required = int(slice_row.get("minimum_required_entries", 0) or 0)
            if minimum_required > 0 and support_entries < minimum_required:
                continue
            candidate_cells += 1
            if bool(slice_row.get("fit_attempted", False)):
                attempted_cells += 1
                continue
            _append_issue(
                issues,
                "python_tof_extension/cell_fit_attempts[{}][{}]".format(delta_index, aero_index),
                "cell had timing model and support but fit was not attempted",
                {
                    "fit_status": slice_row.get("fit_status"),
                    "support_entries": support_entries,
                    "minimum_required_entries": minimum_required,
                    "timing_center_source": slice_row.get("timing_center_source"),
                },
            )
            if str(slice_row.get("fit_status") or "") == "insufficient_support":
                _append_issue(
                    issues,
                    "python_tof_extension/cell_fit_attempts[{}][{}]".format(delta_index, aero_index),
                    "supported cell was incorrectly marked insufficient_support",
                    {
                        "support_entries": support_entries,
                        "minimum_required_entries": minimum_required,
                    },
                )
        if candidate_cells > 0 and attempted_cells <= 0:
            _append_issue(
                issues,
                "python_tof_extension/cell_fit_attempts[{}]".format(delta_index),
                "delta bin has supported timing cells but no cell fits were attempted",
                {
                    "candidate_cells": candidate_cells,
                    "timing_center_source": row.get("selected_timing_center_source"),
                    "cell_fit_attempt_count": row.get("cell_fit_attempt_count"),
                },
            )


def _validate_cell_fit_statuses_and_counters(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
) -> None:
    expected_reason_by_status = {
        "insufficient_support": "insufficient_entries",
        "missing_histogram": "missing_histogram",
        "invalid_global_shape": "invalid_global_shape",
        "invalid_timing_center_model": "invalid_timing_center_model",
    }
    counter_by_status = {
        "insufficient_support": "cell_fit_skipped_insufficient_support_count",
        "missing_histogram": "cell_fit_skipped_missing_histogram_count",
        "invalid_global_shape": "cell_fit_skipped_invalid_global_shape_count",
        "invalid_timing_center_model": "cell_fit_skipped_invalid_timing_center_model_count",
    }
    diag = _diagnostics(payload)
    debug_rows = list(diag.get("delta_support_debug_rows") or [])
    for delta_index, row in enumerate(debug_rows):
        if not isinstance(row, dict):
            continue
        slice_rows = [entry for entry in (row.get("slice_rows") or []) if isinstance(entry, dict)]
        attempted_count = sum(bool(entry.get("fit_attempted", False)) for entry in slice_rows)
        skipped_count = len(slice_rows) - attempted_count
        status_counts = {key: 0 for key in counter_by_status.values()}
        other_count = 0
        for aero_index, slice_row in enumerate(slice_rows):
            if bool(slice_row.get("fit_attempted", False)):
                continue
            status = str(slice_row.get("fit_status") or "")
            reasons = set(str(reason) for reason in (slice_row.get("rejection_reasons") or []))
            if not reasons and slice_row.get("rejection_reason"):
                reasons = {
                    reason.strip()
                    for reason in str(slice_row.get("rejection_reason")).split(";")
                    if reason.strip()
                }
            required_reason = expected_reason_by_status.get(status)
            scope = "python_tof_extension/cell_status[{}][{}]".format(delta_index, aero_index)
            if required_reason and required_reason not in reasons:
                _append_issue(
                    issues,
                    scope,
                    "pre-fit status does not have its required rejection reason",
                    {"fit_status": status, "reasons": sorted(reasons), "required_reason": required_reason},
                )
            counter_key = counter_by_status.get(status)
            if counter_key:
                status_counts[counter_key] += 1
            else:
                other_count += 1
        scope = "python_tof_extension/cell_counters[{}]".format(delta_index)
        if int(row.get("cell_fit_attempt_count", 0) or 0) != attempted_count:
            _append_issue(issues, scope, "cell-fit attempt counter mismatch", {"stored": row.get("cell_fit_attempt_count"), "actual": attempted_count})
        if int(row.get("cell_fit_skipped_count", 0) or 0) != skipped_count:
            _append_issue(issues, scope, "cell-fit skipped counter mismatch", {"stored": row.get("cell_fit_skipped_count"), "actual": skipped_count})
        if attempted_count + skipped_count != len(slice_rows):
            _append_issue(issues, scope, "cell-fit attempts and skips do not cover all aerogel cells")
        for counter_key, actual_count in status_counts.items():
            if int(row.get(counter_key, 0) or 0) != actual_count:
                _append_issue(issues, scope, "cell-fit skipped-status counter mismatch", {"counter": counter_key, "stored": row.get(counter_key), "actual": actual_count})
        if int(row.get("cell_fit_skipped_other_count", 0) or 0) != other_count:
            _append_issue(issues, scope, "cell-fit skipped-other counter mismatch", {"stored": row.get("cell_fit_skipped_other_count"), "actual": other_count})
        expected_skips = sum(status_counts.values()) + other_count
        if skipped_count != expected_skips:
            _append_issue(issues, scope, "cell-fit skipped subcategories do not sum to skipped total", {"skipped": skipped_count, "categorized": expected_skips})


def _validate_weak_high_aero_components(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
    abs_tol: float,
) -> None:
    aero_edges = [float(value) for value in (payload.get("aero_edges") or (0.0, 3.0, 6.0, 10.0, 15.0, 25.0))]
    delta_slices = list(payload.get("delta_slice_fits") or [])
    for delta_index, slice_rows in enumerate(delta_slices):
        for aero_index, row in enumerate(slice_rows or []):
            if not isinstance(row, dict):
                continue
            if aero_index + 1 < len(aero_edges):
                aero_center = 0.5 * (aero_edges[aero_index] + aero_edges[aero_index + 1])
            else:
                aero_center = float(aero_index)
            if aero_center < 10.0:
                continue
            scope = "python_tof_extension/high_aero_weak_component[{}][{}]".format(delta_index, aero_index)
            proton_yield = _float_or_none(row.get("proton_yield"))
            detected = bool(row.get("proton_component_detected", False))
            below = bool(row.get("proton_component_below_significance", False))
            if bool(row.get("valid", False)) and below and detected:
                _append_issue(issues, scope, "cell cannot be both below-significance and proton-detected")
            if bool(row.get("valid", False)) and below and proton_yield is not None and abs(proton_yield) > max(abs_tol, 1.0e-12):
                _append_issue(
                    issues,
                    scope,
                    "below-significance high-aerogel proton component was not zeroed",
                    {"proton_yield": proton_yield},
                )
            if not detected and proton_yield is not None and abs(proton_yield) > max(abs_tol, 1.0e-12):
                _append_issue(
                    issues,
                    scope,
                    "undetected proton component has nonzero effective proton yield",
                    {"proton_yield": proton_yield},
                )
    diag = _diagnostics(payload)
    warnings = [str(value) for value in (diag.get("warnings") or [])]
    warning_name = "high_aero_proton_fraction_exceeds_low_aero"
    if warning_name in warnings:
        fallback_reason = str(payload.get("fallback_reason") or "")
        if warning_name in fallback_reason:
            _append_issue(
                issues,
                "python_tof_extension/high_aero_warning",
                "high-aerogel warning propagated into fallback reason",
                {"fallback_reason": fallback_reason},
            )


def _validate_closure_and_weights(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
    python_weight_payloads: list[dict[str, Any]] | None,
    abs_tol: float,
) -> None:
    diag = _diagnostics(payload)
    for row in diag.get("event_weight_closure_by_delta") or []:
        if not isinstance(row, dict):
            continue
        delta_index = int(row.get("delta_index", -1))
        fitted = _float_or_none(row.get("fitted_proton_yield"))
        summed = _float_or_none(row.get("summed_event_proton_probability"))
        ratio = _float_or_none(row.get("closure_ratio"))
        if fitted is None or summed is None:
            continue
        expected = summed / fitted if fitted != 0.0 else None
        if expected is not None and ratio is not None and abs(ratio - expected) > max(abs_tol, 1.0e-12):
            _append_issue(
                issues,
                "python_tof_extension/event_weight_closure[{}]".format(delta_index),
                "closure ratio does not match summed/fitted algebra",
                {"ratio": ratio, "expected": expected},
            )
    for row in python_weight_payloads or []:
        if not isinstance(row, dict):
            continue
        signature = row.get("signature", "unknown")
        for key in ("proton_weight", "cleaned_factor", "final_cleaned_factor"):
            value = _float_or_none(row.get(key))
            if value is None:
                continue
            if value < -abs_tol or value > 1.0 + abs_tol:
                _append_issue(
                    issues,
                    "python_tof_extension/event_weights",
                    "{} out of [0, 1]".format(key),
                    {"signature": signature, "value": value},
                )


def _validate_branch_alias_metadata(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
) -> None:
    diag = _diagnostics(payload)
    aliases = set(str(value) for value in diag.get("tof_required_aliases") or [])
    forbidden = set(str(value) for value in diag.get("tof_forbidden_replay_names") or [])
    required = {"P_gtr_p", "ssxptar", "ssyptar", "ssdelta"}
    missing = sorted(required - aliases)
    if missing:
        _append_issue(
            issues,
            "python_tof_extension/branch_aliases",
            "missing required alias metadata",
            {"missing": missing},
        )
    if aliases & forbidden:
        _append_issue(
            issues,
            "python_tof_extension/branch_aliases",
            "forbidden replay names listed as required aliases",
            {"overlap": sorted(aliases & forbidden)},
        )


def _validate_k_lambda_reference_persistence(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
) -> None:
    diag = _diagnostics(payload)
    input_loaded = payload.get("k_lambda_simc_input_loaded", diag.get("k_lambda_simc_input_loaded"))
    if input_loaded is None:
        return
    scope = "component_subtraction/k_lambda_reference"
    if not bool(input_loaded):
        if str(payload.get("k_lambda_simc_reference_source") or diag.get("k_lambda_simc_reference_source") or "") == "immutable_aligned_k_lambda_simc":
            _append_issue(issues, scope, "immutable K-Lambda reference source set even though input was not loaded")
        return
    available = payload.get("k_lambda_simc_reference_available", diag.get("k_lambda_simc_reference_available"))
    source = str(payload.get("k_lambda_simc_reference_source") or diag.get("k_lambda_simc_reference_source") or "")
    integral = _float_or_none(payload.get("k_lambda_simc_reference_integral", diag.get("k_lambda_simc_reference_integral")))
    if not bool(available):
        _append_issue(issues, scope, "K-Lambda SIMC input was loaded but the immutable comparison reference is not marked available")
    if source != "immutable_aligned_k_lambda_simc":
        _append_issue(issues, scope, "K-Lambda reference is not using the immutable aligned SIMC source", {"source": source})
    if integral is None or integral <= 0.0:
        _append_issue(issues, scope, "K-Lambda immutable reference has missing or non-positive integral", {"integral": integral})


def _validate_python_tof_extension(
    issues: list[dict[str, Any]],
    payload: dict[str, Any],
    python_weight_payloads: list[dict[str, Any]] | None,
    abs_tol: float,
) -> None:
    _validate_tof_summary_rows(issues, payload)
    _validate_low_aero_offset_metadata(issues, payload)
    _validate_offset_rows(issues, payload)
    _validate_timing_constraints(issues, payload, center_tolerance=1.0e-6)
    _validate_invalid_offset_propagation(issues, payload)
    _validate_cell_fit_attempts(issues, payload)
    _validate_cell_fit_statuses_and_counters(issues, payload)
    _validate_weak_high_aero_components(issues, payload, abs_tol)
    _validate_closure_and_weights(issues, payload, python_weight_payloads, abs_tol)
    _validate_branch_alias_metadata(issues, payload)
    _validate_k_lambda_reference_persistence(issues, payload)


def _write_text_report(path: str, report: dict[str, Any]) -> None:
    lines = [
        "Proton-cleaning parity validation",
        "passed={}".format(report["passed"]),
        "validation_scope={}".format(report.get("validation_scope", "combined")),
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
    parser.add_argument(
        "--validation-scope",
        choices=("legacy_parity", "python_tof_extension", "combined"),
        default="combined",
        help="Run only immutable C-owned parity checks, Python TOF-extension consistency checks, or both.",
    )
    args = parser.parse_args()

    legacy_issues: list[dict[str, Any]] = []
    extension_issues: list[dict[str, Any]] = []
    python_payload = _load_json(args.python_json)
    macro_payload = _load_json(args.macro_json)

    if args.validation_scope in ("legacy_parity", "combined"):
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
    python_weight_payloads = None
    if args.python_weights and args.macro_weights:
        python_weight_payloads = _load_weight_payloads(args.python_weights)
        python_weights = _load_weight_map(args.python_weights)
        macro_weights = _load_weight_map(args.macro_weights)

    if args.validation_scope in ("python_tof_extension", "combined"):
        _validate_python_tof_extension(
            extension_issues,
            python_payload,
            python_weight_payloads,
            args.abs_tol,
        )

    extension_observations = {
        "validation_scope": args.validation_scope,
        "strict_c_checks_applied": args.validation_scope in ("legacy_parity", "combined"),
        "python_tof_extension_checks_applied": args.validation_scope in ("python_tof_extension", "combined"),
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
        "validation_scope": args.validation_scope,
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
