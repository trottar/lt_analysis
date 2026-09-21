"""Render the frozen F.6.2 E.8 presentation-only figure library.

This explicit-input sidecar consumes one accepted aggregate JSON artifact.  It
never constructs analysis state or changes a persisted scientific result.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import datetime as _datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from typing import Any

import matplotlib

matplotlib.use("Agg", force=True)
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np


KINEMATIC = "Q4p4W2p74"
INPUT_SHA256 = "5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1"
ARTIFACT_FINGERPRINT = "ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0"
VALIDATION_FINGERPRINT = "7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b"
SCIENTIFIC_SOURCE_COMMIT = "0b37af2a2927b08bdeaf897c545f290b55329cea"
FIX5_RENDERER_SOURCE_COMMIT = "c88ed65cb18ba6a37358897292b77696016312d1"
FIX5_BUNDLE_PROFILE_COMMIT = "b789d203e11f0927deb59ebcae9dc59fe8add4ae"

ARTIFACT_SCHEMA = "pion_hgcer_method_a_acceptance_refinement_validation_artifact/v1"
VALIDATION_SCHEMA = "pion_hgcer_method_a_acceptance_refinement_validation/v1"
MANIFEST_SCHEMA = "pion_hgcer_f6_2_e8_figure_library_manifest/v1"
PDF_BASENAME = (
    "Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-"
    "figure-library.pdf"
)
MANIFEST_BASENAME = (
    "Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-"
    "figure-library-manifest.json"
)

_COMMIT_RE = re.compile(r"^[0-9a-f]{40}$")
_PARENTS = tuple(
    ("{}-{}".format(phi, epsilon), index)
    for phi, epsilon in (
        ("Left", "lowe"),
        ("Left", "highe"),
        ("Center", "lowe"),
        ("Center", "highe"),
        ("Right", "highe"),
    )
    for index in range(3)
)
_ONE_DIMENSIONAL = ("analysis_MM", "SHMS_xptar", "SHMS_yptar")
_JOINTS = (
    ("analysis_MM__SHMS_xptar", "analysis_MM", "SHMS_xptar"),
    ("analysis_MM__SHMS_yptar", "analysis_MM", "SHMS_yptar"),
    ("SHMS_delta__SHMS_xptar", "SHMS_delta", "SHMS_xptar"),
    ("SHMS_delta__SHMS_yptar", "SHMS_delta", "SHMS_yptar"),
)
_POPULATIONS = ("L", "B", "A")
_FORBIDDEN_PERSISTED_KEYS = frozenset(
    {
        "entry_index",
        "event_id",
        "event_identity",
        "application_records",
        "method_a_training_records",
        "correction_factors",
        "raw_shape_factors",
        "in_support_mask",
        "event_corrections",
        "lookup_table",
        "raw_rows",
        "review_data",
    }
)
_FIXED_PDF_DATE = _datetime.datetime(2000, 1, 1, tzinfo=_datetime.timezone.utc)


class FigureLibraryError(ValueError):
    """The persisted E.8 input or explicit output contract is invalid."""


def _mapping(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise FigureLibraryError("{}_invalid".format(label))
    return value


def _sequence(value: object, label: str) -> Sequence[object]:
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence):
        raise FigureLibraryError("{}_invalid".format(label))
    return value


def _finite(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise FigureLibraryError("{}_nonfinite".format(label))
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise FigureLibraryError("{}_nonfinite".format(label)) from exc
    if not math.isfinite(number):
        raise FigureLibraryError("{}_nonfinite".format(label))
    return number


def _integer(value: object, label: str) -> int:
    if isinstance(value, bool):
        raise FigureLibraryError("{}_invalid".format(label))
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise FigureLibraryError("{}_invalid".format(label)) from exc
    if number != value:
        raise FigureLibraryError("{}_invalid".format(label))
    return number


def _strict_edges(value: object, label: str) -> list[float]:
    edges = [_finite(item, "{}_edge".format(label)) for item in _sequence(value, label)]
    if len(edges) < 2 or any(right <= left for left, right in zip(edges, edges[1:])):
        raise FigureLibraryError("{}_not_strictly_increasing".format(label))
    return edges


def _metric(value: object, label: str) -> Mapping[str, object]:
    metric = _mapping(value, label)
    available = metric.get("available")
    if not isinstance(available, bool):
        raise FigureLibraryError("{}_availability_invalid".format(label))
    reason = metric.get("reason")
    if available:
        _finite(metric.get("value"), "{}_value".format(label))
        if reason is not None:
            raise FigureLibraryError("{}_reason_invalid".format(label))
    elif not isinstance(reason, str) or not reason:
        raise FigureLibraryError("{}_reason_missing".format(label))
    return metric


def _interval(value: object, label: str) -> Mapping[str, object]:
    interval = _mapping(value, label)
    available = interval.get("interval_available")
    if not isinstance(available, bool):
        raise FigureLibraryError("{}_availability_invalid".format(label))
    requested = _integer(interval.get("requested_replica_count"), "{}_requested_count".format(label))
    valid = _integer(interval.get("valid_replica_count"), "{}_valid_count".format(label))
    invalid = _integer(interval.get("invalid_replica_count"), "{}_invalid_count".format(label))
    if min(requested, valid, invalid) < 0 or valid + invalid != requested:
        raise FigureLibraryError("{}_replica_counts_invalid".format(label))
    reason = interval.get("reason")
    if available:
        low = _finite(interval.get("ci_low"), "{}_low".format(label))
        high = _finite(interval.get("ci_high"), "{}_high".format(label))
        if high < low or reason is not None:
            raise FigureLibraryError("{}_bounds_invalid".format(label))
    elif (
        not isinstance(reason, str)
        or not reason
        or interval.get("ci_low") is not None
        or interval.get("ci_high") is not None
    ):
        raise FigureLibraryError("{}_reason_missing".format(label))
    return interval


def _population(value: object, shape: tuple[int, ...], label: str) -> Mapping[str, object]:
    population = _mapping(value, label)
    available = population.get("available")
    if not isinstance(available, bool):
        raise FigureLibraryError("{}_availability_invalid".format(label))
    reason = population.get("reason")
    if available and reason is not None:
        raise FigureLibraryError("{}_reason_invalid".format(label))
    if not available and (not isinstance(reason, str) or not reason):
        raise FigureLibraryError("{}_reason_missing".format(label))
    try:
        contents = np.asarray(population.get("unit_area"), dtype=float)
    except (TypeError, ValueError) as exc:
        raise FigureLibraryError("{}_unit_area_invalid".format(label)) from exc
    if not np.all(np.isfinite(contents)) or np.any(contents < 0.0):
        raise FigureLibraryError("{}_unit_area_invalid".format(label))
    if available:
        if contents.shape != shape:
            raise FigureLibraryError("{}_unit_area_invalid".format(label))
        if not math.isclose(float(np.sum(contents)), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
            raise FigureLibraryError("{}_unit_area_not_normalized".format(label))
    else:
        stored_shape = shape if len(shape) == 1 else (math.prod(shape),)
        if contents.shape != stored_shape or not np.all(contents == 0.0):
            raise FigureLibraryError("{}_unavailable_placeholder_invalid".format(label))
    return population


def _one_dimensional(value: object, label: str) -> Mapping[str, object]:
    payload = _mapping(value, label)
    edges = _strict_edges(payload.get("edges"), "{}_edges".format(label))
    for population in _POPULATIONS:
        _population(payload.get(population), (len(edges) - 1,), "{}_{}".format(label, population))
    metrics = _mapping(payload.get("metrics"), "{}_metrics".format(label))
    _metric(metrics.get("DeltaH"), "{}_DeltaH".format(label))
    _metric(metrics.get("kappa"), "{}_kappa".format(label))
    return payload


def _joint(value: object, expected_x: str, expected_y: str, label: str) -> Mapping[str, object]:
    payload = _mapping(value, label)
    if payload.get("x_variable") != expected_x or payload.get("y_variable") != expected_y:
        raise FigureLibraryError("{}_variables_invalid".format(label))
    x_edges = _strict_edges(payload.get("x_edges"), "{}_x_edges".format(label))
    y_edges = _strict_edges(payload.get("y_edges"), "{}_y_edges".format(label))
    for population in _POPULATIONS:
        _population(
            payload.get(population),
            (len(x_edges) - 1, len(y_edges) - 1),
            "{}_{}".format(label, population),
        )
    metrics = _mapping(payload.get("metrics"), "{}_metrics".format(label))
    _metric(metrics.get("kappa"), "{}_kappa".format(label))
    return payload


def _require_flags(value: Mapping[str, object], expected: Mapping[str, bool], label: str) -> None:
    for name, expectation in expected.items():
        if value.get(name) is not expectation:
            raise FigureLibraryError("{}_{}_invalid".format(label, name))


def _reject_event_payload(value: object) -> None:
    if isinstance(value, Mapping):
        for name, item in value.items():
            if name in _FORBIDDEN_PERSISTED_KEYS:
                raise FigureLibraryError("aggregate_payload_{}_forbidden".format(name))
            _reject_event_payload(item)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        for item in value:
            _reject_event_payload(item)


def _validate_child(child_value: object, parent: Mapping[str, object], position: int) -> Mapping[str, object]:
    child = _mapping(child_value, "child")
    if (
        child.get("setting_id") != parent.get("setting_id")
        or _integer(child.get("canonical_t_index"), "child_t_index") != parent.get("canonical_t_index")
        or _integer(child.get("phi_index"), "child_phi_index") != position
    ):
        raise FigureLibraryError("child_geometry_invalid")
    phi_edges = _strict_edges(parent.get("phi_edges"), "parent_phi_edges")
    if len(phi_edges) != 10:
        raise FigureLibraryError("parent_phi_inventory_invalid")
    if (
        _finite(child.get("phi_low"), "child_phi_low") != phi_edges[position]
        or _finite(child.get("phi_high"), "child_phi_high") != phi_edges[position + 1]
    ):
        raise FigureLibraryError("child_phi_bounds_invalid")
    availability = _mapping(child.get("availability"), "child_availability")
    empty = availability.get("completely_empty")
    if not isinstance(empty, bool):
        raise FigureLibraryError("child_empty_flag_invalid")
    for name in ("has_low_response", "has_prompt_control", "has_full_application"):
        if not isinstance(availability.get(name), bool):
            raise FigureLibraryError("child_{}_invalid".format(name))
    counts = _mapping(child.get("population_counts"), "child_population_counts")
    for name in ("N_low", "N_control", "N_full_application"):
        if _integer(counts.get(name), "child_{}".format(name)) < 0:
            raise FigureLibraryError("child_{}_invalid".format(name))
    support = _mapping(child.get("support"), "child_support")
    for name in ("prompt_control", "full_physical_application"):
        fraction = _mapping(support.get(name), "child_{}_support".format(name)).get("ood_fraction")
        if fraction is not None and not 0.0 <= _finite(fraction, "child_{}_ood".format(name)) <= 1.0:
            raise FigureLibraryError("child_{}_ood_invalid".format(name))
    neff = _mapping(child.get("effective_sample_size"), "child_effective_sample_size")
    _metric(neff.get("baseline_w0"), "child_neff_baseline")
    _metric(neff.get("method_a_w0_times_C"), "child_neff_method_a")
    shapes = _mapping(child.get("one_dimensional"), "child_one_dimensional")
    for name in _ONE_DIMENSIONAL:
        _one_dimensional(shapes.get(name), "child_{}".format(name))
    joint = _mapping(child.get("joint_distributions"), "child_joint_distributions")
    for name, x_name, y_name in _JOINTS:
        _joint(joint.get(name), x_name, y_name, "child_{}".format(name))
    signed = _mapping(child.get("signed_background"), "child_signed_background")
    window = _mapping(signed.get("kaon_window"), "child_kaon_window")
    for name in ("P_B_K", "P_A_K", "DeltaP_K", "f_refine_K"):
        _metric(window.get(name), "child_{}".format(name))
    variance = _mapping(signed.get("variance_proxy"), "child_variance_proxy")
    _metric(variance.get("R_V"), "child_R_V")
    bootstrap = _mapping(child.get("bootstrap"), "child_bootstrap")
    bootstrap_one = _mapping(bootstrap.get("one_dimensional"), "child_bootstrap_one_dimensional")
    bootstrap_mm = _mapping(bootstrap_one.get("analysis_MM"), "child_bootstrap_analysis_MM")
    _interval(bootstrap_mm.get("DeltaH"), "child_bootstrap_DeltaH")
    _interval(bootstrap_mm.get("kappa"), "child_bootstrap_kappa")
    bootstrap_joint = _mapping(bootstrap.get("joint_missing_mass_acceptance"), "child_bootstrap_joint")
    for name in ("analysis_MM__SHMS_xptar", "analysis_MM__SHMS_yptar"):
        _interval(
            _mapping(bootstrap_joint.get(name), "child_bootstrap_{}".format(name)).get("kappa"),
            "child_bootstrap_{}_kappa".format(name),
        )
    bootstrap_window = _mapping(bootstrap.get("kaon_window"), "child_bootstrap_window")
    _interval(bootstrap_window.get("DeltaP_K"), "child_bootstrap_DeltaP_K")
    return child


def validate_artifact(artifact_value: object) -> Mapping[str, object]:
    """Validate the persisted-only F.6.2 payload before any rendering."""

    artifact = _mapping(artifact_value, "artifact")
    if artifact.get("schema_version") != ARTIFACT_SCHEMA:
        raise FigureLibraryError("artifact_schema_invalid")
    _require_flags(
        artifact,
        {
            "non_authoritative": True,
            "validation_only": True,
            "manual_review_required": True,
            "production_application_performed": False,
            "production_objects_mutated": False,
            "yield_constructed": False,
            "cross_section_constructed": False,
            "method_b_numerical_dependency": False,
            "automatic_case_classification": False,
        },
        "artifact",
    )
    if artifact.get("artifact_fingerprint") != ARTIFACT_FINGERPRINT:
        raise FigureLibraryError("artifact_fingerprint_invalid")
    validation = _mapping(artifact.get("validation"), "validation")
    if validation.get("schema_version") != VALIDATION_SCHEMA or validation.get("available") is not True:
        raise FigureLibraryError("validation_schema_or_availability_invalid")
    if validation.get("fingerprint") != VALIDATION_FINGERPRINT:
        raise FigureLibraryError("validation_fingerprint_invalid")
    _require_flags(
        validation,
        {
            "non_authoritative": True,
            "validation_only": True,
            "event_correction_persisted": False,
            "production_application_performed": False,
            "production_objects_mutated": False,
            "yield_constructed": False,
            "cross_section_constructed": False,
            "root_object_constructed": False,
            "child_renormalization_performed": False,
            "smoothing_or_interpolation_performed": False,
            "absolute_probability_constructed": False,
            "method_b_numerical_dependency": False,
            "automatic_case_classification": False,
            "case_thresholds_defined": False,
            "final_yield_uncertainty_claimed": False,
        },
        "validation",
    )
    window = _mapping(validation.get("kaon_window"), "kaon_window")
    if _finite(window.get("mm_max"), "kaon_window_high") <= _finite(window.get("mm_min"), "kaon_window_low"):
        raise FigureLibraryError("kaon_window_invalid")
    parents = _sequence(validation.get("parents"), "parents")
    if len(parents) != len(_PARENTS):
        raise FigureLibraryError("parent_inventory_invalid")
    observed_parents: set[tuple[str, int]] = set()
    children: list[Mapping[str, object]] = []
    for parent_value in parents:
        parent = _mapping(parent_value, "parent")
        setting = parent.get("setting_id")
        if not isinstance(setting, str):
            raise FigureLibraryError("parent_setting_invalid")
        key = (setting, _integer(parent.get("canonical_t_index"), "parent_t_index"))
        if key not in _PARENTS or key in observed_parents:
            raise FigureLibraryError("parent_geometry_invalid")
        observed_parents.add(key)
        parent_children = _sequence(parent.get("children"), "parent_children")
        if len(parent_children) != 9:
            raise FigureLibraryError("parent_child_inventory_invalid")
        for position, child_value in enumerate(parent_children):
            children.append(_validate_child(child_value, parent, position))
    if observed_parents != set(_PARENTS) or len(children) != 135:
        raise FigureLibraryError("canonical_inventory_invalid")
    empty_count = sum(bool(_mapping(child.get("availability"), "child_availability").get("completely_empty")) for child in children)
    if empty_count != 20 or len(children) - empty_count != 115:
        raise FigureLibraryError("rendered_empty_inventory_invalid")
    _reject_event_payload(artifact)
    return validation


def _metric_text(metric: object, interval: object | None = None) -> str:
    value = _mapping(metric, "display_metric")
    if value.get("available") is not True:
        return "unavailable:{}".format(value.get("reason"))
    text = "{:.5g}".format(_finite(value.get("value"), "display_metric_value"))
    if interval is not None:
        bounds = _mapping(interval, "display_interval")
        if bounds.get("interval_available") is True:
            text += " 95% CI [{:.5g}, {:.5g}]".format(
                _finite(bounds.get("ci_low"), "display_interval_low"),
                _finite(bounds.get("ci_high"), "display_interval_high"),
            )
        else:
            text += " ({})".format(bounds.get("reason"))
    return text


def _fraction_text(value: object) -> str:
    if value is None:
        return "unavailable:empty_population"
    return "{:.5g}".format(_finite(value, "display_fraction"))


def _unavailable_text(population: Mapping[str, object]) -> str:
    return "unavailable:{}".format(population.get("reason"))


def _draw_overlay(axis: Any, payload: Mapping[str, object], variable: str, kaon_window: Mapping[str, object]) -> None:
    edges = np.asarray(payload["edges"], dtype=float)
    unavailable: list[str] = []
    for name, color in (("L", "#1f77b4"), ("B", "#333333"), ("A", "#7f3c8d")):
        population = _mapping(payload[name], "plot_population")
        if population["available"] is True:
            axis.stairs(np.asarray(population["unit_area"], dtype=float), edges, label=name, color=color, linewidth=1.4)
        else:
            unavailable.append("{} {}".format(name, _unavailable_text(population)))
    if variable == "analysis_MM":
        axis.axvspan(float(kaon_window["mm_min"]), float(kaon_window["mm_max"]), color="#d9d9d9", alpha=0.65, label="frozen kaon window")
    if unavailable:
        axis.text(0.02, 0.97, "\n".join(unavailable), transform=axis.transAxes, va="top", fontsize=6.5, wrap=True)
    axis.set_title(variable, fontsize=9)
    axis.set_xlabel(variable)
    axis.set_ylabel("Normalized bin fraction")
    axis.legend(fontsize=6.5, loc="best")


def _row_vmax(payload: Mapping[str, object]) -> float:
    values = []
    for name in _POPULATIONS:
        population = _mapping(payload[name], "plot_population")
        if population["available"] is True:
            values.append(float(np.max(np.asarray(population["unit_area"], dtype=float))))
    return max(values or [1.0])


def _draw_joint(axis: Any, payload: Mapping[str, object], population_name: str, vmax: float) -> Any | None:
    population = _mapping(payload[population_name], "plot_population")
    axis.set_title("{} — {}".format(payload["x_variable"] + " × " + payload["y_variable"], population_name), fontsize=8)
    axis.set_xlabel(str(payload["x_variable"]))
    axis.set_ylabel(str(payload["y_variable"]))
    if population["available"] is not True:
        axis.text(0.5, 0.5, _unavailable_text(population), transform=axis.transAxes, ha="center", va="center", fontsize=8, wrap=True)
        return None
    mesh = axis.pcolormesh(
        np.asarray(payload["x_edges"], dtype=float),
        np.asarray(payload["y_edges"], dtype=float),
        np.asarray(population["unit_area"], dtype=float).T,
        shading="auto",
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
    )
    return mesh


def _footer_lines(child: Mapping[str, object]) -> list[str]:
    counts = _mapping(child["population_counts"], "footer_counts")
    support = _mapping(child["support"], "footer_support")
    neff = _mapping(child["effective_sample_size"], "footer_neff")
    shapes = _mapping(child["one_dimensional"], "footer_shapes")
    mm = _mapping(shapes["analysis_MM"], "footer_mm")
    mm_metrics = _mapping(mm["metrics"], "footer_mm_metrics")
    joints = _mapping(child["joint_distributions"], "footer_joints")
    x_joint = _mapping(joints["analysis_MM__SHMS_xptar"], "footer_x_joint")
    y_joint = _mapping(joints["analysis_MM__SHMS_yptar"], "footer_y_joint")
    signed = _mapping(child["signed_background"], "footer_signed")
    kaon = _mapping(signed["kaon_window"], "footer_kaon")
    variance = _mapping(signed["variance_proxy"], "footer_variance")
    bootstrap = _mapping(child["bootstrap"], "footer_bootstrap")
    bootstrap_one = _mapping(bootstrap["one_dimensional"], "footer_bootstrap_one")
    bootstrap_mm = _mapping(bootstrap_one["analysis_MM"], "footer_bootstrap_mm")
    bootstrap_joint = _mapping(bootstrap["joint_missing_mass_acceptance"], "footer_bootstrap_joint")
    bootstrap_window = _mapping(bootstrap["kaon_window"], "footer_bootstrap_window")
    return [
        "N_low/N_control/N_full_application={}/{}/{}; N_eff(B/A)={}/{}; OOD(prompt/full)={}/{}".format(
            counts["N_low"], counts["N_control"], counts["N_full_application"],
            _metric_text(neff["baseline_w0"]), _metric_text(neff["method_a_w0_times_C"]),
            _fraction_text(_mapping(support["prompt_control"], "footer_prompt").get("ood_fraction")),
            _fraction_text(_mapping(support["full_physical_application"], "footer_full").get("ood_fraction")),
        ),
        "analysis_MM: DeltaH {}; kappa {}; MM×xptar kappa {}; MM×yptar kappa {}".format(
            _metric_text(mm_metrics["DeltaH"], bootstrap_mm["DeltaH"]),
            _metric_text(mm_metrics["kappa"], bootstrap_mm["kappa"]),
            _metric_text(_mapping(x_joint["metrics"], "footer_x_metrics")["kappa"], _mapping(bootstrap_joint["analysis_MM__SHMS_xptar"], "footer_x_bootstrap")["kappa"]),
            _metric_text(_mapping(y_joint["metrics"], "footer_y_metrics")["kappa"], _mapping(bootstrap_joint["analysis_MM__SHMS_yptar"], "footer_y_bootstrap")["kappa"]),
        ),
        "Kaon window: P_B^K {}; P_A^K {}; DeltaP^K {}; f_refine^K {}; R_V {}".format(
            _metric_text(kaon["P_B_K"]), _metric_text(kaon["P_A_K"]),
            _metric_text(kaon["DeltaP_K"], bootstrap_window["DeltaP_K"]),
            _metric_text(kaon["f_refine_K"]), _metric_text(variance["R_V"]),
        ),
    ]


def _cover_page(pdf: PdfPages, input_name: str, renderer_commit: str, profile_commit: str) -> None:
    figure, axis = plt.subplots(figsize=(14, 10))
    axis.axis("off")
    lines = [
        "KaonLT E.8 — F.6.2 frozen figure library",
        "Presentation-only sidecar; no production or scientific acceptance is altered.",
        "Input JSON: {}".format(input_name),
        "Input SHA-256: {}".format(INPUT_SHA256),
        "Artifact fingerprint: {}".format(ARTIFACT_FINGERPRINT),
        "Validation fingerprint: {}".format(VALIDATION_FINGERPRINT),
        "F.6.2 scientific source: {}".format(SCIENTIFIC_SOURCE_COMMIT),
        "Fix.5 upstream renderer source: {}".format(FIX5_RENDERER_SOURCE_COMMIT),
        "Fix.5 upstream bundle/profile: {}".format(FIX5_BUNDLE_PROFILE_COMMIT),
        "E.8 renderer source: {}".format(renderer_commit),
        "E.8 bundle/profile: {}".format(profile_commit),
        "Inventory: 15 parents; 135 canonical children; 115 rendered; 20 explicitly empty; 116 PDF pages.",
    ]
    axis.text(0.05, 0.93, lines[0], transform=axis.transAxes, fontsize=19, fontweight="bold", va="top")
    axis.text(0.05, 0.83, "\n\n".join(lines[1:]), transform=axis.transAxes, fontsize=10, va="top", family="monospace")
    pdf.savefig(figure)
    plt.close(figure)


def _child_page(pdf: PdfPages, child: Mapping[str, object], kaon_window: Mapping[str, object]) -> None:
    figure, axes = plt.subplots(3, 3, figsize=(15, 11))
    figure.suptitle(
        "E.8 {} t{} phi{} [{:.0f}, {:.0f})".format(
            child["setting_id"], child["canonical_t_index"], child["phi_index"],
            float(child["phi_low"]), float(child["phi_high"]),
        ),
        fontsize=13,
        fontweight="bold",
    )
    shapes = _mapping(child["one_dimensional"], "plot_shapes")
    for axis, variable in zip(axes[0], _ONE_DIMENSIONAL):
        _draw_overlay(axis, _mapping(shapes[variable], "plot_shape"), variable, kaon_window)
    joints = _mapping(child["joint_distributions"], "plot_joints")
    for row, key in zip(axes[1:], ("SHMS_delta__SHMS_xptar", "SHMS_delta__SHMS_yptar")):
        payload = _mapping(joints[key], "plot_joint")
        vmax = _row_vmax(payload)
        shared_mesh = None
        for axis, population in zip(row, _POPULATIONS):
            mesh = _draw_joint(axis, payload, population, vmax)
            if shared_mesh is None and mesh is not None:
                shared_mesh = mesh
        if shared_mesh is not None:
            colorbar = figure.colorbar(shared_mesh, ax=list(row), pad=0.015)
            colorbar.set_label("Normalized bin fraction", fontsize=6)
            colorbar.ax.tick_params(labelsize=5.5)
    semantic_text = (
        "L: upstream low-HGCer observed reference before the downstream gate, 0 < NPE <= 2.  "
        "B: physical pion-control population, NPE > 2, with w0.  "
        "A: the same physical control population with w0*C.\n"
        "B-to-A is descriptive acceptance refinement only; it is not a uniform-improvement, yield, or promotion claim."
    )
    figure.text(0.015, 0.100, semantic_text, fontsize=6.8, va="top", wrap=True)
    figure.text(0.015, 0.060, "\n".join(_footer_lines(child)), fontsize=6.6, va="top", family="monospace")
    figure.subplots_adjust(left=0.055, right=0.965, top=0.92, bottom=0.165, hspace=0.58, wspace=0.38)
    pdf.savefig(figure)
    plt.close(figure)


def _build_manifest(
    validation: Mapping[str, object], input_name: str, renderer_commit: str, profile_commit: str,
) -> dict[str, object]:
    records: list[dict[str, object]] = [{"pdf_page_number": 1, "page_kind": "cover", "rendered": True, "empty": False}]
    page_number = 2
    for parent_value in _sequence(validation["parents"], "manifest_parents"):
        parent = _mapping(parent_value, "manifest_parent")
        for child_value in _sequence(parent["children"], "manifest_children"):
            child = _mapping(child_value, "manifest_child")
            empty = bool(_mapping(child["availability"], "manifest_availability")["completely_empty"])
            record: dict[str, object] = {
                "page_kind": "child",
                "setting_id": child["setting_id"],
                "canonical_t_index": child["canonical_t_index"],
                "phi_index": child["phi_index"],
                "phi_low": child["phi_low"],
                "phi_high": child["phi_high"],
                "rendered": not empty,
                "empty": empty,
            }
            if empty:
                record["pdf_page_number"] = None
            else:
                record["pdf_page_number"] = page_number
                page_number += 1
            records.append(record)
    return {
        "schema_version": MANIFEST_SCHEMA,
        "input_json_basename": input_name,
        "input_json_sha256": INPUT_SHA256,
        "artifact_fingerprint": ARTIFACT_FINGERPRINT,
        "validation_fingerprint": VALIDATION_FINGERPRINT,
        "scientific_source_commit": SCIENTIFIC_SOURCE_COMMIT,
        "fix5_upstream_renderer_source_commit": FIX5_RENDERER_SOURCE_COMMIT,
        "fix5_upstream_bundle_profile_commit": FIX5_BUNDLE_PROFILE_COMMIT,
        "e8_renderer_source_commit": renderer_commit,
        "e8_bundle_profile_commit": profile_commit,
        "output_pdf_basename": PDF_BASENAME,
        "output_manifest_basename": MANIFEST_BASENAME,
        "pdf_page_count": page_number - 1,
        "parent_count": 15,
        "canonical_child_count": 135,
        "rendered_child_count": 115,
        "empty_child_count": 20,
        "pages": records,
    }


def _render_pdf(path: Path, validation: Mapping[str, object], input_name: str, renderer_commit: str, profile_commit: str) -> None:
    metadata = {
        "Title": "KaonLT E.8 F.6.2 frozen figure library",
        "Author": "KaonLT",
        "Creator": "KaonLT E.8 figure-library renderer",
        "Producer": "Matplotlib",
        "CreationDate": _FIXED_PDF_DATE,
        "ModDate": _FIXED_PDF_DATE,
    }
    window = _mapping(validation["kaon_window"], "render_kaon_window")
    with PdfPages(path, metadata=metadata) as pdf:
        _cover_page(pdf, input_name, renderer_commit, profile_commit)
        for parent_value in _sequence(validation["parents"], "render_parents"):
            parent = _mapping(parent_value, "render_parent")
            for child_value in _sequence(parent["children"], "render_children"):
                child = _mapping(child_value, "render_child")
                if not bool(_mapping(child["availability"], "render_availability")["completely_empty"]):
                    _child_page(pdf, child, window)


def _write_json(path: Path, value: Mapping[str, object]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n")


def _temporary_output(destination: Path, suffix: str) -> Path:
    descriptor, name = tempfile.mkstemp(prefix=".e8_", suffix=suffix, dir=destination.parent)
    os.close(descriptor)
    return Path(name)


def _promote_without_overwrite(temporary: Path, destination: Path) -> None:
    """Publish a same-directory private file without replacing another output."""

    os.link(temporary, destination)
    temporary.unlink()


def _write_figure_library(
    validation: Mapping[str, object], *, input_name: str, renderer_commit: str, profile_commit: str,
    output_pdf: Path, output_manifest: Path,
) -> dict[str, object]:
    """Render an already validated fixture or production payload atomically."""

    temporary_pdf: Path | None = None
    temporary_manifest: Path | None = None
    promoted: list[Path] = []
    try:
        temporary_pdf = _temporary_output(output_pdf, ".pdf")
        temporary_manifest = _temporary_output(output_manifest, ".json")
        _render_pdf(temporary_pdf, validation, input_name, renderer_commit, profile_commit)
        manifest = _build_manifest(validation, input_name, renderer_commit, profile_commit)
        _write_json(temporary_manifest, manifest)
        _promote_without_overwrite(temporary_pdf, output_pdf)
        promoted.append(output_pdf)
        temporary_pdf = None
        _promote_without_overwrite(temporary_manifest, output_manifest)
        promoted.append(output_manifest)
        temporary_manifest = None
        return manifest
    except Exception:
        for path in (temporary_pdf, temporary_manifest, *promoted):
            if path is not None and path.exists():
                path.unlink()
        raise


def _commit(value: object, label: str) -> str:
    if not isinstance(value, str) or _COMMIT_RE.fullmatch(value) is None:
        raise FigureLibraryError("{}_invalid".format(label))
    return value


def _resolve_invocation(arguments: argparse.Namespace) -> tuple[Path, Path, Path, str, str]:
    input_path = Path(arguments.input_json).resolve(strict=True)
    output_pdf = Path(arguments.output_pdf).resolve()
    output_manifest = Path(arguments.output_manifest).resolve()
    if output_pdf.name != PDF_BASENAME or output_manifest.name != MANIFEST_BASENAME:
        raise FigureLibraryError("output_basename_invalid")
    if not output_pdf.parent.is_dir() or not output_manifest.parent.is_dir():
        raise FigureLibraryError("output_directory_invalid")
    if len({input_path, output_pdf, output_manifest}) != 3:
        raise FigureLibraryError("input_output_paths_collide")
    if output_pdf.exists() or output_manifest.exists():
        raise FigureLibraryError("output_path_already_exists")
    return (
        input_path,
        output_pdf,
        output_manifest,
        _commit(arguments.renderer_source_commit, "renderer_source_commit"),
        _commit(arguments.bundle_profile_commit, "bundle_profile_commit"),
    )


def _parse_input_bytes(input_bytes: bytes, *, expected_sha256: str = INPUT_SHA256) -> object:
    """Hash before parsing; the injectable expected value is test-only plumbing."""

    if hashlib.sha256(input_bytes).hexdigest() != expected_sha256:
        raise FigureLibraryError("input_sha256_mismatch")
    try:
        return json.loads(input_bytes.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FigureLibraryError("input_json_invalid") from exc


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-json", required=True, type=Path)
    parser.add_argument("--output-pdf", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--renderer-source-commit", required=True)
    parser.add_argument("--bundle-profile-commit", required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    try:
        input_path, output_pdf, output_manifest, renderer_commit, profile_commit = _resolve_invocation(arguments)
        input_bytes = input_path.read_bytes()
        artifact = _parse_input_bytes(input_bytes)
        validation = validate_artifact(artifact)
        _write_figure_library(
            validation,
            input_name=input_path.name,
            renderer_commit=renderer_commit,
            profile_commit=profile_commit,
            output_pdf=output_pdf,
            output_manifest=output_manifest,
        )
    except (FigureLibraryError, OSError, ValueError, TypeError) as exc:
        print("error: {}".format(exc), file=sys.stderr)
        return 1
    print("wrote {}".format(output_pdf))
    print("wrote {}".format(output_manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
