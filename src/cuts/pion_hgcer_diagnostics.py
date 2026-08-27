"""Detached Part-1 HGCer-versus-delta pion diagnostics.

This module is deliberately downstream of neither pion fitting nor yield
production.  It reads the same accepted experimental source events as the
production paths, bypasses only a *downstream* HGCer NPE threshold, and returns
detached diagnostic objects plus JSON-safe accounting.  Nothing here produces
a particle-subtraction weight or mutates a production histogram.
"""

from __future__ import annotations

from array import array
from copy import deepcopy
import hashlib
import json
import math
import os
from types import MappingProxyType

try:  # Pure-Python contract tests intentionally run without PyROOT.
    import ROOT
except ImportError:  # pragma: no cover - exercised in non-ROOT test runners
    ROOT = None

from canonical_binning import find_canonical_bin
from data_coordinates import analysis_event_coordinates, validate_kaon_data_coordinate_contract
from root_histogram_ownership import unique_root_object_name


DIAGNOSTIC_LABEL = "PION HGCer t-DELTA DIAGNOSTICS — PART 1 — NON-AUTHORITATIVE"
_SIDES = ("kaon", "pion")


class _PionHGCerDiagnosticBuildFailure(RuntimeError):
    """Tag a recoverable Part-1 build failure with its diagnostic stage."""

    def __init__(self, diagnostic_stage, original_exception):
        self.diagnostic_stage = str(diagnostic_stage)
        self.original_exception = original_exception
        super().__init__(str(original_exception))


def _iter_or_empty(value):
    """Return an iterable without applying boolean coercion to NumPy arrays."""
    return () if value is None else value


def _float_edges(edges):
    """Copy a possibly NumPy-backed edge sequence into plain Python floats."""
    return [float(edge) for edge in _iter_or_empty(edges)]


def canonical_t_delta_index(value, edges):
    """Use the shared [low, high) / final-edge-inclusive bin membership."""
    try:
        scalar = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(scalar):
        return None
    index = find_canonical_bin(scalar, tuple(_float_edges(edges)))
    return int(index) if index >= 0 else None


def resolve_pion_hgcer_delta_edges(config, proton_cleaning_result=None):
    """Resolve diagnostic delta edges without changing the proton fit binning."""
    resolved = dict(config or {})
    proton_edges = []
    if bool(resolved.get("reuse_proton_delta_edges", True)):
        proton_edges = _float_edges(
            (proton_cleaning_result or {}).get("delta_edges")
        )
    if len(proton_edges) >= 2 and all(
        proton_edges[index] < proton_edges[index + 1]
        for index in range(len(proton_edges) - 1)
    ):
        return proton_edges, "proton_cleaning_result.delta_edges"
    lower, upper = [float(value) for value in (resolved.get("delta_range") or (-10.0, 20.0))]
    bins = int(resolved.get("delta_bins", 10))
    if bins <= 0 or not lower < upper:
        raise ValueError("invalid pion HGCer diagnostic delta configuration")
    width = (upper - lower) / float(bins)
    return [lower + width * index for index in range(bins + 1)], "diagnostic_config.delta_range"


def _finite(value, default=None):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def _json_ready(value):
    if isinstance(value, MappingProxyType):
        value = dict(value)
    if isinstance(value, dict):
        return {str(key): _json_ready(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(child) for child in value]
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    try:
        if hasattr(value, "tolist"):
            return _json_ready(value.tolist())
        if hasattr(value, "item"):
            return _json_ready(value.item())
    except Exception:
        pass
    return value


def _config_fingerprint(config):
    encoded = json.dumps(_json_ready(config or {}), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _empty_metrics():
    return {
        "record_count": 0,
        "allcuts_count": 0,
        "nommcuts_count": 0,
        "signed_weight_sum": 0.0,
        "absolute_weight_support": 0.0,
        "sum_weight_squared": 0.0,
        "effective_entries": 0.0,
    }


def _finalize_metrics(metrics):
    result = dict(metrics or {})
    absolute_support = float(result.get("absolute_weight_support", 0.0) or 0.0)
    sum_weight_squared = float(result.get("sum_weight_squared", 0.0) or 0.0)
    result["effective_entries"] = (
        float((absolute_support * absolute_support) / sum_weight_squared)
        if sum_weight_squared > 0.0
        else 0.0
    )
    return result


def classify_pion_hgcer_support(kaon_metrics, pion_metrics, thresholds):
    """Classify a cell from non-cancelling support on both observed sides."""
    thresholds = dict(thresholds or {})
    supported_weight = float(thresholds.get("supported_absolute_weight", 0.0) or 0.0)
    marginal_weight = float(thresholds.get("marginal_absolute_weight", 0.0) or 0.0)
    supported_entries = float(thresholds.get("supported_effective_entries", 0.0) or 0.0)
    marginal_entries = float(thresholds.get("marginal_effective_entries", 0.0) or 0.0)

    def meets(metrics, weight_floor, entry_floor):
        metrics = _finalize_metrics(metrics)
        return (
            float(metrics.get("absolute_weight_support", 0.0) or 0.0) >= weight_floor
            and float(metrics.get("effective_entries", 0.0) or 0.0) >= entry_floor
        )

    if meets(kaon_metrics, supported_weight, supported_entries) and meets(pion_metrics, supported_weight, supported_entries):
        return "supported"
    if meets(kaon_metrics, marginal_weight, marginal_entries) and meets(pion_metrics, marginal_weight, marginal_entries):
        return "marginal"
    return "unsupported"


def _require_root():
    if ROOT is None:
        raise RuntimeError("PyROOT is required to build pion HGCer diagnostic histograms")


def _new_hist_1d(name, title, bins, lower, upper):
    _require_root()
    histogram = ROOT.TH1D(unique_root_object_name(name, scope="pion_hgcer", role="part1"), title, int(bins), float(lower), float(upper))
    histogram.SetDirectory(0)
    histogram.Sumw2()
    return histogram


def _new_hist_2d(name, title, x_bins, x_lower, x_upper, y_bins, y_lower, y_upper):
    _require_root()
    histogram = ROOT.TH2D(
        unique_root_object_name(name, scope="pion_hgcer", role="part1"),
        title,
        int(x_bins), float(x_lower), float(x_upper),
        int(y_bins), float(y_lower), float(y_upper),
    )
    histogram.SetDirectory(0)
    histogram.Sumw2()
    return histogram


def _new_variable_x_hist_1d(name, title, x_edges):
    _require_root()
    histogram = ROOT.TH1D(
        unique_root_object_name(name, scope="pion_hgcer", role="part1"),
        title,
        len(x_edges) - 1, array("d", [float(edge) for edge in x_edges]),
    )
    histogram.SetDirectory(0)
    histogram.Sumw2()
    return histogram


def _new_variable_x_hist_2d(name, title, x_edges, y_bins, y_lower, y_upper):
    _require_root()
    histogram = ROOT.TH2D(
        unique_root_object_name(name, scope="pion_hgcer", role="part1"),
        title,
        len(x_edges) - 1, array("d", [float(edge) for edge in x_edges]),
        int(y_bins), float(y_lower), float(y_upper),
    )
    histogram.SetDirectory(0)
    histogram.Sumw2()
    return histogram


def _new_delta_t_hist(name, title, delta_edges, t_edges):
    _require_root()
    histogram = ROOT.TH2D(
        unique_root_object_name(name, scope="pion_hgcer", role="part1"),
        title,
        len(delta_edges) - 1, array("d", [float(edge) for edge in delta_edges]),
        len(t_edges) - 1, array("d", [float(edge) for edge in t_edges]),
    )
    histogram.SetDirectory(0)
    histogram.Sumw2()
    return histogram


def _make_histograms(t_edges, delta_edges, config):
    lower_npe, upper_npe = [float(value) for value in config["hgcer_npe_range"]]
    lower_mm, upper_mm = [float(value) for value in config["mm_range"]]
    npe_bins = int(config["hgcer_npe_bins"])
    mm_bins = int(config["mm_bins"])
    histograms = {}
    for side in _SIDES:
        side_title = "proton-cleaned kaon" if side == "kaon" else "pion control"
        for weighting, title_suffix in (("weighted", "signed weighted yield"), ("absolute", "absolute support")):
            histograms["H_hgcer_{}_{}".format(side, weighting)] = _new_hist_1d(
                "H_hgcer_{}_{}".format(side, weighting),
                f"{side_title} HGCer NPE ({title_suffix});P_hgcer_npeSum;{title_suffix}",
                npe_bins, lower_npe, upper_npe,
            )
            histograms["H_delta_{}_{}".format(side, weighting)] = _new_variable_x_hist_1d(
                "H_delta_{}_{}".format(side, weighting),
                f"{side_title} SHMS #delta ({title_suffix});ssdelta [%];{title_suffix}",
                delta_edges,
            )
            histograms["H_mm_{}_{}".format(side, weighting)] = _new_hist_1d(
                "H_mm_{}_{}".format(side, weighting),
                f"{side_title} shifted MM ({title_suffix});shifted MM [GeV];{title_suffix}",
                mm_bins, lower_mm, upper_mm,
            )
            histograms["H_hgcer_vs_delta_{}_{}".format(side, weighting)] = _new_variable_x_hist_2d(
                "H_hgcer_vs_delta_{}_{}".format(side, weighting),
                "{} HGCer versus #delta;ssdelta [%];P_hgcer_npeSum".format(side_title),
                delta_edges,
                npe_bins, lower_npe, upper_npe,
            )
            histograms["H_hgcer_vs_t_{}_{}".format(side, weighting)] = _new_variable_x_hist_2d(
                "H_hgcer_vs_t_{}_{}".format(side, weighting),
                "{} HGCer versus canonical |t|;|t| [GeV^{2}];P_hgcer_npeSum".format(side_title),
                t_edges,
                npe_bins, lower_npe, upper_npe,
            )
            histograms["H_mm_vs_delta_{}_{}".format(side, weighting)] = _new_variable_x_hist_2d(
                "H_mm_vs_delta_{}_{}".format(side, weighting),
                "{} shifted MM versus #delta;ssdelta [%];shifted MM [GeV]".format(side_title),
                delta_edges,
                mm_bins, lower_mm, upper_mm,
            )
            histograms["H_mm_vs_hgcer_{}_{}".format(side, weighting)] = _new_hist_2d(
                "H_mm_vs_hgcer_{}_{}".format(side, weighting),
                "{} shifted MM versus HGCer;P_hgcer_npeSum;shifted MM [GeV]".format(side_title),
                npe_bins, lower_npe, upper_npe,
                mm_bins, lower_mm, upper_mm,
            )
        histograms["H_support_absolute_{}".format(side)] = _new_delta_t_hist(
            "H_support_absolute_{}".format(side),
            "{} absolute support;SHMS #delta [%];canonical |t| [GeV^{2}]".format(side_title),
            delta_edges, t_edges,
        )
        histograms["H_support_effective_{}".format(side)] = _new_delta_t_hist(
            "H_support_effective_{}".format(side),
            "{} effective entries;SHMS #delta [%];canonical |t| [GeV^{2}]".format(side_title),
            delta_edges, t_edges,
        )
    histograms["H_support_class"] = _new_delta_t_hist(
        "H_support_class",
        "Part-1 HGCer support class (0 unsupported, 1 marginal, 2 supported);SHMS #delta [%];canonical |t| [GeV^{2}]",
        delta_edges, t_edges,
    )
    return histograms


def _new_source_audit(side, source_label, source_spec, *, threshold, selection_state, coordinate_fingerprint):
    return {
        "side": str(side),
        "source_label": str(source_label),
        "tree_name": (source_spec or {}).get("tree_name"),
        "rf_state": "noRF" if side == "kaon" else "production_pion_source_tree",
        "coefficient": float((source_spec or {}).get("coefficient", 0.0) or 0.0),
        "source_selection_state": str(selection_state),
        "production_hgcer_threshold": float(threshold),
        "downstream_hgcer_threshold_bypassed": True,
        "coordinate_fingerprint": str(coordinate_fingerprint),
        "tree_entries_seen": 0,
        "accepted_pre_hgcer_records": 0,
        "records_inside_canonical_t": 0,
        "records_inside_delta": 0,
        "nonfinite_hgcer": 0,
        "nonfinite_coordinate": 0,
        "missing_frozen_proton_factor": 0,
        "proton_factor_records": 0,
        "proton_factor_sum": 0.0,
        "proton_factor_mean": None,
        "hgc_below_threshold_records": 0,
        "hgc_at_or_above_threshold_records": 0,
        "raw_to_analysis_t_migration": {},
        "coordinate_closure": {
            "checked_records": 0,
            "maximum_abs_mm_residual": 0.0,
            "maximum_abs_t_residual": 0.0,
        },
        **_empty_metrics(),
    }


def _record_metric(metrics, record):
    weight = float(record["diagnostic_weight"])
    metrics["record_count"] += 1
    metrics["allcuts_count"] += int(bool(record["allcuts"]))
    metrics["nommcuts_count"] += int(bool(record["nommcuts"]))
    metrics["signed_weight_sum"] += weight
    metrics["absolute_weight_support"] += abs(weight)
    metrics["sum_weight_squared"] += weight * weight


def _record_t_migration(audit, raw_t, analysis_t, t_edges):
    raw_index = canonical_t_delta_index(raw_t, t_edges)
    analysis_index = canonical_t_delta_index(analysis_t, t_edges)
    key = "{}->{}".format("out" if raw_index is None else raw_index, "out" if analysis_index is None else analysis_index)
    migration = audit["raw_to_analysis_t_migration"]
    migration[key] = int(migration.get(key, 0)) + 1


def _freeze_records(records):
    return tuple(MappingProxyType(dict(record)) for record in records)


def build_pion_hgcer_tdelta_diagnostic(
    *,
    kaon_source_bundle,
    pion_tree_bundle,
    proton_cleaning_result,
    proton_coordinate_fingerprint=None,
    pion_control_cache,
    coordinate_contract,
    t_edges,
    config,
    hole_contains,
    evaluate_pion_event,
    mm_min,
    mm_max,
):
    """Build the isolated pre-HGCer diagnostic population and ROOT products.

    Kaon inputs come from the immutable prepared noRF source bundle and use the
    committed *pre-RF* ``cleaned_factor``.  Pion inputs re-read the same source
    trees used by the authoritative cache but deliberately omit its downstream
    HGCer PID gate.  Neither record collection is returned to production.
    """
    _require_root()
    coordinate = validate_kaon_data_coordinate_contract(coordinate_contract)
    coordinate_fingerprint = str(coordinate["coordinate_fingerprint"])
    t_edges = _float_edges(t_edges)
    if len(t_edges) < 2 or any(t_edges[index] >= t_edges[index + 1] for index in range(len(t_edges) - 1)):
        raise ValueError("pion HGCer diagnostic requires resolved increasing canonical t edges")
    config = deepcopy(config or {})
    delta_edges, delta_source = resolve_pion_hgcer_delta_edges(config, proton_cleaning_result)
    threshold = float(config.get("production_hgcer_threshold", 2.0))
    cache_fingerprint = str((pion_control_cache or {}).get("coordinate_fingerprint") or "")
    if cache_fingerprint and cache_fingerprint != coordinate_fingerprint:
        raise RuntimeError("pion_hgcer_diagnostic_coordinate_mismatch_with_authoritative_pion_cache")
    application_coordinate = str(
        proton_coordinate_fingerprint
        or (proton_cleaning_result or {}).get("coordinate_fingerprint")
        or ""
    )
    if application_coordinate and application_coordinate != coordinate_fingerprint:
        raise RuntimeError("pion_hgcer_diagnostic_coordinate_mismatch_with_proton_result")

    try:
        histograms = _make_histograms(t_edges, delta_edges, config)
    except Exception as exc:
        raise _PionHGCerDiagnosticBuildFailure(
            "histogram_construction", exc
        ) from exc
    records = {side: [] for side in _SIDES}
    source_audit = {side: {} for side in _SIDES}
    cells = {
        (t_index, delta_index): {
            "t_index": int(t_index),
            "t_low": float(t_edges[t_index]),
            "t_high": float(t_edges[t_index + 1]),
            "delta_index": int(delta_index),
            "delta_low": float(delta_edges[delta_index]),
            "delta_high": float(delta_edges[delta_index + 1]),
            "kaon": _empty_metrics(),
            "pion": _empty_metrics(),
        }
        for t_index in range(len(t_edges) - 1)
        for delta_index in range(len(delta_edges) - 1)
    }

    def accept_record(side, source_label, source_spec, entry_index, evt, *, allcuts, nommcuts, proton_factor=1.0, factor_payload=None, selection_state):
        audit = source_audit[side].setdefault(
            str(source_label),
            _new_source_audit(
                side, source_label, source_spec,
                threshold=threshold,
                selection_state=selection_state,
                coordinate_fingerprint=coordinate_fingerprint,
            ),
        )
        if not (bool(allcuts) or bool(nommcuts)):
            return
        audit["accepted_pre_hgcer_records"] += 1
        npe = _finite(getattr(evt, "P_hgcer_npeSum", None))
        if npe is None:
            audit["nonfinite_hgcer"] += 1
            return
        if npe < threshold:
            audit["hgc_below_threshold_records"] += 1
        else:
            audit["hgc_at_or_above_threshold_records"] += 1
        coordinates = analysis_event_coordinates(evt, coordinate)
        raw_mm = _finite(coordinates.get("raw_mm"))
        raw_t = _finite(coordinates.get("raw_t"))
        analysis_mm = _finite(coordinates.get("analysis_mm"))
        analysis_t = _finite(coordinates.get("analysis_t"))
        if None in (raw_mm, raw_t, analysis_mm, analysis_t):
            audit["nonfinite_coordinate"] += 1
            return
        expected_mm = float(raw_mm) + float(coordinate["mm_shift"])
        expected_t = float(raw_t) + float(coordinate["t_shift"])
        coordinate_closure = audit["coordinate_closure"]
        coordinate_closure["checked_records"] += 1
        coordinate_closure["maximum_abs_mm_residual"] = max(
            float(coordinate_closure["maximum_abs_mm_residual"]),
            abs(float(analysis_mm) - expected_mm),
        )
        coordinate_closure["maximum_abs_t_residual"] = max(
            float(coordinate_closure["maximum_abs_t_residual"]),
            abs(float(analysis_t) - expected_t),
        )
        _record_t_migration(audit, raw_t, analysis_t, t_edges)
        t_index = canonical_t_delta_index(analysis_t, t_edges)
        if t_index is None:
            return
        audit["records_inside_canonical_t"] += 1
        delta_value = _finite(getattr(evt, "ssdelta", None))
        delta_index = canonical_t_delta_index(delta_value, delta_edges)
        if delta_index is None:
            return
        audit["records_inside_delta"] += 1
        coefficient = float((source_spec or {}).get("coefficient", 0.0) or 0.0)
        factor = _finite(proton_factor)
        if factor is None:
            audit["missing_frozen_proton_factor"] += 1
            return
        if side == "kaon":
            audit["proton_factor_records"] += 1
            audit["proton_factor_sum"] += float(factor)
        diagnostic_weight = coefficient * factor
        record = {
            "side": str(side),
            "source_label": str(source_label),
            "entry_index": int(entry_index),
            "coefficient": coefficient,
            "proton_cleaning_factor": float(factor) if side == "kaon" else None,
            "diagnostic_weight": float(diagnostic_weight),
            "raw_MM": float(raw_mm),
            "raw_t": float(raw_t),
            "analysis_MM": float(analysis_mm),
            "analysis_t": float(analysis_t),
            "coordinate_fingerprint": coordinate_fingerprint,
            "canonical_t_index": int(t_index),
            "ssdelta": float(delta_value),
            "delta_index": int(delta_index),
            "P_hgcer_npeSum": float(npe),
            "allcuts": bool(allcuts),
            "nommcuts": bool(nommcuts),
            "Q2": _finite(getattr(evt, "Q2", None)),
            "W": _finite(getattr(evt, "W", None)),
            "epsilon": _finite(getattr(evt, "epsilon", None)),
            "phi": _finite(getattr(evt, "ph_q", None)),
            "frozen_proton_payload_present": bool(factor_payload is not None) if side == "kaon" else False,
            "rf_applied_to_diagnostic": False,
        }
        records[side].append(record)
        _record_metric(audit, record)
        _record_metric(cells[(t_index, delta_index)][side], record)
        absolute_weight = abs(diagnostic_weight)
        for weighting, fill_weight in (("weighted", diagnostic_weight), ("absolute", absolute_weight)):
            histograms["H_hgcer_{}_{}".format(side, weighting)].Fill(npe, fill_weight)
            histograms["H_delta_{}_{}".format(side, weighting)].Fill(delta_value, fill_weight)
            histograms["H_mm_{}_{}".format(side, weighting)].Fill(analysis_mm, fill_weight)
            histograms["H_hgcer_vs_delta_{}_{}".format(side, weighting)].Fill(delta_value, npe, fill_weight)
            histograms["H_hgcer_vs_t_{}_{}".format(side, weighting)].Fill(analysis_t, npe, fill_weight)
            histograms["H_mm_vs_delta_{}_{}".format(side, weighting)].Fill(delta_value, analysis_mm, fill_weight)
            histograms["H_mm_vs_hgcer_{}_{}".format(side, weighting)].Fill(npe, analysis_mm, fill_weight)

    # The kaon side comes from the source records already admitted by the
    # proton preparer.  Those records share the exact accepted/hole/coordinate
    # state and have never received the kaon HGCer PID gate.
    active_kaon_sources = (kaon_source_bundle or {}).get("sources") or {}
    prepared_kaon_sources = (kaon_source_bundle or {}).get("prepared_sources") or {}
    frozen_lookup = (proton_cleaning_result or {}).get("_prepared_event_weight_lookup") or {}
    for source_label, prepared_spec in prepared_kaon_sources.items():
        active_spec = active_kaon_sources.get(source_label) or prepared_spec or {}
        tree = active_spec.get("tree")
        entry_payloads = (prepared_spec or {}).get("entries") or {}
        if tree is None:
            source_audit["kaon"][str(source_label)] = _new_source_audit(
                "kaon", source_label, active_spec,
                threshold=threshold,
                selection_state="prepared_noRF_pre_particle_subtraction",
                coordinate_fingerprint=coordinate_fingerprint,
            )
            continue
        kaon_audit = source_audit["kaon"].setdefault(
            str(source_label),
            _new_source_audit(
                "kaon", source_label, active_spec,
                threshold=threshold,
                selection_state="prepared_noRF_pre_particle_subtraction",
                coordinate_fingerprint=coordinate_fingerprint,
            ),
        )
        try:
            kaon_audit["tree_entries_seen"] = int(tree.GetEntries())
        except Exception:
            try:
                kaon_audit["tree_entries_seen"] = int(len(tree))
            except Exception:
                kaon_audit["tree_entries_seen"] = 0
        for entry_index, evt in enumerate(tree):
            prepared = entry_payloads.get(int(entry_index))
            if prepared is None:
                continue
            frozen_payload = frozen_lookup.get("{}:{}".format(source_label, int(entry_index)))
            factor = (frozen_payload or {}).get("cleaned_factor")
            accept_record(
                "kaon", source_label, active_spec, entry_index, evt,
                allcuts=prepared.get("allcuts", False),
                nommcuts=prepared.get("nommcuts", False),
                proton_factor=factor,
                factor_payload=frozen_payload,
                selection_state="prepared_noRF_pre_particle_subtraction",
            )

    # The pion side deliberately reuses the exact source trees opened for the
    # authoritative control cache.  Its production-only HGCer threshold is
    # omitted here, while every acceptance/hole/MM/t decision remains intact.
    pion_source_specs = (
        ("prompt", (pion_tree_bundle or {}).get("prompt_tree"), float(((kaon_source_bundle or {}).get("sources") or {}).get("prompt", {}).get("coefficient", 0.0) or 0.0), "prompt_tree"),
        ("rand", (pion_tree_bundle or {}).get("rand_tree"), float(((kaon_source_bundle or {}).get("sources") or {}).get("rand", {}).get("coefficient", 0.0) or 0.0), "rand_tree"),
        ("dummy", (pion_tree_bundle or {}).get("dummy_prompt_tree"), float(((kaon_source_bundle or {}).get("sources") or {}).get("dummy_prompt", {}).get("coefficient", 0.0) or 0.0), "dummy_prompt_tree"),
        ("dummy_rand", (pion_tree_bundle or {}).get("dummy_rand_tree"), float(((kaon_source_bundle or {}).get("sources") or {}).get("dummy_rand", {}).get("coefficient", 0.0) or 0.0), "dummy_rand_tree"),
    )
    for source_label, tree, coefficient, tree_role in pion_source_specs:
        source_spec = {
            "tree": tree,
            "tree_name": (pion_tree_bundle or {}).get(
                "prompt_tree_name" if source_label in ("prompt", "dummy") else "rand_tree_name"
            ),
            "coefficient": coefficient,
        }
        if tree is None:
            source_audit["pion"][str(source_label)] = _new_source_audit(
                "pion", source_label, source_spec,
                threshold=threshold,
                selection_state="same_source_tree_downstream_hgcer_bypassed:{}".format(tree_role),
                coordinate_fingerprint=coordinate_fingerprint,
            )
            continue
        pion_audit = source_audit["pion"].setdefault(
            str(source_label),
            _new_source_audit(
                "pion", source_label, source_spec,
                threshold=threshold,
                selection_state="same_source_tree_downstream_hgcer_bypassed:{}".format(tree_role),
                coordinate_fingerprint=coordinate_fingerprint,
            ),
        )
        try:
            pion_audit["tree_entries_seen"] = int(tree.GetEntries())
        except Exception:
            try:
                pion_audit["tree_entries_seen"] = int(len(tree))
            except Exception:
                pion_audit["tree_entries_seen"] = 0
        for entry_index, evt in enumerate(tree):
            base_allcuts, base_nommcuts, _ = evaluate_pion_event(evt, mm_min, mm_max, mm_offset=0.0)
            hole_rejected = bool(hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer)) if hole_contains is not None else False
            accept_record(
                "pion", source_label, source_spec, entry_index, evt,
                allcuts=bool(base_allcuts and not hole_rejected),
                nommcuts=bool(base_nommcuts and not hole_rejected),
                proton_factor=1.0,
                selection_state="same_source_tree_downstream_hgcer_bypassed:{}".format(tree_role),
            )

    serialized_cells = []
    for key, cell in cells.items():
        cell["kaon"] = _finalize_metrics(cell["kaon"])
        cell["pion"] = _finalize_metrics(cell["pion"])
        cell["support_class"] = classify_pion_hgcer_support(
            cell["kaon"], cell["pion"], config.get("support_thresholds")
        )
        support_number = {"unsupported": 0.0, "marginal": 1.0, "supported": 2.0}[cell["support_class"]]
        x_bin = int(cell["delta_index"]) + 1
        y_bin = int(cell["t_index"]) + 1
        for side in _SIDES:
            histograms["H_support_absolute_{}".format(side)].SetBinContent(
                x_bin, y_bin, float(cell[side]["absolute_weight_support"])
            )
            histograms["H_support_effective_{}".format(side)].SetBinContent(
                x_bin, y_bin, float(cell[side]["effective_entries"])
            )
        histograms["H_support_class"].SetBinContent(x_bin, y_bin, support_number)
        serialized_cells.append(cell)
    for side in _SIDES:
        for audit in source_audit[side].values():
            finalized = _finalize_metrics(audit)
            audit.update(finalized)
            audit["tree_hgcer_truncation_evidence"] = (
                "pre_hgcer_events_observed"
                if int(audit.get("hgc_below_threshold_records", 0)) > 0
                else "no_below_threshold_events_observed_tree_may_be_truncated_or_statistically_empty"
            )
            if int(audit.get("proton_factor_records", 0)) > 0:
                audit["proton_factor_mean"] = float(
                    float(audit["proton_factor_sum"])
                    / float(audit["proton_factor_records"])
                )

    frozen_records = {side: _freeze_records(section) for side, section in records.items()}
    return {
        "status": "available",
        "diagnostic_label": DIAGNOSTIC_LABEL,
        "non_authoritative": True,
        "production_side_effect_free": True,
        "production_hgcer_pid_unchanged": True,
        "proton_factor_scope": "kaon_cleaned_factor_pre_rf_only",
        "rf_restoration_applied": False,
        "coordinate_contract": dict(coordinate),
        "coordinate_fingerprint": coordinate_fingerprint,
        "config": config,
        "config_fingerprint": _config_fingerprint(config),
        "t_edges": t_edges,
        "delta_edges": delta_edges,
        "delta_edge_source": delta_source,
        "source_audit": source_audit,
        "records": frozen_records,
        "cells": tuple(MappingProxyType(dict(cell)) for cell in serialized_cells),
        "histograms": histograms,
    }


def serialize_pion_hgcer_tdelta_diagnostic(payload, include_records=True):
    """Return a JSON-safe Part-1 diagnostic payload without ROOT objects."""
    payload = payload or {}
    records = payload.get("records")
    records = {} if records is None else records
    result = {
        "status": payload.get("status", "unavailable"),
        "reason": payload.get("reason"),
        "exception_type": payload.get("exception_type"),
        "exception_message": payload.get("exception_message"),
        "diagnostic_stage": payload.get("diagnostic_stage"),
        "diagnostic_label": payload.get("diagnostic_label", DIAGNOSTIC_LABEL),
        "non_authoritative": True,
        "production_side_effect_free": True,
        "production_hgcer_pid_unchanged": bool(payload.get("production_hgcer_pid_unchanged", True)),
        "proton_factor_scope": payload.get("proton_factor_scope"),
        "rf_restoration_applied": bool(payload.get("rf_restoration_applied", False)),
        "coordinate_contract": payload.get("coordinate_contract") or {},
        "coordinate_fingerprint": payload.get("coordinate_fingerprint"),
        "config": payload.get("config") or {},
        "config_fingerprint": payload.get("config_fingerprint"),
        "t_edges": _iter_or_empty(payload.get("t_edges")),
        "delta_edges": _iter_or_empty(payload.get("delta_edges")),
        "delta_edge_source": payload.get("delta_edge_source"),
        "source_audit": payload.get("source_audit") or {},
        "cells": _iter_or_empty(payload.get("cells")),
        "page_manifest": _iter_or_empty(payload.get("page_manifest")),
        "record_counts": {
            side: len(_iter_or_empty(records.get(side))) for side in _SIDES
        },
    }
    if include_records:
        result["records"] = {
            side: list(_iter_or_empty(records.get(side))) for side in _SIDES
        }
    return _json_ready(result)


def write_pion_hgcer_tdelta_json(path, payload):
    serialized = serialize_pion_hgcer_tdelta_diagnostic(payload)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(serialized, handle, sort_keys=True, indent=2, allow_nan=False)
    return path


def _normalized_clone(histogram, name):
    clone = histogram.Clone(unique_root_object_name(name, scope="pion_hgcer", role="plot"))
    clone.SetDirectory(0)
    integral = float(clone.Integral())
    if integral > 0.0:
        clone.Scale(1.0 / integral)
    return clone


def _style_overlay(kaon_hist, pion_hist, *, shape=False):
    kaon_hist.SetLineColor(ROOT.kBlue + 1)
    kaon_hist.SetLineWidth(2)
    pion_hist.SetLineColor(ROOT.kRed + 1)
    pion_hist.SetLineWidth(2)
    maximum = max(float(kaon_hist.GetMaximum()), float(pion_hist.GetMaximum()), 1.0e-12)
    kaon_hist.SetMaximum(1.20 * maximum)
    kaon_hist.SetMinimum(0.0 if shape else min(0.0, 1.10 * min(float(kaon_hist.GetMinimum()), float(pion_hist.GetMinimum()))))
    kaon_hist.Draw("hist")
    pion_hist.Draw("hist same")
    legend = ROOT.TLegend(0.57, 0.74, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(kaon_hist, "proton-cleaned kaon", "l")
    legend.AddEntry(pion_hist, "pion control", "l")
    legend.Draw()
    return legend


def _page_text_canvas(title, lines):
    canvas = ROOT.TCanvas(unique_root_object_name("C_pion_hgcer_status", scope="pion_hgcer", role="page"), title, 1100, 750)
    box = ROOT.TPaveText(0.06, 0.06, 0.94, 0.94, "NDC")
    box.SetFillStyle(0)
    box.SetBorderSize(0)
    box.SetTextAlign(12)
    box.AddText(title)
    for line in lines:
        box.AddText(str(line))
    box.Draw()
    canvas._pion_hgcer_box = box
    return canvas


def _projection(records, side, *, t_index=None, delta_index=None, quantity="npe", weight_mode="absolute", config=None, name="projection"):
    config = config or {}
    if quantity == "npe":
        lower, upper = [float(value) for value in config["hgcer_npe_range"]]
        histogram = _new_hist_1d(name, "HGCer NPE;P_hgcer_npeSum;{}".format("signed weighted yield" if weight_mode == "weighted" else "absolute support"), int(config["hgcer_npe_bins"]), lower, upper)
        field = "P_hgcer_npeSum"
    else:
        lower, upper = [float(value) for value in config["mm_range"]]
        histogram = _new_hist_1d(name, "Shifted MM;shifted MM [GeV];{}".format("signed weighted yield" if weight_mode == "weighted" else "absolute support"), int(config["mm_bins"]), lower, upper)
        field = "analysis_MM"
    for record in _iter_or_empty(records.get(side)):
        if t_index is not None and int(record["canonical_t_index"]) != int(t_index):
            continue
        if delta_index is not None and int(record["delta_index"]) != int(delta_index):
            continue
        weight = float(record["diagnostic_weight"])
        histogram.Fill(float(record[field]), weight if weight_mode == "weighted" else abs(weight))
    return histogram


def _npe_vs_delta(records, side, t_index, delta_edges, config, name):
    lower, upper = [float(value) for value in config["hgcer_npe_range"]]
    histogram = _new_variable_x_hist_2d(
        name, "HGCer NPE versus #delta;ssdelta [%];P_hgcer_npeSum",
        delta_edges,
        int(config["hgcer_npe_bins"]), lower, upper,
    )
    for record in _iter_or_empty(records.get(side)):
        if int(record["canonical_t_index"]) == int(t_index):
            histogram.Fill(float(record["ssdelta"]), float(record["P_hgcer_npeSum"]), abs(float(record["diagnostic_weight"])))
    return histogram


def _should_emit_cell(status, mode):
    mode = str(mode or "none")
    return mode == "all" or (mode == "supported" and status == "supported") or (mode == "supported_marginal" and status in {"supported", "marginal"})


def expected_pion_hgcer_page_manifest(payload):
    """Return the deterministic Part-1 page IDs without requiring PyROOT."""
    prefix = "pion_hgcer.part1"
    payload = payload or {}
    if str(payload.get("status") or "unavailable") != "available":
        return [{"page_id": "{}.status".format(prefix), "scope": "setting-wide", "authoritative": False}]
    t_edges = _iter_or_empty(payload.get("t_edges"))
    config = payload.get("config") or {}
    manifest = [
        {"page_id": "{}.audit".format(prefix), "scope": "setting-wide", "authoritative": False},
        {"page_id": "{}.setting_wide".format(prefix), "scope": "setting-wide", "authoritative": False},
        {"page_id": "{}.support".format(prefix), "scope": "setting-wide", "authoritative": False},
    ]
    for t_index in range(max(0, len(t_edges) - 1)):
        manifest.append({"page_id": "{}.t{}".format(prefix, t_index + 1), "scope": "t{}".format(t_index + 1), "authoritative": False})
    for cell in _iter_or_empty(payload.get("cells")):
        status = str(cell.get("support_class") or "unsupported")
        t_index = int(cell["t_index"])
        delta_index = int(cell["delta_index"])
        scope = "t{} delta{}".format(t_index + 1, delta_index + 1)
        if _should_emit_cell(status, config.get("emit_cell_pages")):
            page_id = "{}.t{}.delta{}".format(prefix, t_index + 1, delta_index + 1)
        elif bool(config.get("emit_status_pages", True)):
            page_id = "{}.t{}.delta{}.status".format(prefix, t_index + 1, delta_index + 1)
        else:
            continue
        manifest.append({"page_id": page_id, "scope": scope, "authoritative": False})
    return manifest


def render_pion_hgcer_tdelta_pages(pdf_name, payload, *, title_prefix="", page_manifest=None, close_pdf=False):
    """Append the Part-1 pages to an already-open Kaon ``rand_sub`` PDF."""
    if ROOT is None or not isinstance(payload, dict):
        return []
    page_manifest = page_manifest if page_manifest is not None else []
    prefix = "pion_hgcer.part1"
    if str(payload.get("status") or "unavailable") != "available":
        status_page = _page_text_canvas(
            DIAGNOSTIC_LABEL,
            (
                "{}".format(title_prefix),
                "status: unavailable",
                "reason: {}".format(payload.get("reason") or "not built"),
                "Production pion/proton selections and outputs were not changed.",
            ),
        )
        status_page.Print(str(pdf_name) + (")" if bool(close_pdf) else ""))
        entry = {"page_id": "{}.status".format(prefix), "scope": "setting-wide", "authoritative": False}
        page_manifest.append(entry)
        payload["page_manifest"] = [entry]
        return [entry]
    config = payload.get("config") or {}
    records = payload.get("records")
    records = {} if records is None else records
    histograms = payload.get("histograms") or {}
    t_edges = _iter_or_empty(payload.get("t_edges"))
    delta_edges = _iter_or_empty(payload.get("delta_edges"))
    pages = []

    audit_lines = [
        DIAGNOSTIC_LABEL,
        "{}".format(title_prefix),
        "coordinate fingerprint: {}".format(str(payload.get("coordinate_fingerprint") or "missing")[-16:]),
        "delta edges from {}: {}".format(payload.get("delta_edge_source"), delta_edges),
        "support thresholds: {}".format((payload.get("config") or {}).get("support_thresholds") or {}),
        "No production HGCer PID, pion fit, weight, or yield was changed.",
    ]
    for side in _SIDES:
        for label, audit in sorted((payload.get("source_audit") or {}).get(side, {}).items()):
            audit_lines.append(
                "{} {}: selected={} below NPE>{:.1f}={} factor-missing={} [{}]".format(
                    side, label, audit.get("accepted_pre_hgcer_records", 0),
                    float(audit.get("production_hgcer_threshold", 2.0)),
                    audit.get("hgc_below_threshold_records", 0),
                    audit.get("missing_frozen_proton_factor", 0),
                    audit.get("tree_hgcer_truncation_evidence", "unknown"),
                )
            )
            audit_lines.append(
                "  t migration={} coordinate max residuals=({:.2e}, {:.2e}) factor mean={}".format(
                    audit.get("raw_to_analysis_t_migration", {}),
                    float((audit.get("coordinate_closure") or {}).get("maximum_abs_mm_residual", 0.0)),
                    float((audit.get("coordinate_closure") or {}).get("maximum_abs_t_residual", 0.0)),
                    audit.get("proton_factor_mean"),
                )
            )
    pages.append(("{}.audit".format(prefix), _page_text_canvas(DIAGNOSTIC_LABEL, audit_lines), "setting-wide"))

    wide = ROOT.TCanvas(unique_root_object_name("C_pion_hgcer_wide", scope="pion_hgcer", role="page"), DIAGNOSTIC_LABEL, 1200, 850)
    wide.Divide(2, 2)
    wide.cd(1)
    _style_overlay(histograms["H_hgcer_kaon_weighted"], histograms["H_hgcer_pion_weighted"])
    wide.cd(2)
    kaon_shape = _normalized_clone(histograms["H_hgcer_kaon_absolute"], "H_hgcer_kaon_shape")
    pion_shape = _normalized_clone(histograms["H_hgcer_pion_absolute"], "H_hgcer_pion_shape")
    _style_overlay(kaon_shape, pion_shape, shape=True)
    wide.cd(3)
    histograms["H_hgcer_vs_delta_kaon_absolute"].Draw("colz")
    wide.cd(4)
    histograms["H_hgcer_vs_delta_pion_absolute"].Draw("colz")
    wide._pion_hgcer_shapes = (kaon_shape, pion_shape)
    pages.append(("{}.setting_wide".format(prefix), wide, "setting-wide"))

    support = ROOT.TCanvas(unique_root_object_name("C_pion_hgcer_support", scope="pion_hgcer", role="page"), DIAGNOSTIC_LABEL, 1200, 850)
    support.Divide(2, 2)
    for pad, key in enumerate(("H_support_absolute_kaon", "H_support_absolute_pion", "H_support_effective_kaon", "H_support_class"), start=1):
        support.cd(pad)
        histograms[key].Draw("colz text")
    pages.append(("{}.support".format(prefix), support, "setting-wide"))

    for t_index in range(max(0, len(t_edges) - 1)):
        page = ROOT.TCanvas(unique_root_object_name("C_pion_hgcer_t{}".format(t_index + 1), scope="pion_hgcer", role="page"), DIAGNOSTIC_LABEL, 1350, 850)
        page.Divide(3, 2)
        kaon_delta = _npe_vs_delta(records, "kaon", t_index, delta_edges, config, "H_hgcer_delta_kaon_t{}".format(t_index + 1))
        pion_delta = _npe_vs_delta(records, "pion", t_index, delta_edges, config, "H_hgcer_delta_pion_t{}".format(t_index + 1))
        page.cd(1); kaon_delta.Draw("colz")
        page.cd(2); pion_delta.Draw("colz")
        kaon_npe_weighted = _projection(records, "kaon", t_index=t_index, quantity="npe", weight_mode="weighted", config=config, name="H_hgcer_kaon_weighted_t{}".format(t_index + 1))
        pion_npe_weighted = _projection(records, "pion", t_index=t_index, quantity="npe", weight_mode="weighted", config=config, name="H_hgcer_pion_weighted_t{}".format(t_index + 1))
        page.cd(3); _style_overlay(kaon_npe_weighted, pion_npe_weighted)
        kaon_npe = _projection(records, "kaon", t_index=t_index, quantity="npe", config=config, name="H_hgcer_kaon_t{}".format(t_index + 1))
        pion_npe = _projection(records, "pion", t_index=t_index, quantity="npe", config=config, name="H_hgcer_pion_t{}".format(t_index + 1))
        page.cd(4); _style_overlay(kaon_npe, pion_npe, shape=True)
        kaon_mm = _projection(records, "kaon", t_index=t_index, quantity="mm", weight_mode="weighted", config=config, name="H_mm_kaon_t{}".format(t_index + 1))
        pion_mm = _projection(records, "pion", t_index=t_index, quantity="mm", weight_mode="weighted", config=config, name="H_mm_pion_t{}".format(t_index + 1))
        page.cd(5); _style_overlay(kaon_mm, pion_mm)
        page.cd(6)
        t_info = ROOT.TPaveText(0.08, 0.15, 0.92, 0.85, "NDC")
        t_info.SetFillStyle(0); t_info.SetBorderSize(0); t_info.SetTextAlign(12)
        t_info.AddText("canonical t{} [{:.4f}, {:.4f}]".format(t_index + 1, float(t_edges[t_index]), float(t_edges[t_index + 1])))
        t_info.AddText("Pre-HGCer diagnostic only; no #delta merge or P_{#pi} fit.")
        t_info.Draw()
        page._pion_hgcer_local = (kaon_delta, pion_delta, kaon_npe_weighted, pion_npe_weighted, kaon_npe, pion_npe, kaon_mm, pion_mm, t_info)
        pages.append(("{}.t{}".format(prefix, t_index + 1), page, "t{}".format(t_index + 1)))

    for cell in _iter_or_empty(payload.get("cells")):
        status = str(cell.get("support_class") or "unsupported")
        if not _should_emit_cell(status, config.get("emit_cell_pages")):
            if bool(config.get("emit_status_pages", True)):
                t_index = int(cell["t_index"])
                delta_index = int(cell["delta_index"])
                status_canvas = _page_text_canvas(
                    "{} — compact cell status".format(DIAGNOSTIC_LABEL),
                    (
                        "t{} [{:.4f}, {:.4f}]  delta{} [{:.2f}, {:.2f}]".format(
                            t_index + 1, cell["t_low"], cell["t_high"],
                            delta_index + 1, cell["delta_low"], cell["delta_high"],
                        ),
                        "support class: {}".format(status),
                        "kaon abs={:.3e}, Neff={:.2f}".format(
                            float((cell.get("kaon") or {}).get("absolute_weight_support", 0.0)),
                            float((cell.get("kaon") or {}).get("effective_entries", 0.0)),
                        ),
                        "pion abs={:.3e}, Neff={:.2f}".format(
                            float((cell.get("pion") or {}).get("absolute_weight_support", 0.0)),
                            float((cell.get("pion") or {}).get("effective_entries", 0.0)),
                        ),
                        "No response estimate, subtraction, or bin merge was performed.",
                    ),
                )
                pages.append((
                    "{}.t{}.delta{}.status".format(prefix, t_index + 1, delta_index + 1),
                    status_canvas,
                    "t{} delta{}".format(t_index + 1, delta_index + 1),
                ))
            continue
        t_index = int(cell["t_index"])
        delta_index = int(cell["delta_index"])
        cell_page = ROOT.TCanvas(unique_root_object_name("C_pion_hgcer_t{}_d{}".format(t_index + 1, delta_index + 1), scope="pion_hgcer", role="page"), DIAGNOSTIC_LABEL, 1100, 750)
        cell_page.Divide(2, 2)
        kaon_npe = _projection(records, "kaon", t_index=t_index, delta_index=delta_index, quantity="npe", config=config, name="H_hgcer_kaon_t{}_d{}".format(t_index + 1, delta_index + 1))
        pion_npe = _projection(records, "pion", t_index=t_index, delta_index=delta_index, quantity="npe", config=config, name="H_hgcer_pion_t{}_d{}".format(t_index + 1, delta_index + 1))
        cell_page.cd(1); _style_overlay(kaon_npe, pion_npe, shape=True)
        kaon_mm = _projection(records, "kaon", t_index=t_index, delta_index=delta_index, quantity="mm", config=config, name="H_mm_kaon_t{}_d{}".format(t_index + 1, delta_index + 1))
        pion_mm = _projection(records, "pion", t_index=t_index, delta_index=delta_index, quantity="mm", config=config, name="H_mm_pion_t{}_d{}".format(t_index + 1, delta_index + 1))
        cell_page.cd(2); _style_overlay(kaon_mm, pion_mm, shape=True)
        cell_page.cd(3)
        info = ROOT.TPaveText(0.08, 0.15, 0.92, 0.86, "NDC")
        info.SetFillStyle(0); info.SetBorderSize(0); info.SetTextAlign(12)
        info.AddText("t{} [{:.4f}, {:.4f}]  delta{} [{:.2f}, {:.2f}]".format(t_index + 1, cell["t_low"], cell["t_high"], delta_index + 1, cell["delta_low"], cell["delta_high"]))
        info.AddText("support class: {}".format(status))
        for side in _SIDES:
            metrics = cell.get(side) or {}
            info.AddText("{}: abs={:.3e}, Neff={:.2f}, signed={:.3e}".format(side, float(metrics.get("absolute_weight_support", 0.0)), float(metrics.get("effective_entries", 0.0)), float(metrics.get("signed_weight_sum", 0.0))))
        info.AddText("Diagnostic only: no P_pi extraction or cell merge.")
        info.Draw()
        cell_page.cd(4)
        status_box = ROOT.TPaveText(0.08, 0.30, 0.92, 0.70, "NDC")
        status_box.SetFillStyle(0); status_box.SetBorderSize(0); status_box.AddText(DIAGNOSTIC_LABEL); status_box.AddText("pre-HGCer view when the source tree permits it")
        status_box.Draw()
        cell_page._pion_hgcer_cell = (kaon_npe, pion_npe, kaon_mm, pion_mm, info, status_box)
        pages.append(("{}.t{}.delta{}".format(prefix, t_index + 1, delta_index + 1), cell_page, "t{} delta{}".format(t_index + 1, delta_index + 1)))

    emitted = []
    for page_index, (page_id, canvas, scope) in enumerate(pages):
        suffix = ")" if bool(close_pdf) and page_index == len(pages) - 1 else ""
        canvas.Print(str(pdf_name) + suffix)
        entry = {"page_id": page_id, "scope": scope, "authoritative": False}
        page_manifest.append(entry)
        emitted.append(entry)
    expected_manifest = expected_pion_hgcer_page_manifest(payload)
    if [entry["page_id"] for entry in emitted] != [entry["page_id"] for entry in expected_manifest]:
        raise RuntimeError("pion_hgcer_diagnostic_renderer_manifest_mismatch")
    payload["page_manifest"] = list(emitted)
    return emitted
