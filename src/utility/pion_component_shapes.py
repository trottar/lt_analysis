#! /usr/bin/python

from __future__ import annotations

import math
import os
import re
import sys
from copy import deepcopy
from pathlib import Path
from uuid import uuid4

import ROOT
import numpy as np
from ltsep import Misc


MODULE_DIR = Path(__file__).resolve().parent
SRC_DIR = MODULE_DIR.parent
for _path in (MODULE_DIR, SRC_DIR / "cuts", SRC_DIR / "binning"):
    path_text = str(_path)
    if path_text not in sys.path:
        sys.path.append(path_text)

from background_config import (
    BG_OPT_MM_PLOT_MAX,
    BG_OPT_MM_PLOT_MIN,
    BG_OPT_MM_PLOT_NBINS,
    PARTICLE_SUBTRACTION_MODE_COMPONENTS,
    SIMC_PION_COMPONENT_BACKGROUND_MAP,
    resolve_particle_subtraction_mode,
    resolve_simc_pion_component_files,
    resolve_simc_tree_name,
)
from utility import normalize_hist_to_unit_area
from apply_cuts import apply_simc_cuts, apply_simc_sub_cuts, set_val
from binning_helpers import find_2d_bin_indices
from hgcer_hole import apply_HGCer_hole_cut


_SIMC_COMPONENT_SHAPE_CACHE = {}

_KAON_SIMC_REQUIRED_BRANCHES = (
    "missmass",
    "t",
    "ssdelta",
    "ssxptar",
    "ssyptar",
    "hsdelta",
    "hsxptar",
    "hsyptar",
    "Q2",
    "W",
    "phgcer_x_det",
    "phgcer_y_det",
)
_KAON_SIMC_BINNED_REQUIRED_BRANCHES = ("phipq",)


def _sigma0_source_identity(inpDict, phi_setting):
    return {
        "Q2": str((inpDict or {}).get("Q2") or ""),
        "W": str((inpDict or {}).get("W") or ""),
        "EPSSET": str((inpDict or {}).get("EPSSET") or ""),
        "phi_setting": str(phi_setting or ""),
    }


def _sigma0_manifest_entry(inpDict, phi_setting):
    return (
        ((((inpDict or {}).get("background_samples") or {}).get("by_phi") or {})
         .get(phi_setting, {}) or {})
        .get("sigma0", {})
        or {}
    )


def resolve_kaon_simc_sigma0_root_filename(root_filename, inpDict, phi_setting):
    """Resolve only an explicitly configured K-Sigma0 SIMC source.

    K-Sigma0 is an external required input.  In particular, this resolver does
    not search generated-background directories or infer a path from another
    kaon/pion component when no Sigma0 path was supplied.
    """
    sample_entry = _sigma0_manifest_entry(inpDict, phi_setting)
    requested_root = str(root_filename or sample_entry.get("root") or "").strip() or None
    source_identity = deepcopy(
        sample_entry.get("source_identity") or _sigma0_source_identity(inpDict, phi_setting)
    )
    source_strategy = str(sample_entry.get("source_strategy") or "external_required")
    candidate_roots = [requested_root] if requested_root else []
    existing_roots = [candidate for candidate in candidate_roots if os.path.isfile(candidate)]
    if not requested_root:
        resolution_source = "no_source_configured"
        resolved_root = None
    elif existing_roots:
        resolution_source = "explicit_configured_root"
        resolved_root = existing_roots[0]
    else:
        resolution_source = "configured_path_does_not_exist"
        resolved_root = requested_root

    return resolved_root, {
        "requested_root": requested_root,
        "resolved_root": resolved_root,
        "candidate_roots": candidate_roots,
        "existing_roots": existing_roots,
        "rejected_roots": [],
        "fallback_used": False,
        "configured": bool(requested_root),
        "source_strategy": source_strategy,
        "source_identity": source_identity,
        "resolution_source": resolution_source,
    }


def resolve_kaon_simc_signal_root_filename(root_filename, inpDict, phi_setting):
    """Resolve the iter-weighted K-Lambda SIMC product for one phi setting.

    ``Prod_Coin`` is the archived SIMC input.  The K-Lambda comparison must
    instead use the setting-specific kaon file after ``iter_weight`` applies
    the active LT-separation model.
    """
    requested_root = str(root_filename or "").strip()
    candidate_roots = []
    rejected_roots = []

    def add_candidate(candidate):
        candidate = str(candidate or "").strip()
        if candidate and candidate not in candidate_roots:
            candidate_roots.append(candidate)

    def is_archived_prod_coin(candidate):
        return os.path.basename(str(candidate or "")).lower().startswith("prod_coin_")

    output_dir = os.path.dirname(requested_root)
    if not output_dir:
        output_dir = str((inpDict or {}).get("OUTPATH") or "").strip()

    q2 = str((inpDict or {}).get("Q2") or "").strip()
    w = str((inpDict or {}).get("W") or "").strip()
    epsset = str((inpDict or {}).get("EPSSET") or "").strip()
    phi_token = str(phi_setting or "").strip()

    if requested_root:
        if is_archived_prod_coin(requested_root):
            rejected_roots.append(requested_root)
        else:
            add_candidate(requested_root)

    if output_dir and q2 and w and epsset and phi_token:
        add_candidate(
            os.path.join(
                output_dir,
                "{}_kaon_Simc_Q{}W{}_{}e.root".format(
                    phi_token,
                    q2,
                    w,
                    epsset,
                ),
            )
        )

    existing_roots = [candidate for candidate in candidate_roots if os.path.isfile(candidate)]
    resolved_root = existing_roots[0] if existing_roots else (candidate_roots[0] if candidate_roots else None)
    resolution = {
        "requested_root": requested_root or None,
        "resolved_root": resolved_root,
        "candidate_roots": candidate_roots,
        "existing_roots": existing_roots,
        "fallback_used": bool(existing_roots and existing_roots[0] != requested_root),
        "rejected_roots": rejected_roots,
    }
    return resolved_root, resolution


def _sanitize_token(value):
    token = re.sub(r"[^A-Za-z0-9_]+", "_", str(value or "").strip())
    return token.strip("_") or "component"


def _build_hist_name(prefix, component_name, phi_setting, context):
    return "{}_{}_{}_{}_{}".format(
        prefix,
        _sanitize_token(component_name),
        _sanitize_token(phi_setting),
        _sanitize_token(context),
        uuid4().hex[:10],
    )


def _warn_component_load(component_name, phi_setting, reason, **details):
    print("WARNING: SIMC template shape load issue")
    print("  component_name = {}".format(component_name))
    print("  phi_setting = {}".format(phi_setting))
    print("  reason = {}".format(reason))
    for key, value in details.items():
        print("  {} = {}".format(key, value))


def _display_component_name(component_name):
    return str(component_name).replace("_", "-")


def _make_component_hist(component_name, phi_setting, context, prefix, n_bins, x_min, x_max):
    hist = ROOT.TH1D(
        _build_hist_name(prefix, component_name, phi_setting, context),
        "{} {} {}".format(prefix, component_name, phi_setting),
        int(n_bins),
        float(x_min),
        float(x_max),
    )
    hist.SetDirectory(0)
    return hist


def _clone_detached_hist(hist):
    if hist is None:
        return None
    cloned = hist.Clone("{}_{}".format(hist.GetName(), uuid4().hex[:10]))
    if hasattr(cloned, "SetDirectory"):
        cloned.SetDirectory(0)
    return cloned


def _clone_binned_shape_map(binned_shapes):
    cloned_shapes = {}
    for t_key, phi_shape_map in (binned_shapes or {}).items():
        cloned_shapes[t_key] = {}
        for phi_key, hist in (phi_shape_map or {}).items():
            cloned_shapes[t_key][phi_key] = _clone_detached_hist(hist)
    return cloned_shapes


def _clone_component_payload(payload):
    if not isinstance(payload, dict):
        return payload
    return {
        "setting_shape": _clone_detached_hist(payload.get("setting_shape")),
        "setting_shape_full": _clone_detached_hist(payload.get("setting_shape_full")),
        "binned_shapes": _clone_binned_shape_map(payload.get("binned_shapes") or {}),
        "diagnostics": deepcopy(payload.get("diagnostics") or {}),
    }


def _freeze_bin_edges_for_cache(edges):
    if edges is None:
        return None
    return tuple(round(float(value), 8) for value in np.asarray(edges, dtype=float))


def _freeze_cache_value(value):
    if isinstance(value, dict):
        return tuple(
            (str(key), _freeze_cache_value(item))
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        )
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_cache_value(item) for item in value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def _get_root_file_cache_token(root_filename):
    try:
        stat_result = os.stat(root_filename)
        return (
            os.path.abspath(root_filename),
            int(getattr(stat_result, "st_mtime_ns", int(stat_result.st_mtime * 1e9))),
            int(stat_result.st_size),
        )
    except OSError:
        return (os.path.abspath(root_filename), None, None)


def _build_component_shape_cache_key(
    root_filename,
    tree_name,
    phi_setting,
    particle_type,
    component_name,
    mm_min,
    mm_max,
    mm_plot_min,
    mm_plot_max,
    mm_plot_nbins,
    t_bins,
    phi_bins,
    use_full_mm_range,
    cache_identity=None,
):
    return (
        _get_root_file_cache_token(root_filename),
        str(tree_name),
        str(phi_setting),
        str(particle_type).strip().lower(),
        str(component_name),
        round(float(mm_min), 8),
        round(float(mm_max), 8),
        round(float(mm_plot_min), 8),
        round(float(mm_plot_max), 8),
        int(mm_plot_nbins),
        _freeze_bin_edges_for_cache(t_bins),
        _freeze_bin_edges_for_cache(phi_bins),
        bool(use_full_mm_range),
        _freeze_cache_value(cache_identity),
    )


def _build_binned_shape_map(component_name, phi_setting, context, n_bins, x_min, x_max, t_bins, phi_bins):
    shapes = {}
    if t_bins is None or phi_bins is None or len(t_bins) < 2 or len(phi_bins) < 2:
        return shapes

    for t_index in range(len(t_bins) - 1):
        t_key = "t_bin{}".format(t_index + 1)
        shapes[t_key] = {}
        for phi_index in range(len(phi_bins) - 1):
            phi_key = "phi_bin{}".format(phi_index + 1)
            shapes[t_key][phi_key] = _make_component_hist(
                component_name,
                phi_setting,
                "{}_{}_{}".format(context, t_key, phi_key),
                "H_MM_component_shape_full",
                n_bins,
                x_min,
                x_max,
            )
    return shapes


def _empty_component_payload(
    component_name,
    root_filename,
    tree_name,
    phi_setting,
    context,
    setting_hist,
    setting_full_hist,
    fallback_reason,
    t_bins=None,
    phi_bins=None,
    diagnostic_updates=None,
):
    diagnostics = {
        "component_name": component_name,
        "root_filename": root_filename,
        "tree_name": tree_name,
        "path_exists": bool(root_filename and os.path.isfile(root_filename)),
        "root_open_success": None,
        "tree_exists": None,
        "tree_entries": None,
        "required_branches": [],
        "missing_required_branches": [],
        "n_events_seen": 0,
        "n_events_passed": 0,
        "n_events_passed_mm_window": 0,
        "weighted_integral_before_norm": 0.0,
        "weighted_integral_after_norm": 0.0,
        "setting_shape_integral_before_norm": 0.0,
        "setting_shape_integral_after_norm": 0.0,
        "setting_shape_normalized": False,
        "normalized": False,
        "fallback_used": True,
        "fallback_reason": fallback_reason,
        "n_binned_shapes": 0 if t_bins is None or phi_bins is None else (len(t_bins) - 1) * (len(phi_bins) - 1),
        "n_normalized_binned_shapes": 0,
    }
    if diagnostic_updates:
        diagnostics.update(deepcopy(diagnostic_updates))
    return {
        "setting_shape": setting_hist,
        "setting_shape_full": setting_full_hist,
        "binned_shapes": _build_binned_shape_map(
            component_name,
            phi_setting,
            context,
            setting_full_hist.GetNbinsX(),
            setting_full_hist.GetXaxis().GetXmin(),
            setting_full_hist.GetXaxis().GetXmax(),
            t_bins,
            phi_bins,
        ),
        "diagnostics": diagnostics,
    }


def load_pion_simc_component_shape(
    root_filename,
    inpDict,
    phi_setting,
    particle_type,
    component_name,
    t_bins=None,
    phi_bins=None,
    hgcer_cutg=None,
    use_full_mm_range=True,
    context="",
    source_provenance=None,
    required_branch_names=None,
    failure_reasons=None,
    cache_identity=None,
):
    tree_name = resolve_simc_tree_name(inpDict)
    mm_min = float(inpDict["mm_min"])
    mm_max = float(inpDict["mm_max"])
    mm_plot_min = float(inpDict.get("bg_opt_mm_plot_min", BG_OPT_MM_PLOT_MIN))
    mm_plot_max = float(inpDict.get("bg_opt_mm_plot_max", BG_OPT_MM_PLOT_MAX))
    mm_plot_nbins = int(inpDict.get("bg_opt_mm_plot_nbins", BG_OPT_MM_PLOT_NBINS))
    source_provenance = deepcopy(source_provenance or {})
    failure_reasons = dict(failure_reasons or {})
    required_branch_names = tuple(required_branch_names or ())
    resolution_fallback_used = source_provenance.pop("fallback_used", None)
    path_exists = bool(root_filename and os.path.isfile(root_filename))
    base_diagnostic_updates = {
        "path_exists": path_exists,
        "root_open_success": None,
        "tree_exists": None,
        "tree_entries": None,
        "required_branches": list(required_branch_names),
        "missing_required_branches": [],
    }
    if source_provenance:
        base_diagnostic_updates.update(source_provenance)
    if resolution_fallback_used is not None:
        base_diagnostic_updates["resolution_fallback_used"] = bool(resolution_fallback_used)

    cache_key = None
    if root_filename and hgcer_cutg is None:
        cache_key = _build_component_shape_cache_key(
            root_filename,
            tree_name,
            phi_setting,
            particle_type,
            component_name,
            mm_min,
            mm_max,
            mm_plot_min,
            mm_plot_max,
            mm_plot_nbins,
            t_bins,
            phi_bins,
            use_full_mm_range,
            cache_identity=cache_identity,
        )
        cached_payload = _SIMC_COMPONENT_SHAPE_CACHE.get(cache_key)
        if cached_payload is not None:
            print(
                "[SIMC TEMPLATE CACHE] {} {} file={}".format(
                    phi_setting,
                    component_name,
                    root_filename,
                )
            )
            return _clone_component_payload(cached_payload)

    full_nbins = mm_plot_nbins if use_full_mm_range else 100
    full_xmin = mm_plot_min if use_full_mm_range else mm_min
    full_xmax = mm_plot_max if use_full_mm_range else mm_max

    setting_shape = _make_component_hist(
        component_name,
        phi_setting,
        context,
        "H_MM_component_shape",
        100,
        mm_min,
        mm_max,
    )
    setting_shape_full = _make_component_hist(
        component_name,
        phi_setting,
        context,
        "H_MM_component_shape_full",
        full_nbins,
        full_xmin,
        full_xmax,
    )
    binned_shapes = _build_binned_shape_map(
        component_name,
        phi_setting,
        context,
        full_nbins,
        full_xmin,
        full_xmax,
        t_bins,
        phi_bins,
    )

    if not root_filename or not path_exists:
        reason_key = "no_source_configured" if not root_filename else "configured_path_does_not_exist"
        fallback_reason = failure_reasons.get(reason_key, "missing ROOT file")
        _warn_component_load(
            component_name,
            phi_setting,
            fallback_reason,
            root_filename=root_filename,
        )
        return _empty_component_payload(
            component_name,
            root_filename,
            tree_name,
            phi_setting,
            context,
            setting_shape,
            setting_shape_full,
            fallback_reason,
            t_bins=t_bins,
            phi_bins=phi_bins,
            diagnostic_updates=base_diagnostic_updates,
        )

    input_file = ROOT.TFile.Open(root_filename, "READ")
    if not input_file or input_file.IsZombie():
        fallback_reason = failure_reasons.get("root_open_failed", "unable to open ROOT file")
        _warn_component_load(
            component_name,
            phi_setting,
            fallback_reason,
            root_filename=root_filename,
        )
        return _empty_component_payload(
            component_name,
            root_filename,
            tree_name,
            phi_setting,
            context,
            setting_shape,
            setting_shape_full,
            fallback_reason,
            t_bins=t_bins,
            phi_bins=phi_bins,
            diagnostic_updates=dict(base_diagnostic_updates, root_open_success=False),
        )

    tree_simc = input_file.Get(tree_name)
    if not tree_simc:
        input_file.Close()
        fallback_reason = failure_reasons.get("missing_simc_tree", "missing SIMC tree")
        _warn_component_load(
            component_name,
            phi_setting,
            fallback_reason,
            root_filename=root_filename,
            tree_name=tree_name,
        )
        return _empty_component_payload(
            component_name,
            root_filename,
            tree_name,
            phi_setting,
            context,
            setting_shape,
            setting_shape_full,
            fallback_reason,
            t_bins=t_bins,
            phi_bins=phi_bins,
            diagnostic_updates=dict(
                base_diagnostic_updates,
                root_open_success=True,
                tree_exists=False,
            ),
        )

    tree_entries = int(tree_simc.GetEntries())
    missing_required_branches = [
        branch_name
        for branch_name in required_branch_names
        if not tree_simc.GetBranch(branch_name) and not tree_simc.GetLeaf(branch_name)
    ]
    tree_diagnostic_updates = dict(
        base_diagnostic_updates,
        root_open_success=True,
        tree_exists=True,
        tree_entries=tree_entries,
        missing_required_branches=missing_required_branches,
    )
    if missing_required_branches:
        input_file.Close()
        fallback_reason = failure_reasons.get(
            "incompatible_tree_missing_branches",
            "incompatible SIMC tree missing required branches",
        )
        _warn_component_load(
            component_name,
            phi_setting,
            fallback_reason,
            root_filename=root_filename,
            tree_name=tree_name,
            missing_required_branches=missing_required_branches,
        )
        return _empty_component_payload(
            component_name,
            root_filename,
            tree_name,
            phi_setting,
            context,
            setting_shape,
            setting_shape_full,
            fallback_reason,
            t_bins=t_bins,
            phi_bins=phi_bins,
            diagnostic_updates=tree_diagnostic_updates,
        )
    if required_branch_names and tree_entries == 0:
        input_file.Close()
        fallback_reason = failure_reasons.get("zero_entry_tree", "SIMC tree has zero entries")
        _warn_component_load(
            component_name,
            phi_setting,
            fallback_reason,
            root_filename=root_filename,
            tree_name=tree_name,
        )
        return _empty_component_payload(
            component_name,
            root_filename,
            tree_name,
            phi_setting,
            context,
            setting_shape,
            setting_shape_full,
            fallback_reason,
            t_bins=t_bins,
            phi_bins=phi_bins,
            diagnostic_updates=tree_diagnostic_updates,
        )

    set_val(inpDict)

    hole_cut = hgcer_cutg
    if str(particle_type).strip().lower() == "kaon" and hole_cut is None:
        hole_cut = apply_HGCer_hole_cut(
            inpDict["Q2"],
            inpDict["W"],
            inpDict["EPSSET"],
            phi_setting,
        )

    t_bin_edges = None if t_bins is None else np.asarray(t_bins, dtype=float)
    phi_bin_edges = None if phi_bins is None else np.asarray(phi_bins, dtype=float)
    n_events_seen = 0
    n_events_passed = 0
    n_events_passed_mm_window = 0
    total_entries = tree_entries

    print(
        "\nGrabbing {} {} SIMC template...".format(
            phi_setting,
            _display_component_name(component_name),
        )
    )

    for evt in tree_simc:
        Misc.progressBar(n_events_seen, total_entries, bar_length=25)
        n_events_seen += 1
        adj_missmass = evt.missmass
        adj_t = -evt.t
        weight = float(getattr(evt, "iter_weight", 1.0))
        if not math.isfinite(weight):
            continue

        base_cuts = apply_simc_cuts(evt, mm_min, mm_max)
        base_sub_cuts = apply_simc_sub_cuts(evt, mm_min, mm_max)
        if str(particle_type).strip().lower() == "kaon":
            outside_hole = not hole_cut.IsInside(evt.phgcer_x_det, evt.phgcer_y_det)
            allcuts = base_cuts and outside_hole
            nommcuts = base_sub_cuts and outside_hole
        else:
            allcuts = base_cuts
            nommcuts = base_sub_cuts

        if nommcuts:
            n_events_passed += 1
            setting_shape_full.Fill(adj_missmass, weight)
            if t_bin_edges is not None and phi_bin_edges is not None:
                phi_shift_deg = float(evt.phipq) * (180.0 / math.pi)
                t_index, phi_index = find_2d_bin_indices(
                    adj_t,
                    phi_shift_deg,
                    t_bin_edges,
                    phi_bin_edges,
                )
                if t_index is not None and phi_index is not None:
                    binned_shapes["t_bin{}".format(t_index + 1)][
                        "phi_bin{}".format(phi_index + 1)
                    ].Fill(adj_missmass, weight)

        if allcuts:
            n_events_passed_mm_window += 1
            setting_shape.Fill(adj_missmass, weight)

    input_file.Close()

    setting_integral_before_norm = float(setting_shape.Integral())
    full_integral_before_norm = float(setting_shape_full.Integral())
    setting_normalized = normalize_hist_to_unit_area(
        setting_shape,
        quiet=True,
        context="{} {} setting shape".format(phi_setting, component_name),
    )
    full_normalized = normalize_hist_to_unit_area(
        setting_shape_full,
        quiet=True,
        context="{} {} full shape".format(phi_setting, component_name),
    )

    normalized_binned_shapes = 0
    total_binned_shapes = 0
    for t_key, phi_shape_map in binned_shapes.items():
        for phi_key, binned_hist in phi_shape_map.items():
            total_binned_shapes += 1
            if normalize_hist_to_unit_area(
                binned_hist,
                quiet=True,
                context="{} {} {} {}".format(phi_setting, component_name, t_key, phi_key),
            ):
                normalized_binned_shapes += 1

    fallback_reason = ""
    fallback_used = False
    if n_events_passed == 0:
        fallback_used = True
        fallback_reason = failure_reasons.get(
            "no_events_passed_component_shape_cuts",
            "no SIMC events passed component-shape cuts",
        )
    elif full_integral_before_norm <= 0.0:
        fallback_used = True
        fallback_reason = failure_reasons.get(
            "weighted_integral_non_positive",
            "component full-shape integral was non-positive",
        )
    elif not full_normalized:
        fallback_used = True
        fallback_reason = failure_reasons.get(
            "normalization_failed",
            "component full-shape normalization failed",
        )

    if fallback_used:
        _warn_component_load(
            component_name,
            phi_setting,
            fallback_reason,
            root_filename=root_filename,
            tree_name=tree_name,
            n_events_seen=n_events_seen,
            n_events_passed=n_events_passed,
            n_events_passed_mm_window=n_events_passed_mm_window,
            weighted_integral_before_norm=full_integral_before_norm,
        )

    diagnostics = {
        "component_name": component_name,
        "root_filename": root_filename,
        "tree_name": tree_name,
        "n_events_seen": n_events_seen,
        "n_events_passed": n_events_passed,
        "n_events_passed_mm_window": n_events_passed_mm_window,
        "weighted_integral_before_norm": full_integral_before_norm,
        "weighted_integral_after_norm": float(setting_shape_full.Integral()),
        "setting_shape_integral_before_norm": setting_integral_before_norm,
        "setting_shape_integral_after_norm": float(setting_shape.Integral()),
        "setting_shape_normalized": bool(setting_normalized),
        "normalized": bool(full_normalized),
        "fallback_used": fallback_used,
        "fallback_reason": fallback_reason,
        "n_binned_shapes": total_binned_shapes,
        "n_normalized_binned_shapes": normalized_binned_shapes,
    }
    diagnostics.update(tree_diagnostic_updates)

    print(
        "[SIMC TEMPLATE] {} {} tree={} seen={} passed={} mm_passed={} full_before={:.6e} full_after={:.6e}".format(
            phi_setting,
            component_name,
            tree_name,
            n_events_seen,
            n_events_passed,
            n_events_passed_mm_window,
            full_integral_before_norm,
            diagnostics["weighted_integral_after_norm"],
        )
    )

    payload = {
        "setting_shape": setting_shape,
        "setting_shape_full": setting_shape_full,
        "binned_shapes": binned_shapes,
        "diagnostics": diagnostics,
    }
    if cache_key is not None:
        _SIMC_COMPONENT_SHAPE_CACHE[cache_key] = _clone_component_payload(payload)
    return _clone_component_payload(payload) if cache_key is not None else payload


def load_setting_pion_component_shapes(
    inpDict,
    phi_setting,
    particle_type=None,
    t_bins=None,
    phi_bins=None,
    hgcer_cutg=None,
    use_full_mm_range=True,
    context="",
):
    mode = resolve_particle_subtraction_mode(inpDict)
    component_files = resolve_simc_pion_component_files(inpDict, phi_setting)
    if mode != PARTICLE_SUBTRACTION_MODE_COMPONENTS:
        return {
            "mode": mode,
            "tree_name": resolve_simc_tree_name(inpDict),
            "component_files": component_files,
            "components": {},
            "diagnostics": {},
        }

    resolved_particle_type = particle_type or inpDict.get("ParticleType", "")
    component_payloads = {}
    diagnostics = {}
    for component_name in SIMC_PION_COMPONENT_BACKGROUND_MAP:
        payload = load_pion_simc_component_shape(
            component_files.get(component_name),
            inpDict,
            phi_setting,
            resolved_particle_type,
            component_name,
            t_bins=t_bins,
            phi_bins=phi_bins,
            hgcer_cutg=hgcer_cutg,
            use_full_mm_range=use_full_mm_range,
            context=context,
        )
        component_payloads[component_name] = payload
        diagnostics[component_name] = deepcopy(payload["diagnostics"])

    return {
        "mode": mode,
        "tree_name": resolve_simc_tree_name(inpDict),
        "component_files": component_files,
        "components": component_payloads,
        "diagnostics": diagnostics,
    }


def load_kaon_simc_extra_shape(
    root_filename,
    inpDict,
    phi_setting,
    component_name,
    t_bins=None,
    phi_bins=None,
    hgcer_cutg=None,
    use_full_mm_range=True,
    context="",
):
    return load_pion_simc_component_shape(
        root_filename,
        inpDict,
        phi_setting,
        "kaon",
        component_name,
        t_bins=t_bins,
        phi_bins=phi_bins,
        hgcer_cutg=hgcer_cutg,
        use_full_mm_range=use_full_mm_range,
        context=context,
    )


def load_kaon_simc_signal_shape(
    root_filename,
    inpDict,
    phi_setting,
    t_bins=None,
    phi_bins=None,
    hgcer_cutg=None,
    use_full_mm_range=True,
    context="",
):
    resolved_root_filename, root_resolution = resolve_kaon_simc_signal_root_filename(
        root_filename,
        inpDict,
        phi_setting,
    )
    print(
        "[SIMC K-Lambda] {} requested={} resolved={} fallback={}".format(
            phi_setting,
            root_resolution.get("requested_root"),
            root_resolution.get("resolved_root"),
            root_resolution.get("fallback_used"),
        )
    )
    payload = load_kaon_simc_extra_shape(
        resolved_root_filename,
        inpDict,
        phi_setting,
        "k_lambda_signal",
        t_bins=t_bins,
        phi_bins=phi_bins,
        hgcer_cutg=hgcer_cutg,
        use_full_mm_range=use_full_mm_range,
        context=context,
    )
    if isinstance(payload, dict):
        diagnostics = dict(payload.get("diagnostics") or {})
        diagnostics["root_resolution"] = deepcopy(root_resolution)
        payload["diagnostics"] = diagnostics
    return payload


def load_kaon_simc_sigma0_shape(
    root_filename,
    inpDict,
    phi_setting,
    t_bins=None,
    phi_bins=None,
    hgcer_cutg=None,
    use_full_mm_range=True,
    context="",
):
    resolved_root_filename, root_resolution = resolve_kaon_simc_sigma0_root_filename(
        root_filename,
        inpDict,
        phi_setting,
    )
    required_branches = list(_KAON_SIMC_REQUIRED_BRANCHES)
    if t_bins is not None and phi_bins is not None:
        required_branches.extend(_KAON_SIMC_BINNED_REQUIRED_BRANCHES)
    payload = load_pion_simc_component_shape(
        resolved_root_filename,
        inpDict,
        phi_setting,
        "kaon",
        "k_sigma0_signal",
        t_bins=t_bins,
        phi_bins=phi_bins,
        hgcer_cutg=hgcer_cutg,
        use_full_mm_range=use_full_mm_range,
        context=context,
        source_provenance=root_resolution,
        required_branch_names=required_branches,
        failure_reasons={
            "no_source_configured": "no_source_configured",
            "configured_path_does_not_exist": "configured_path_does_not_exist",
            "root_open_failed": "root_open_failed",
            "missing_simc_tree": "missing_simc_tree",
            "incompatible_tree_missing_branches": "incompatible_tree_missing_branches",
            "zero_entry_tree": "zero_entry_tree",
            "no_events_passed_component_shape_cuts": "no_events_passed_component_shape_cuts",
            "weighted_integral_non_positive": "weighted_integral_non_positive",
            "normalization_failed": "normalization_failed",
        },
        cache_identity=root_resolution.get("source_identity"),
    )
    diagnostics = dict(payload.get("diagnostics") or {}) if isinstance(payload, dict) else {}
    if diagnostics:
        diagnostics["root_resolution"] = deepcopy(root_resolution)
        payload["diagnostics"] = diagnostics
        identity = diagnostics.get("source_identity") or {}
        prefix = "[SIMC K-SIGMA0] phi={} Q2={} W={} EPSSET={}".format(
            identity.get("phi_setting", phi_setting),
            identity.get("Q2", ""),
            identity.get("W", ""),
            identity.get("EPSSET", ""),
        )
        if diagnostics.get("fallback_used"):
            print(
                "{} UNAVAILABLE requested={} resolved={} source={} exists={} reason={}".format(
                    prefix,
                    diagnostics.get("requested_root"),
                    diagnostics.get("resolved_root"),
                    diagnostics.get("resolution_source"),
                    diagnostics.get("path_exists"),
                    diagnostics.get("fallback_reason"),
                )
            )
        else:
            print(
                "{} requested={} resolved={} source={} tree={} entries={} seen={} passed={} mm_passed={} normalized={}".format(
                    prefix,
                    diagnostics.get("requested_root"),
                    diagnostics.get("resolved_root"),
                    diagnostics.get("resolution_source"),
                    diagnostics.get("tree_name"),
                    diagnostics.get("tree_entries"),
                    diagnostics.get("n_events_seen"),
                    diagnostics.get("n_events_passed"),
                    diagnostics.get("n_events_passed_mm_window"),
                    diagnostics.get("normalized"),
                )
            )
    return payload


def attach_pion_component_payload(hist_dict, component_payload):
    if not isinstance(hist_dict, dict) or not isinstance(component_payload, dict):
        return hist_dict

    hist_dict["_simc_pion_component_payload"] = component_payload
    hist_dict["simc_pion_component_files"] = deepcopy(
        component_payload.get("component_files") or {}
    )
    hist_dict["simc_pion_component_diagnostics"] = deepcopy(
        component_payload.get("diagnostics") or {}
    )
    return hist_dict
