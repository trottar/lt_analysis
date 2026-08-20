#! /usr/bin/python
"""Small, ROOT-free ownership helpers for retained histogram pipelines.

The helper deliberately keeps no reference to a ROOT object.  ROOT ownership is
made explicit at every clone boundary: clone, detach from the current ROOT
directory, then hand the object to exactly one result container.
"""

from __future__ import annotations

import itertools
import hashlib
import math
import re
import struct


_CLONE_COUNTER = itertools.count(1)
_DEBUG_ENABLED = False
_DEBUG_RECORDS = []


def fingerprint_histogram_content_error(hist):
    """Return a stable SHA-256 fingerprint for a TH1-like detached product.

    Histogram ownership is intentionally independent of this helper.  The
    fingerprint is a scalar producer/consumer contract and covers the axis,
    every bin (including under/overflow), and its stored uncertainty.
    """
    if hist is None:
        raise RuntimeError("cannot fingerprint missing histogram")
    try:
        axis = hist.GetXaxis()
        nbins = int(hist.GetNbinsX())
        values = [float(nbins), float(axis.GetXmin()), float(axis.GetXmax())]
        # Xmin/Xmax alone does not identify variable-width ROOT axes.
        for bin_index in range(1, nbins + 2):
            try:
                edge = float(axis.GetBinLowEdge(bin_index))
            except Exception:
                edge = float(axis.GetXmin()) + (
                    float(axis.GetXmax()) - float(axis.GetXmin())
                ) * float(bin_index - 1) / float(max(nbins, 1))
            values.append(edge)
        for bin_index in range(0, nbins + 2):
            values.extend((
                float(hist.GetBinContent(bin_index)),
                float(hist.GetBinError(bin_index)),
            ))
    except Exception as exc:
        raise RuntimeError("cannot fingerprint histogram: {}".format(exc))
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError("cannot fingerprint histogram with nonfinite content or error")
    digest = hashlib.sha256()
    for value in values:
        digest.update(struct.pack("!d", value))
    return digest.hexdigest()


def _safe_name(value, fallback):
    text = str(value or fallback).strip()
    text = re.sub(r"[^A-Za-z0-9_]+", "_", text)
    return text.strip("_") or fallback


def configure_particle_subtraction_root_ownership_debug(enabled=False):
    """Enable lightweight clone provenance records without retaining ROOT objects."""
    global _DEBUG_ENABLED
    _DEBUG_ENABLED = bool(enabled)
    _DEBUG_RECORDS.clear()


def get_particle_subtraction_root_ownership_debug_records():
    """Return scalar debug records suitable for tests or terminal diagnostics."""
    return [dict(record) for record in _DEBUG_RECORDS]


def unique_root_object_name(name, *, scope="scope", role="object"):
    """Allocate a collision-free ROOT object name without retaining an object."""
    name_base = _safe_name(name, "root_object")
    return "{}_{}_{}_{}".format(
        name_base,
        _safe_name(scope, "scope"),
        _safe_name(role, "object"),
        next(_CLONE_COUNTER),
    )


def clone_root_histogram(
    template_hist,
    *,
    scope="scope",
    role="histogram",
    optional=False,
    reset=False,
    sumw2=None,
    name=None,
    title=None,
):
    """Create a uniquely named, detached ROOT histogram clone.

    ``optional=True`` is the sole silent-none path.  Required histogram
    ownership violations fail at the clone boundary instead of surfacing later
    as a ROOT or PyROOT lifetime error.  ``sumw2`` is opt-in so callers retain
    their existing statistical behavior exactly.
    """
    if template_hist is None:
        if optional:
            return None
        raise ValueError("Cannot clone required ROOT histogram for {}:{}".format(scope, role))

    source_name = _safe_name(getattr(template_hist, "GetName", lambda: "source")(), "source")
    clone_base = _safe_name(name, "{}_{}_{}".format(scope, role, source_name))
    unique_name = unique_root_object_name(
        clone_base,
        scope=scope,
        role=role,
    )
    cloned = template_hist.Clone(unique_name)
    if cloned is None:
        raise RuntimeError("ROOT failed to clone histogram '{}'".format(source_name))
    if hasattr(cloned, "SetDirectory"):
        cloned.SetDirectory(0)
    if sumw2 is True and hasattr(cloned, "Sumw2"):
        cloned.Sumw2()
    if reset:
        cloned.Reset()
    if title is not None and hasattr(cloned, "SetTitle"):
        cloned.SetTitle(str(title))

    if _DEBUG_ENABLED:
        directory_status = None
        if hasattr(cloned, "GetDirectory"):
            try:
                directory_status = cloned.GetDirectory() is None
            except Exception:
                directory_status = None
        _DEBUG_RECORDS.append(
            {
                "python_id": int(id(cloned)),
                "root_name": str(cloned.GetName()) if hasattr(cloned, "GetName") else unique_name,
                "scope": str(scope),
                "role": str(role),
                "source_name": source_name,
                "detached": directory_status,
            }
        )
    return cloned
