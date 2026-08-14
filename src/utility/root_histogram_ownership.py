#! /usr/bin/python
"""Small, ROOT-free ownership helpers for retained histogram pipelines.

The helper deliberately keeps no reference to a ROOT object.  ROOT ownership is
made explicit at every clone boundary: clone, detach from the current ROOT
directory, then hand the object to exactly one result container.
"""

from __future__ import annotations

import itertools
import re


_CLONE_COUNTER = itertools.count(1)
_DEBUG_ENABLED = False
_DEBUG_RECORDS = []


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
    unique_name = "{}_{}_{}_{}".format(
        clone_base,
        _safe_name(scope, "scope"),
        _safe_name(role, "histogram"),
        next(_CLONE_COUNTER),
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
