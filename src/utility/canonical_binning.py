#! /usr/bin/python
"""Canonical interval membership shared across particle-stage consumers."""

from __future__ import annotations

import math


def find_canonical_bin(value, edges):
    """Return the canonical bin for ``[low, high)`` intervals.

    The final configured endpoint belongs to the final bin.  Invalid values,
    malformed edges, and values outside the closed overall domain return -1.
    """
    try:
        numeric_value = float(value)
        numeric_edges = [float(edge) for edge in edges]
    except (TypeError, ValueError):
        return -1
    if (
        len(numeric_edges) < 2
        or not math.isfinite(numeric_value)
        or not all(math.isfinite(edge) for edge in numeric_edges)
    ):
        return -1
    if numeric_value < numeric_edges[0] or numeric_value > numeric_edges[-1]:
        return -1
    if numeric_value == numeric_edges[-1]:
        return len(numeric_edges) - 2
    for index, (low_edge, high_edge) in enumerate(zip(numeric_edges[:-1], numeric_edges[1:])):
        if low_edge <= numeric_value < high_edge:
            return index
    return -1
