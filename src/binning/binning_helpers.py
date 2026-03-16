#! /usr/bin/python

import math

import numpy as np


def find_bin_index(value, bin_edges):
    if not math.isfinite(value):
        return None

    index = int(np.searchsorted(bin_edges, value, side="right") - 1)
    if index < 0 or index >= len(bin_edges) - 1:
        return None

    if bin_edges[index] <= value < bin_edges[index + 1]:
        return index

    return None


def find_2d_bin_indices(primary_value, secondary_value, primary_edges, secondary_edges):
    primary_index = find_bin_index(primary_value, primary_edges)
    if primary_index is None:
        return None, None

    secondary_index = find_bin_index(secondary_value, secondary_edges)
    if secondary_index is None:
        return primary_index, None

    return primary_index, secondary_index
