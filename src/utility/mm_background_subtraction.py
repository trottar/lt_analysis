#! /usr/bin/python

#
# Utilities for propagating MM-fit background estimates into non-MM histograms.
#

import math
import numpy as np


def clone_reset_hist(hist, suffix):
    template = hist.Clone("{}{}".format(hist.GetName(), suffix))
    if hasattr(template, "SetDirectory"):
        template.SetDirectory(0)
    template.Reset()
    return template


def build_mm_background_weights(sample_hist, background_hist):
    nbins = sample_hist.GetNbinsX()
    weights = np.zeros(nbins + 2, dtype=np.float64)

    for bin_index in range(0, nbins + 2):
        sample_content = float(sample_hist.GetBinContent(bin_index))
        background_content = float(background_hist.GetBinContent(bin_index))

        if (
            sample_content > 0.0
            and background_content > 0.0
            and math.isfinite(sample_content)
            and math.isfinite(background_content)
        ):
            weight = background_content / sample_content
            if not math.isfinite(weight):
                weight = 0.0
            weights[bin_index] = min(max(weight, 0.0), 1.0)

    return weights


def build_mm_residual_weights(background_weights):
    return np.clip(1.0 - np.asarray(background_weights, dtype=np.float64), 0.0, None)


def mm_background_weight_from_value(adj_mm, reference_hist, background_weights, residual_weights=None):
    axis = reference_hist.GetXaxis()
    bin_index = axis.FindBin(adj_mm)
    bin_index = max(0, min(reference_hist.GetNbinsX() + 1, bin_index))

    weight = float(background_weights[bin_index])
    if residual_weights is not None:
        weight *= float(residual_weights[bin_index])
    return weight
