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


def _resolve_lambda_window(sample_hist, lambda_window=None):
    axis = sample_hist.GetXaxis()
    if lambda_window is None:
        return float(axis.GetXmin()), float(axis.GetXmax())
    return float(lambda_window[0]), float(lambda_window[1])


def _bin_center_from_index(axis, nbins, bin_index):
    if bin_index <= 0:
        return float(axis.GetBinLowEdge(1))
    if bin_index > nbins:
        return float(axis.GetBinUpEdge(nbins))
    return float(axis.GetBinCenter(bin_index))


def build_mm_background_weights_with_diagnostics(sample_hist, background_hist, lambda_window=None):
    nbins = sample_hist.GetNbinsX()
    weights = np.zeros(nbins + 2, dtype=np.float64)
    axis = sample_hist.GetXaxis()
    lambda_lo, lambda_hi = _resolve_lambda_window(sample_hist, lambda_window=lambda_window)
    lambda_sample_yield = 0.0
    oversub_lambda_integral = 0.0
    oversub_count = 0
    max_ratio = 0.0
    affected_mm_centers = []

    for bin_index in range(0, nbins + 2):
        sample_content = float(sample_hist.GetBinContent(bin_index))
        background_content = float(background_hist.GetBinContent(bin_index))
        bin_center = _bin_center_from_index(axis, nbins, bin_index)
        in_lambda_window = lambda_lo <= bin_center <= lambda_hi

        if in_lambda_window and sample_content > 0.0 and math.isfinite(sample_content):
            lambda_sample_yield += sample_content
        if (
            sample_content > 0.0
            and background_content > 0.0
            and math.isfinite(sample_content)
            and math.isfinite(background_content)
        ):
            weight = background_content / sample_content
            if not math.isfinite(weight):
                weight = 0.0
            max_ratio = max(max_ratio, float(weight))
            if weight > 1.0:
                oversub_count += 1
                affected_mm_centers.append(bin_center)
                if in_lambda_window:
                    oversub_lambda_integral += max(background_content - sample_content, 0.0)
            weights[bin_index] = min(max(weight, 0.0), 1.0)

    diagnostics = {
        "oversub_bin_count": int(oversub_count),
        "max_unclamped_ratio": float(max_ratio),
        "oversub_integral": float(oversub_lambda_integral),
        "lambda_window_yield": float(lambda_sample_yield),
        "affected_lambda_fraction": float(oversub_lambda_integral / lambda_sample_yield)
        if lambda_sample_yield > 0.0
        else 0.0,
        "affected_mm_range": None,
    }
    if affected_mm_centers:
        diagnostics["affected_mm_range"] = [
            float(min(affected_mm_centers)),
            float(max(affected_mm_centers)),
        ]
    return weights, diagnostics


def build_mm_background_weights(sample_hist, background_hist):
    weights, _ = build_mm_background_weights_with_diagnostics(sample_hist, background_hist)
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
