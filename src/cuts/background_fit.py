#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-09-18 02:57:48 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
##################################################################################################################################################
# Import relevant packages

import uproot as up
import numpy as np
import root_numpy as rnp
import ROOT
from ROOT import TF1
import sys, os, math

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

################################################################################################################################################

def get_fit_histogram_padded(fit_func,
                             hist,
                             mm_min,
                             mm_max,
                             n_pad: int = 0,
                             allow_negative: bool = False,
                             sb_left: tuple = (1.070, 1.095),
                             sb_right: tuple = (1.135, 1.180)):
    """
    Build a background histogram from *fit_func* that matches the statistics
    of the *CUT* histogram *hist* inside [mm_min, mm_max] (± n_pad bins).

    Parameters
    ----------
    fit_func   : ROOT.TF1
        Function already fitted on the **uncut** MM spectrum.
    hist       : ROOT.TH1
        Histogram *after* all analysis cuts (same binning as the fit).
    mm_min, mm_max : float
        MM window we keep.  Outside that (plus padding) the output is zero.
    n_pad : int, default 0
        Extra *bins* on each side that are also kept.
    allow_negative : bool, default False
        Clamp negative predictions to zero if False.
    sb_left, sb_right : tuple(float,float)
        Ranges that define the two side-bands (must lie inside hist axis).
    """
    # ------------------------------------------------------------------
    # 0.  Re-normalise the fit so its _current_ side-band integral
    #     equals the counts that survive in the cut histogram
    # ------------------------------------------------------------------
    def sideband_counts(h, rng):
        i1 = h.FindBin(rng[0])
        i2 = h.FindBin(rng[1])
        return h.Integral(i1, i2)

    # counts in the CUT histogram
    obs_SB = sideband_counts(hist, sb_left) + sideband_counts(hist, sb_right)

    # expected counts from the fit (density → counts-per-bin)
    bw      = hist.GetBinWidth(1)  # constant in your spectra
    exp_SB  = (fit_func.Integral(*sb_left) + fit_func.Integral(*sb_right)) / bw

    if exp_SB > 0.0:
        scale = obs_SB / exp_SB
        for ip in range(fit_func.GetNpar()):
            fit_func.SetParameter(ip, fit_func.GetParameter(ip) * scale)
    # If exp_SB==0 just leave the parameters untouched (nothing to scale to)

    # ------------------------------------------------------------------
    # 1.  Prepare an empty clone (wipe contents, errors & stats)
    # ------------------------------------------------------------------
    h_bg = hist.Clone(hist.GetName() + "_bg_fit_pad")
    h_bg.Reset("ICES")               # I C E S  → clean slate incl. Sumw2

    ax = h_bg.GetXaxis()
    nb = h_bg.GetNbinsX()

    i_lo = max(1,  ax.FindBin(mm_min) - n_pad)
    i_hi = min(nb, ax.FindBin(mm_max) + n_pad)

    # ------------------------------------------------------------------
    # 2.  Fill the background histogram
    # ------------------------------------------------------------------
    for ib in range(1, nb + 1):
        if i_lo <= ib <= i_hi:
            x      = ax.GetBinCenter(ib)
            y_pred = fit_func.Eval(x)             # counts *per bin* (no bw factor)

            if (not allow_negative) and (y_pred < 0.0):
                y_pred = 0.0

            h_bg.SetBinContent(ib, y_pred)
            h_bg.SetBinError  (ib, 0.0)           # propagate param errors if needed
        else:
            h_bg.SetBinContent(ib, 0.0)
            h_bg.SetBinError  (ib, 0.0)

    return h_bg

# ----------------------------------------------------------------------
# compute a single scale factor that brings the TF1 integral in the
# side-bands into agreement with the *cut* histogram’s side-band counts
# ----------------------------------------------------------------------
def sideband_scale(sb_ranges, tf1, h_cut):
    """Return multiplicative scale factor for *tf1*."""
    bw   = h_cut.GetBinWidth(1)           # constant in your spectra
    obs  = 0.0
    pred = 0.0
    for lo, hi in sb_ranges:
        # observed counts (already per-bin, no bin-width factor)
        obs  += h_cut.Integral(h_cut.FindBin(lo), h_cut.FindBin(hi))
        # predicted counts from the fit (density → counts)
        pred += tf1.Integral(lo, hi) / bw
    return 1.0 if pred <= 0.0 else obs / pred


# ----------------------------------------------------------------------
# convenience: apply that scale factor to *all* TF1 parameters
# ----------------------------------------------------------------------
def rescale_tf1(tf1, scale):
    """Multiply every parameter of *tf1* by *scale* (in-place)."""
    for ip in range(tf1.GetNpar()):
        tf1.SetParameter(ip, tf1.GetParameter(ip) * scale)

################################################################################################################################################

#no_bg_subtract=True
no_bg_subtract=False

def bg_fit(phi_setting, inpDict, hist):

    if no_bg_subtract:
        fit_func = TF1("fit_func_zero", "0", inpDict["mm_min"], inpDict["mm_max"])
        fit_vis  = fit_func.Clone(f"{hist.GetName()}_bg_vis")
        bg_par   = 0
        return fit_func, fit_vis, bg_par    
    
    mm_min = inpDict["mm_min"]
    mm_max = inpDict["mm_max"]

    # --- Use physics-motivated wide sidebands for fitting, NOT just mm_min/mm_max ---
    sb_left = inpDict.get("sb_left", (1.070, 1.095))
    sb_right = inpDict.get("sb_right", (1.135, 1.180))
    sig_lo = max(1.107, mm_min)
    sig_hi = min(1.123, mm_max)

    print(f"Signal window: [{sig_lo:.3f}, {sig_hi:.3f}]")
    print(f"Sideband left: [{sb_left[0]:.3f}, {sb_left[1]:.3f}]")
    print(f"Sideband right: [{sb_right[0]:.3f}, {sb_right[1]:.3f}]")

    h_sb = hist.Clone(hist.GetName() + "_sb")

    def is_valid_sideband(sb):
        return sb[1] > sb[0]

    n_sb_bins = 0
    for ib in range(1, h_sb.GetNbinsX() + 1):
        x = h_sb.GetBinCenter(ib)
        in_sb_left  = is_valid_sideband(sb_left)  and (sb_left[0]  <= x <= sb_left[1])
        in_sb_right = is_valid_sideband(sb_right) and (sb_right[0] <= x <= sb_right[1])
        if not (in_sb_left or in_sb_right):
            h_sb.SetBinContent(ib, 0.0)
            h_sb.SetBinError(ib, 0.0)
        else:
            n_sb_bins += 1

    print("Nonzero bins in h_sb after masking:", n_sb_bins)

    if n_sb_bins < 2:
        print("[WARNING] All sideband bins are empty after masking! Fit will be flat zero.")
        fit_func = TF1("fit_func_zero", "0", mm_min, mm_max)
        fit_vis  = fit_func.Clone(f"{hist.GetName()}_bg_vis")
        bg_par   = 0
        fit_hist_inrange = get_fit_histogram_in_range(fit_func, hist, mm_min, mm_max)
        return fit_hist_inrange, fit_vis, bg_par

    fit_min = sb_left[0]
    fit_max = sb_right[1]
    fit_func = TF1("fit_func", "pol1", fit_min, fit_max)
    h_sb.Fit(fit_func, "Q0")

    # ---------------------------------------------------------------
    # renormalise the already-fitted function to the *cut* histogram
    # ---------------------------------------------------------------
    sb_list = [sb_left, sb_right]
    scale   = sideband_scale(sb_list, fit_func, hist)   # <-- hist is the *cut* one
    rescale_tf1(fit_func, scale)

    bg_par = fit_func.Integral(sig_lo, sig_hi) / hist.GetBinWidth(1)
    bg_par = max(0.0, bg_par)  # physical prior
    fit_vis = fit_func.Clone(f"{hist.GetName()}_bg_vis")

    scale = inpDict.get("bg_stat_scale", 1.0)
    if scale != 1.0:
        for ip in range(fit_func.GetNpar()):
            fit_func.SetParameter(ip, fit_func.GetParameter(ip) * scale)
        bg_par *= scale

    fit_hist_inrange = get_fit_histogram_padded(fit_func, hist, mm_min, mm_max, n_pad=0)
    return fit_hist_inrange, fit_vis, bg_par
