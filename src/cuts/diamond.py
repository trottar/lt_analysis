#! /usr/bin/python
###########################################################################################################################
# Created - 26/July/22, Author - Jacob Murphy
# Based on script created - 20/July/21, Author - Muhammad Junaid (mjo147@uregina.ca), University of Regina, Canada (Copyright (c) junaid) #
###########################################################################################################################
# Python version of the pion plotting script. Now utilises uproot to select event of each type and writes them to a root file.
# Python should allow for easier reading of databases storing diferent variables.
# This version of script is for shift workers at JLab
# To run this script, execute: python3 scriptname runnumber

###################################################################################################################################################

# Import relevant packages
import ROOT
import numpy as np
import sys, math, os, subprocess
import array
import re # Regexp package - for string manipulation
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine, TCutG
from ROOT import TExec
from ROOT import kBlack, kBlue, kRed
from array import array
import pandas as pd
from scipy.optimize import curve_fit
import glob

###############################################################################################################################################
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

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import open_root_file, apply_bin_threshold

################################################################################################################################################
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
#################################################################################################################################################

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

def diamond_fit_legacy(Q2vsW_hist, Q2Val, fitrange=10):
    """Legacy (curve_fit/projection) implementation kept for reference."""

    def line(x, a, b):
        return a * x + b

    if fitrange <= 1:
        print("Fit range too small, setting to default of 10.")
        fitrange = 10

    lol, hil, lor, hir = [], [], [], []
    xvl, xvr = [], []

    fitl = Q2vsW_hist.GetXaxis().FindBin(Q2Val) - int(fitrange * 5)
    fitr = Q2vsW_hist.GetXaxis().FindBin(Q2Val) + int(fitrange)

    print("\nFinding diamond fits (legacy)...")
    for b in range(fitrange * 2):
        # Left side
        proj_l = Q2vsW_hist.ProjectionY("y", b + fitl, b + fitl)
        minYl = proj_l.GetBinCenter(proj_l.FindFirstBinAbove(0))
        maxYl = proj_l.GetBinCenter(proj_l.FindLastBinAbove(0))
        lol.append(minYl)
        hil.append(maxYl)

        # Right side
        proj_r = Q2vsW_hist.ProjectionY("y", b + fitr, b + fitr + 1)
        minYr = proj_r.GetBinCenter(proj_r.FindFirstBinAbove(0))
        maxYr = proj_r.GetBinCenter(proj_r.FindLastBinAbove(0))
        lor.append(minYr)
        hir.append(maxYr)

        xl = Q2vsW_hist.GetXaxis().GetBinCenter(b + fitl)
        xr = Q2vsW_hist.GetXaxis().GetBinCenter(b + fitr)
        xvl.append(xl)
        xvr.append(xr)

    # Fit all four sides
    popt_l_low, _ = curve_fit(line, xvl, lol)
    popt_l_high, _ = curve_fit(line, xvl, hil)
    popt_r_low, _ = curve_fit(line, xvr, lor)
    popt_r_high, _ = curve_fit(line, xvr, hir)

    return popt_l_low, popt_l_high, popt_r_low, popt_r_high


def _th2_to_numpy(th2):
    """Convert ROOT.TH2* into (histo[xbin,ybin], x_centers, y_centers)."""
    nbx = th2.GetNbinsX()
    nby = th2.GetNbinsY()

    x = np.array([th2.GetXaxis().GetBinCenter(i) for i in range(1, nbx + 1)], dtype=float)
    y = np.array([th2.GetYaxis().GetBinCenter(j) for j in range(1, nby + 1)], dtype=float)

    histo = np.empty((nbx, nby), dtype=float)
    for ix in range(1, nbx + 1):
        for iy in range(1, nby + 1):
            histo[ix - 1, iy - 1] = th2.GetBinContent(ix, iy)

    return histo, x, y


def _diamond_fit_numpy(
    histo,
    x,
    y,
    control_left=None,
    control_right=None,
    threshold=0.0,
):
    """
    Ported from diamond_new.py:
      - Threshold bins
      - Extract upper/lower W edges vs Q2 and left/right edges vs Q2
      - Fit each edge with a straight line: W = a*Q2 + b
    Returns: (down, up, left, right) where each is (a, b).
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    histo = np.asarray(histo, dtype=float)

    if histo.ndim != 2:
        raise ValueError("_diamond_fit_numpy: histo must be 2D")
    if len(x) != histo.shape[0] or len(y) != histo.shape[1]:
        raise ValueError("_diamond_fit_numpy: x/y sizes do not match histo shape")

    if control_left is None or control_right is None:
        # Use occupancy-weighted center and a small control window.
        mask = np.where(histo > threshold, 1.0, 0.0)
        if mask.sum() == 0:
            raise ValueError("_diamond_fit_numpy: histogram has no bins above threshold.")
        dx = x[1] - x[0] if len(x) > 1 else 1.0
        dy = y[1] - y[0] if len(y) > 1 else 1.0
        xc = np.average(x, weights=np.sum(mask, axis=1))
        yc = np.average(y, weights=np.sum(mask, axis=0))
        control_left = (xc - dx * 10.0, yc - dy * 10.0)
        control_right = (xc + dx * 10.0, yc + dy * 10.0)

    histo = np.where(histo > threshold, histo, 0.0)

    # Upper/Lower edges as W(Q2)
    start_ud = int(np.digitize(control_left[0], x, right=True))
    end_ud = int(np.digitize(control_right[0], x, right=True)) + 1
    start_ud = max(0, min(start_ud, len(x)))
    end_ud = max(0, min(end_ud, len(x)))

    x_u, up = [], []
    x_d, down = [], []

    for b in range(start_ud, end_ud):
        proj_ud = histo[b, :]
        nz = proj_ud > 0
        if not np.any(nz):
            continue
        x_u.append(x[b])
        up.append(y[nz][-1])
        x_d.append(x[b])
        down.append(y[nz][0])

    for b in range(end_ud, len(x)):
        proj_ud = histo[b, :]
        nz = proj_ud > 0
        if not np.any(nz):
            continue
        x_d.append(x[b])
        down.append(y[nz][0])

    for b in range(0, start_ud):
        proj_ud = histo[b, :]
        nz = proj_ud > 0
        if not np.any(nz):
            continue
        x_u.append(x[b])
        up.append(y[nz][-1])

    # Left/Right edges as W(Q2) by first finding Q2(W), then fitting W(Q2)
    start_lr = int(np.digitize(control_left[1], y, right=True))
    end_lr = int(np.digitize(control_right[1], y, right=True)) + 1
    start_lr = max(0, min(start_lr, len(y)))
    end_lr = max(0, min(end_lr, len(y)))

    y_l, left = [], []
    y_r, right = [], []

    for b in range(start_lr, end_lr):
        proj_lr = histo[:, b]
        nz = proj_lr > 0
        if not np.any(nz):
            continue
        y_l.append(y[b])
        left.append(x[nz][0])
        y_r.append(y[b])
        right.append(x[nz][-1])

    for b in range(end_lr, len(y)):
        proj_lr = histo[:, b]
        nz = proj_lr > 0
        if not np.any(nz):
            continue
        y_l.append(y[b])
        left.append(x[nz][0])

    for b in range(0, start_lr):
        proj_lr = histo[:, b]
        nz = proj_lr > 0
        if not np.any(nz):
            continue
        y_r.append(y[b])
        right.append(x[nz][-1])

    if len(down) < 2 or len(up) < 2 or len(left) < 2 or len(right) < 2:
        raise ValueError("Not enough edge points to perform the fit (increase stats or lower threshold).")

    popt_down = np.polyfit(x_d, down, 1)
    popt_up = np.polyfit(x_u, up, 1)
    popt_left = np.polyfit(left, y_l, 1)
    popt_right = np.polyfit(right, y_r, 1)

    return popt_down, popt_up, popt_left, popt_right


def diamond_fit(Q2vsW_hist, Q2Val, fitrange=10, threshold=0.0, use_legacy=False, auto_control=False):
    """
    Updated diamond-fit procedure (ported from diamond_new.py):
      - Fits *four* lines in (Q2, W): down/up/left/right in the form W = a*Q2 + b.

    Parameters
    ----------
    Q2vsW_hist : ROOT.TH2*
        Histogram with x=Q2, y=W.
    Q2Val : float
        Central Q2 value (used to define an x-control window matching the legacy behavior).
    fitrange : int
        Defines the x-control window. Legacy used: [Q2bin - 5*fitrange, Q2bin + fitrange].
    threshold : float
        Bin-content threshold used to define occupied bins.
        Use threshold=0.0 if you already applied apply_bin_threshold().
    use_legacy : bool
        If True, uses the original projection+curve_fit implementation.

    Returns
    -------
    (down, up, left, right) : tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]
        Each is (a, b) for W = a*Q2 + b.
    """

    if use_legacy:
        return diamond_fit_legacy(Q2vsW_hist, Q2Val, fitrange)

    if fitrange <= 1:
        print("Fit range too small, setting to default of 10.")
        fitrange = 10

    histo, x, y = _th2_to_numpy(Q2vsW_hist)

    mask = np.where(histo > threshold, 1.0, 0.0)
    if mask.sum() == 0:
        raise ValueError("diamond_fit: histogram has no bins above threshold.")

    # Compute a y-center from occupancy.
    dy = y[1] - y[0] if len(y) > 1 else 1.0
    yc = np.average(y, weights=np.sum(mask, axis=0))

    control_left = None
    control_right = None

    # Legacy-style control window (centered on Q2Val) if requested.
    if not auto_control:
        nbx = len(x)
        q2bin = Q2vsW_hist.GetXaxis().FindBin(Q2Val)
        fitl_bin = max(1, q2bin - int(fitrange * 5))
        fitr_bin = min(nbx, q2bin + int(fitrange))

        control_left = (x[fitl_bin - 1], yc - dy * 10.0)
        control_right = (x[fitr_bin - 1], yc + dy * 10.0)

    try:
        return _diamond_fit_numpy(
            histo,
            x,
            y,
            control_left=control_left,
            control_right=control_right,
            threshold=threshold,
        )
    except Exception as e:
        print(f"diamond_fit: new procedure failed ({e}); falling back to legacy.")
        return diamond_fit_legacy(Q2vsW_hist, Q2Val, fitrange)

###################################################################################################
# Convex polygon helpers (used for "common overlap" cuts)

def _sort_ccw_points(points):
    pts = [(float(p[0]), float(p[1])) for p in points]
    if len(pts) < 3:
        return pts
    cx = sum(p[0] for p in pts) / float(len(pts))
    cy = sum(p[1] for p in pts) / float(len(pts))
    pts = sorted(pts, key=lambda p: math.atan2(p[1] - cy, p[0] - cx))
    return pts

def _dedup_points(points, tol=1e-6):
    out = []
    for p in points:
        if not out:
            out.append(p)
            continue
        if abs(p[0] - out[-1][0]) > tol or abs(p[1] - out[-1][1]) > tol:
            out.append(p)
    if len(out) >= 2 and abs(out[0][0] - out[-1][0]) <= tol and abs(out[0][1] - out[-1][1]) <= tol:
        out.pop()
    return out

def _segment_line_intersection(p1, p2, q1, q2):
    x1, y1 = float(p1[0]), float(p1[1])
    x2, y2 = float(p2[0]), float(p2[1])
    x3, y3 = float(q1[0]), float(q1[1])
    x4, y4 = float(q2[0]), float(q2[1])

    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(den) < 1.0e-12:
        return None

    px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / den
    py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / den
    if (not np.isfinite(px)) or (not np.isfinite(py)):
        return None
    return (px, py)

def _inside_halfplane(p, a, b, eps=1e-12):
    return ((b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])) >= -eps

def _convex_polygon_intersection(subject, clip):
    if subject is None or clip is None:
        return None
    subj = _sort_ccw_points(subject)
    clip = _sort_ccw_points(clip)
    output = subj

    for i in range(len(clip)):
        if not output:
            return None
        cp1 = clip[i]
        cp2 = clip[(i + 1) % len(clip)]

        input_list = output
        output = []
        S = input_list[-1]

        for E in input_list:
            if _inside_halfplane(E, cp1, cp2):
                if not _inside_halfplane(S, cp1, cp2):
                    inter = _segment_line_intersection(S, E, cp1, cp2)
                    if inter is not None:
                        output.append(inter)
                output.append(E)
            elif _inside_halfplane(S, cp1, cp2):
                inter = _segment_line_intersection(S, E, cp1, cp2)
                if inter is not None:
                    output.append(inter)
            S = E

        output = _dedup_points(output)

    if len(output) < 3:
        return None
    return _sort_ccw_points(output)

def _point_in_convex_poly(x, y, poly, eps=1e-12):
    if poly is None or len(poly) < 3:
        return False
    px, py = float(x), float(y)
    pts = _sort_ccw_points(poly)
    for i in range(len(pts)):
        a = pts[i]
        b = pts[(i + 1) % len(pts)]
        if ((b[0] - a[0]) * (py - a[1]) - (b[1] - a[1]) * (px - a[0])) < -eps:
            return False
    return True

def _line_intersection_ab(p1, p2):
    a1, b1 = float(p1[0]), float(p1[1])
    a2, b2 = float(p2[0]), float(p2[1])
    den = (a1 - a2)
    if abs(den) < 1.0e-12:
        return None
    q2 = (b2 - b1) / den
    w = a1 * q2 + b1
    if (not np.isfinite(q2)) or (not np.isfinite(w)):
        return None
    return (q2, w)

def _poly_from_diamond_fits(fits):
    down, up, left, right = fits
    corners = [
        _line_intersection_ab(down, left),
        _line_intersection_ab(down, right),
        _line_intersection_ab(up, right),
        _line_intersection_ab(up, left),
    ]
    pts = [p for p in corners if p is not None]
    if len(pts) < 3:
        return None
    return _sort_ccw_points(pts)


def _parse_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    v = str(value).strip().lower()
    if v in ("1", "true", "t", "yes", "y", "on"):
        return True
    if v in ("0", "false", "f", "no", "n", "off", ""):
        return False
    return default


HARDCODED_DIAMOND_AB_PARAMS = {
    "a1": -0.22478443202783813,
    "b1": 3.6863942419482063,
    "a2": -0.1780992912706183,
    "b2": 3.586856595515984,
    "a3": -0.3516380890076293,
    "b3": 4.210272623820992,
    "a4": -0.3664439092360068,
    "b4": 4.427889346144653,
}


def _hardcoded_diamond_fit_tuple():
    # Order must match: (down, up, left, right)
    return (
        np.array([HARDCODED_DIAMOND_AB_PARAMS["a1"], HARDCODED_DIAMOND_AB_PARAMS["b1"]], dtype=float),
        np.array([HARDCODED_DIAMOND_AB_PARAMS["a2"], HARDCODED_DIAMOND_AB_PARAMS["b2"]], dtype=float),
        np.array([HARDCODED_DIAMOND_AB_PARAMS["a3"], HARDCODED_DIAMOND_AB_PARAMS["b3"]], dtype=float),
        np.array([HARDCODED_DIAMOND_AB_PARAMS["a4"], HARDCODED_DIAMOND_AB_PARAMS["b4"]], dtype=float),
    )

def _tgraph_from_poly(poly, name):
    if poly is None or len(poly) < 3:
        return None
    pts = _sort_ccw_points(poly)
    g = ROOT.TGraph(len(pts) + 1)
    g.SetName(name)
    for i, (qx, wy) in enumerate(pts):
        g.SetPoint(i, qx, wy)
    g.SetPoint(len(pts), pts[0][0], pts[0][1])
    g.SetFillStyle(0)
    return g

###################################################################################################

def DiamondPlot(ParticleType, Q2Val, Q2min, Q2max, WVal, Wmin, Wmax, phi_setting, tmin, tmax, inpDict):

    inpDict["use_hardcoded_diamond_fits"] = True

    ##############
    # HARD CODED #
    ##############
    if Q2Val == 0.4:
        Q2Val = 0.385
        WVal = 2.20
        Qs = "0p4"
        Ws = "2p20"
    elif Q2Val == "2p1":
        Q2Val = 2.115
        WVal = 2.95
        Qs = "2p1"
        Ws = "2p95"
    else:
        Qs = str(Q2Val).replace('.','p')
        Ws = str(WVal).replace('.','p')
    ##############
    ##############
    ##############

    use_hardcoded_diamond_fits = _parse_bool(inpDict.get("use_hardcoded_diamond_fits", False), default=False)
    hardcoded_fit_results = _hardcoded_diamond_fit_tuple() if use_hardcoded_diamond_fits else None
    hardcoded_cut_poly = _poly_from_diamond_fits(hardcoded_fit_results) if hardcoded_fit_results is not None else None
    if use_hardcoded_diamond_fits:
        print("INFO: Using hardcoded diamond fit coefficients (a1..b4).")
        if hardcoded_cut_poly is None:
            print("!!!!! ERROR !!!!!\n Hardcoded diamond coefficients do not form a valid polygon.\n!!!!! ERROR !!!!!")
            sys.exit(2)
    
    FilenameOverride = 'Q'+Qs+'W'+Ws
    
    Analysis_Distributions = OUTPATH+"/{}_{}_diamond_{}.pdf".format(phi_setting, ParticleType, FilenameOverride)

    nbins = 200
    
    lowe_input = False
    mide_input = False
    highe_input = False

    # Arbitrary lengths far longer than any file name
    lenh = 10000
    lenm = 10000
    lenl = 10000
    if(phi_setting == '0'): phi_setting = ""
    print("\n\nKinematics: ",FilenameOverride,"\nPhi Setting: ",phi_setting)
    for file in glob.glob(OUTPATH+'/*'+phi_setting+'*'+ParticleType+'*'+FilenameOverride+'*.root'):
	# Searches through OUTPUT recursively for files matching the wild card format, taking the shortest one
        # Shortest file assumed to be full analyisis as it will not have "part" or "week" or "dummy" labels
        #print(file)
        if "Simc" in file:
            continue
        if "high" in file:
            if (len(file) < lenh):
                highe_input = file
                lenh = len(file)
        
        if "mid" in file:
            if (len(file) < lenm):
                mide_input = file
                lenm = len(file)

        if "low" in file:
            if (len(file) < lenl):
                lowe_input = file
                lenl = len(file)

    if (highe_input == False and mide_input == False and lowe_input == False):
        print("!!!!! ERROR !!!!!\n No valid file found! \n!!!!! ERROR !!!!!")
        sys.exit(1)

    ##############################################################################################################################################
    # Common-cut polygon:
    # Intersect all available epsilon/phi diamond polygons so the chosen cut is
    # guaranteed to lie inside the full set of setting diamonds.
    common_poly = None
    common_poly_sources = []

    # Use the same hard-coded threshold rule as the main hist thresholding.
    if Q2Val == 3.0 and WVal == 3.14:
        event_threshold_common = 15
    else:
        event_threshold_common = 5

    if (phi_setting == "Center") and (not use_hardcoded_diamond_fits):
        def _find_shortest_root_any(phi_tokens, eps_tag):
            best = None
            best_len = 10**9
            for tok in phi_tokens:
                for f in glob.glob(OUTPATH + '/*' + tok + '*' + ParticleType + '*' + FilenameOverride + '*.root'):
                    fl = f.lower()
                    if eps_tag in fl:
                        if len(f) < best_len:
                            best = f
                            best_len = len(f)
            return best

        def _find_center_fallback(eps_tag):
            cand = None
            cand_len = 10**9
            for f in glob.glob(OUTPATH + '/*' + ParticleType + '*' + FilenameOverride + '*.root'):
                fl = f.lower()
                if (eps_tag in fl) and ("left" not in fl) and ("right" not in fl):
                    if len(f) < cand_len:
                        cand = f
                        cand_len = len(f)
            return cand

        def _build_q2w_hist_inmem(root_path, hname):
            if root_path is None:
                return None
            infile = open_root_file(root_path, "READ")
            if (not infile) or infile.IsZombie():
                try:
                    if infile:
                        infile.Close()
                except Exception:
                    pass
                return None
            tree = infile.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))
            if not tree:
                infile.Close()
                return None
            ROOT.gROOT.cd()
            h = TH2D(hname, "", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
            h.SetDirectory(0)
            for ev in tree:
                h.Fill(ev.Q2, ev.W)
            apply_bin_threshold(h, event_threshold_common)
            infile.Close()
            return h

        phi_defs_common = [
            ("Center", ["Center", "center"]),
            ("Left",   ["Left", "left"]),
            ("Right",  ["Right", "right"]),
        ]

        eps_defs_common = [
            ("Low", "low"),
            ("Mid", "mid"),
            ("High", "high"),
        ]

        def _find_source_file(eps_tag, phi_lbl, toks):
            # For center, prefer the file already selected earlier in this call.
            if phi_lbl == "Center":
                if eps_tag == "low" and lowe_input:
                    return lowe_input
                if eps_tag == "mid" and mide_input:
                    return mide_input
                if eps_tag == "high" and highe_input:
                    return highe_input
            fpath = _find_shortest_root_any(toks, eps_tag)
            if fpath is None and phi_lbl == "Center":
                fpath = _find_center_fallback(eps_tag)
            return fpath

        polys = []
        seen_sources = set()
        for eps_lbl, eps_tag in eps_defs_common:
            for phi_lbl, toks in phi_defs_common:
                src_file = _find_source_file(eps_tag, phi_lbl, toks)
                if not src_file:
                    continue
                src_key = (eps_lbl, phi_lbl, os.path.normpath(src_file))
                if src_key in seen_sources:
                    continue
                seen_sources.add(src_key)

                htmp = _build_q2w_hist_inmem(
                    src_file,
                    "Q2vsW_{}_{}_common_{}".format(eps_tag, phi_lbl.lower(), FilenameOverride),
                )
                if not htmp:
                    continue
                try:
                    fits_tmp = diamond_fit(htmp, Q2Val, fitrange=10, threshold=0.0, auto_control=True)
                    poly_tmp = _poly_from_diamond_fits(fits_tmp)
                    if poly_tmp:
                        polys.append(poly_tmp)
                        common_poly_sources.append((eps_lbl, phi_lbl, src_file))
                except Exception as e:
                    print("WARNING: common-cut: failed {}-{} fit ({})".format(eps_lbl, phi_lbl, e))

        if len(polys) >= 1:
            common_poly = polys[0]
            for p in polys[1:]:
                common_poly = _convex_polygon_intersection(common_poly, p)
                if common_poly is None:
                    break
            if common_poly is not None:
                print("Common-cut polygon built with {} vertices from {} sources".format(len(common_poly), len(common_poly_sources)))
            else:
                print("WARNING: common-cut: no overlap polygon found across {} sources".format(len(common_poly_sources)))

    labelh = ""
    labelm = ""
    labell = ""
    if (highe_input !=False):
        #print("test high")
        labelh = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{High (Blue) Epsilon}; Q2; W"
        if (mide_input !=False):
            labelh = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{High (Blue) and Mid (Red) Epsilon}; Q2; W"
            #print("test high mid")
            if (lowe_input !=False):
                labelh = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{High (Blue), Mid (Red), and Low (Green) Epsilon}; Q2; W"
                #print("test high mid low")
        elif (lowe_input !=False):
            #print("test high low")
            labelh = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{High (Blue) and Low (Red) Epsilon}; Q2; W"
    elif (mide_input !=False):
        #print("test mid")
        labelm = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{Mid (Blue) Epsilon}; Q2; W"
        if (lowe_input !=False):
            #print("test mid low")
            labelm = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{Mid (Blue) and Low (Red) Epsilon}; Q2; W"
    elif (lowe_input !=False):
        #print("test low")
        labell = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{Low (Blue) Epsilon}; Q2; W"


    Title = ""
    Q2vsW_cut = TH2D("Q2vsW_cut", labelh, nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_mide_cut = TH2D("Q2vsW_mide_cut",labelm, nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_lowe_cut = TH2D("Q2vsW_lowe_cut", labell, nbins, Q2min, Q2max, nbins, Wmin, Wmax)    

    Q2vsW_hi_cut = TH2D("Q2vsW_high_cut", "High Epsilon Q2 vs W Dist for Prompt Events (Prompt Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_mi_cut = TH2D("Q2vsW_middle_cut","Mid Epsilon Q2 vs W Dist for Prompt Events (Prompt Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_lo_cut = TH2D("Q2vsW_low_cut", "Low Epsilon Q2 vs W Dist for Prompt Events (Prompt Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)

    W_cut = TH1D("W_cut", "High Epsilon W Dist for Prompt Events (Prompt Cut); W", nbins, Wmin, Wmax)
    Q2_cut = TH1D("Q2_cut", "High Epsilon Q2 Dist for Prompt Events (Prompt  Cut); Q2", nbins, Q2min, Q2max)
    t_cut = TH1D("t_cut", "High Epsilon -t Dist for Prompt Events (t-Range  Cut); -t", nbins, tmin, tmax)
    t_mi_cut = TH1D("t_mi_cut", "Mid Epsilon -t Dist for Prompt Events (t-Range  Cut); -t", nbins, tmin, tmax)

    Q2vsW_lolo_cut = TH2D("Q2vsW_low_lowcut", "Low Epsilon Q2 vs W Dist for Prompt Events (Diamond Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_hilo_cut = TH2D("Q2vsW_high_lowcut", "High Epsilon Q2 vs W Dist for Prompt Events (Diamond and t Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_milo_cut = TH2D("Q2vsW_mid_lowcut","Mid Epsilon Q2 vs W Dist for Prompt Events (Diamond and t Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_himi_cut = TH2D("Q2vsW_high_midcut", "High Epsilon Q2 vs W Dist for Prompt Events (Mid-Diamond and t Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)

    if hardcoded_cut_poly is not None:
        cut_poly = _sort_ccw_points(hardcoded_cut_poly)
    else:
        cut_poly = _sort_ccw_points(common_poly) if common_poly is not None else None

    def _passes_active_cut(q2, w):
        if cut_poly is None:
            return False
        return _point_in_convex_poly(q2, w, cut_poly)

    for k in range(0, 3):

	# Construct the name of the rootfile based upon the info we provided
        if (k==2):
            if (highe_input == False): 
                continue
            elif (highe_input != False): # Special condition, with 5th arg, use 5th arg as file name
                rootName = highe_input
                print("\n***** High Epsilon File Found! *****\n")
        if (k==1):
            if (mide_input == False): 
                continue
            elif (mide_input != False): # Special condition, with 5th arg, use 5th arg as file name
                rootName = mide_input
                print("\n*****Mid Epsilon File Found! *****\n")
        if (k==0):
            if (lowe_input == False): 
                continue
            elif (lowe_input != False): # Special condition, with 5th arg, use 5th arg as file name
                rootName = lowe_input
                print("\n***** Low Epsilon File Found! *****\n")
        print ("Attempting to process %s" %(rootName))

	###############################################################################################################################################
        ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen

	###############################################################################################################################################

	# Read stuff from the main event tree
        infile = open_root_file(rootName, "READ")

	# Assumes 2021 trees do not have Prompt MM cut, as some do not right now. *** NEED TO BE REPLAYED AGAIN WITH THIS BRANCH ***
        Cut_Events_all_noRF_tree = infile.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))

	##############################################################################################################################################
        countB = 0
        countA = 0
        badfit = True
        if (k==2): # High
            # for event in Cut_Events_Prompt_tree:
            for event in Cut_Events_all_noRF_tree:
                Q2vsW_cut.Fill(event.Q2, event.W)
                Q2vsW_hi_cut.Fill(event.Q2, event.W)
        elif (k==1): # Mid
            # for event in Cut_Events_Prompt_tree:
            for event in Cut_Events_all_noRF_tree:
                Q2vsW_mide_cut.Fill(event.Q2, event.W)
                Q2vsW_mi_cut.Fill(event.Q2, event.W)
        elif (k==0): # Low
            # for event in Cut_Events_Prompt_tree:
            for event in Cut_Events_all_noRF_tree:
                Q2vsW_lowe_cut.Fill(event.Q2, event.W)
                Q2vsW_lo_cut.Fill(event.Q2, event.W)
                W_cut.Fill(event.W)
                Q2_cut.Fill(event.Q2)
                countB +=1

        # Apply the threshold to the histograms
        ##############
        # HARD CODED #
        ##############
        if Q2Val == 3.0 and WVal == 3.14:
            event_threshold = 15 # Q2=3.0,W=3.14
        else:
            event_threshold = 5 # Q2=2.115,W=2.95 | Q2=3.0,W=2.32 | Q2=4.4,W=2.74
        ##############
        ##############
        ##############
        if (k==2): # High
            apply_bin_threshold(Q2vsW_cut, event_threshold)
            apply_bin_threshold(Q2vsW_hi_cut, event_threshold)
        elif (k==1): # Mid
            apply_bin_threshold(Q2vsW_mide_cut, event_threshold)
            apply_bin_threshold(Q2vsW_mi_cut, event_threshold)
        elif (k==0): # Low
            apply_bin_threshold(Q2vsW_lowe_cut, event_threshold)
            apply_bin_threshold(Q2vsW_lo_cut, event_threshold)
                
        # Does assume nbins bins for Q2 and W, centered at kinematic values
        minQ = Q2_cut.FindFirstBinAbove(0)
        maxQ = Q2_cut.FindLastBinAbove(0)
        fitrange = int((maxQ-minQ)/100)
        print("fitrange: ",fitrange)
        if (k == 0):  # Low epsilon
            # Build low-epsilon polygon fallback if common overlap polygon is unavailable.
            if hardcoded_fit_results is not None:
                fit_results = hardcoded_fit_results
            else:
                fit_results = diamond_fit(Q2vsW_lowe_cut, Q2Val, fitrange)
            low_fit_poly = _poly_from_diamond_fits(fit_results)
            if (cut_poly is None) and (low_fit_poly is not None) and (phi_setting != "Center"):
                cut_poly = _sort_ccw_points(low_fit_poly)

            for event in Cut_Events_all_noRF_tree:
                if _passes_active_cut(event.Q2, event.W):
                    Q2vsW_lolo_cut.Fill(event.Q2, event.W)
                    countA += 1
                    
        if (lowe_input != False and k>0):
            print("\n\n")
            if (k==2):
                for event in Cut_Events_all_noRF_tree:
                    if _passes_active_cut(event.Q2, event.W):
                        if (tmax != False):
                            Q2vsW_hilo_cut.Fill(event.Q2, event.W)
                            t_cut.Fill(-event.MandelT)
                        else:
                            print("!!!!! Error! tmax not found! Skipping t-range cut !!!!!")
                            Q2vsW_hilo_cut.Fill(event.Q2, event.W)
            elif (k==1):
                for event in Cut_Events_all_noRF_tree:
                    if _passes_active_cut(event.Q2, event.W):
                        if (tmax != False):
                            Q2vsW_milo_cut.Fill(event.Q2, event.W)
                            t_mi_cut.Fill(-event.MandelT)
                        else:
                            print("!!!!! Error! tmax not found! Skipping t-range cut !!!!!")
                            Q2vsW_milo_cut.Fill(event.Q2, event.W)        

        print("Histograms filled")

        infile.Close()

    if phi_setting == "Center":
        if cut_poly is None:
            print("!!!!! ERROR !!!!!\n No valid cut polygon available.\n!!!!! ERROR !!!!!")
            sys.exit(2)

        paramDict = {
            "cut_mode": "poly",
            "poly_points": [[float(p[0]), float(p[1])] for p in _sort_ccw_points(cut_poly)],
        }
        if use_hardcoded_diamond_fits:
            paramDict["diamond_fit_source"] = "hardcoded"
            paramDict["diamond_fit_params"] = dict(HARDCODED_DIAMOND_AB_PARAMS)
            paramDict["diamond_fit_lines"] = {
                "down": [float(hardcoded_fit_results[0][0]), float(hardcoded_fit_results[0][1])],
                "up": [float(hardcoded_fit_results[1][0]), float(hardcoded_fit_results[1][1])],
                "left": [float(hardcoded_fit_results[2][0]), float(hardcoded_fit_results[2][1])],
                "right": [float(hardcoded_fit_results[3][0]), float(hardcoded_fit_results[3][1])],
            }
        if common_poly is not None:
            paramDict["poly_sources"] = [
                {"epsilon": eps, "phi": phi, "file": fpath} for (eps, phi, fpath) in common_poly_sources
            ]

    else:

        paramDict = {}

    ##############################################################################################################################################
    c1_kin = TCanvas("c1_kin", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
    gStyle.SetTitleFontSize(0.03)
    gStyle.SetPalette(86)
    ex1 = TExec("ex1","gStyle->SetPalette(86)")
    ex2 = TExec("ex2","gStyle->SetPalette(75)")
    ex3 = TExec("ex3","gStyle->SetPalette(68)")
    gStyle.SetOptStat(0)
    pages = 2
    if (highe_input !=False):
        #print("test high")
        Q2vsW_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
        Q2vsW_cut.Draw("col")
        ex1.Draw()
        Q2vsW_cut.Draw("col same")
        if (mide_input !=False):
            #print("test high mid")
            Q2vsW_mide_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_mide_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
            ex2.Draw()
            Q2vsW_mide_cut.Draw("col same")
            pages = 3
            if (lowe_input !=False):
                #print("test high mid low")
                Q2vsW_lowe_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
                Q2vsW_lowe_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)                
                ex3.Draw()
                Q2vsW_lowe_cut.Draw("col same")
                pages = 6
        elif (lowe_input !=False):
            #print("test high low")
            Q2vsW_lowe_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_lowe_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
            ex2.Draw()
            Q2vsW_lowe_cut.Draw("col same")
            pages = 4
    elif (mide_input !=False):
	#print("test mid")
        Q2vsW_mide_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_mide_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_mide_cut.Draw("colz")
        ex1.Draw()
        Q2vsW_mide_cut.Draw("col same")
        if (lowe_input !=False):
            #print("test mid low")
            Q2vsW_lowe_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_lowe_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)            
            ex2.Draw()
            Q2vsW_lowe_cut.Draw("col same")
            pages = 4
    elif (lowe_input !=False):
        #print("test low")
        Q2vsW_lowe_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_lowe_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_lowe_cut.Draw("colz")
        ex1.Draw()
        Q2vsW_lowe_cut.Draw("col same")

    c1_kin.Print(Analysis_Distributions + "(")

    #############################################################################################################################

    end = ""
    endm = ""
    endc = ""
    endf = ""
    if (pages==2): end = ")"
    if (pages==3): endm = ")"
    if (pages==4): endc = ")"
    if (pages==6): endf = ")"

    gStyle.SetOptStat(1)

    # Defer PDF closing to the final overlay plot
    end = ""
    endm = ""
    endc = ""
    endf = ""
    close_pdf = ")"

    gStyle.SetOptStat(1)
    gStyle.SetPalette(55)

    if (tmax != False):
        c1_kint = TCanvas("c1_kint", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        t_cut.Draw("colz")
        c1_kint.Print(Analysis_Distributions)
    if (highe_input != False):
        c1_kinh = TCanvas("c1_kinh", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        Q2vsW_hi_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_hi_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
        Q2vsW_hi_cut.Draw("colz")
        c1_kinh.Print(Analysis_Distributions+end)
    if (mide_input != False):
        c1_kinm = TCanvas("c1_kinm", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        Q2vsW_mi_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_mi_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_mi_cut.Draw("colz")
        c1_kinm.Print(Analysis_Distributions+end+endm)
    if (lowe_input != False):
        c1_kinl = TCanvas("c1_kinl", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        Q2vsW_lo_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_lo_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_lo_cut.Draw("colz")
        c1_kinl.Print(Analysis_Distributions)
        c1_kinll = TCanvas("c1_kinll", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        Q2vsW_lolo_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_lolo_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_lolo_cut.Draw("colz")
        #    lol.clear()
        #   lor.clear()
        #  hil.clear()
        # hir.clear()
        # xvl.clear()
        # xvr.clear()
        gcut_lolo = _tgraph_from_poly(cut_poly, "gCutLolo_{}".format(FilenameOverride))
        if gcut_lolo:
            gcut_lolo.SetLineColor(ROOT.kBlack)
            gcut_lolo.SetLineWidth(5)
            gcut_lolo.SetLineStyle(1)
            gcut_lolo.Draw("L SAME")
        c1_kinll.Print(Analysis_Distributions+end)

        if (mide_input != False):
            c1_kinml = TCanvas("c1_kinml", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
            Q2vsW_milo_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_milo_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)            
            Q2vsW_milo_cut.Draw("colz")
            c1_kinml.Print(Analysis_Distributions+endc)

        if (highe_input != False):
            c1_kinhl = TCanvas("c1_kinhl", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
            Q2vsW_hilo_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_hilo_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
            Q2vsW_hilo_cut.Draw("colz")
            c1_kinhl.Print(Analysis_Distributions+endc+endf)
	
            

    #############################################################################################################################
    # Overlay plot: all diamonds (Q2 vs W acceptance) superimposed, color-coded by epsilon
    #
    # Use contours at ~0.5 counts after thresholding so that bins with content > 0 define the diamond boundary.
    #
    gStyle.SetOptStat(0)

    

    # Overlay plot: all three phi settings (Center/Left/Right) for High and Low epsilon on one plot
    #
    # Draw ONLY the fitted diamond outlines (no distributions). Each (epsilon, phi) combination gets a unique color.
    #
    gStyle.SetOptStat(0)

    c1_overlay_allphi = TCanvas(
        "c1_overlay_allphi_{}".format(FilenameOverride),
        "%s All-Phi Diamond Outline Overlay" % ParticleType,
        100, 0, 1100, 950
    )

    frame_overlay_allphi = TH2D(
        "frame_overlay_allphi_{}".format(FilenameOverride),
        "All Diamond Outlines Overlay; Q2; W",
        nbins, Q2min, Q2max,
        nbins, Wmin, Wmax
    )
    frame_overlay_allphi.GetXaxis().SetRangeUser(Q2min - Q2min * 0.1, Q2max + Q2max * 0.1)
    frame_overlay_allphi.GetYaxis().SetRangeUser(Wmin - Wmin * 0.1, Wmax + Wmax * 0.1)
    frame_overlay_allphi.Draw()

    leg = TLegend(0.12, 0.70, 0.58, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    # Keep Python references to TGraph objects so ROOT legend pointers remain valid through Print()
    overlay_keepalive = []

    # Same threshold rule as the main plot block (used ONLY to build the hist for edge finding)
    if Q2Val == 3.0 and WVal == 3.14:
        event_threshold_overlay = 15
    else:
        event_threshold_overlay = 5

    def _find_shortest_root(phi_tokens, eps_tag):
        # Choose the shortest matching file path (same heuristic as main selection)
        best = None
        best_len = 10**9
        for tok in phi_tokens:
            for f in glob.glob(OUTPATH + '/*' + tok + '*' + ParticleType + '*' + FilenameOverride + '*.root'):
                fl = f.lower()
                if eps_tag in fl:
                    if len(f) < best_len:
                        best = f
                        best_len = len(f)
        return best

    def _build_q2w_hist_for_file(root_path, hname):
        # Build a Q2vsW TH2D in memory (SetDirectory(0)) so it survives infile.Close().
        if root_path is None:
            return None

        infile = open_root_file(root_path, "READ")
        if (not infile) or infile.IsZombie():
            print("WARNING: Could not open ROOT file: {}".format(root_path))
            try:
                if infile:
                    infile.Close()
            except Exception:
                pass
            return None

        tree = infile.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))
        if not tree:
            print("!!!!! ERROR !!!!!\n Missing tree in file: {}\n!!!!! ERROR !!!!!".format(root_path))
            infile.Close()
            return None

        ROOT.gROOT.cd()
        h = TH2D(hname, "", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
        h.SetDirectory(0)

        for ev in tree:
            h.Fill(ev.Q2, ev.W)

        apply_bin_threshold(h, event_threshold_overlay)
        infile.Close()
        return h

    def _line_intersection(p1, p2):
        # Each p is [a, b] for W = a*Q2 + b
        a1, b1 = float(p1[0]), float(p1[1])
        a2, b2 = float(p2[0]), float(p2[1])
        den = (a1 - a2)
        if abs(den) < 1.0e-12:
            return None
        q2 = (b2 - b1) / den
        w = a1 * q2 + b1
        if (not np.isfinite(q2)) or (not np.isfinite(w)):
            return None
        return (q2, w)

    def _diamond_graph_from_fits(fits, gname):
        # fits = (down, up, left, right), each [a, b] for W = a*Q2 + b
        down, up, left, right = fits

        corners = [
            _line_intersection(down, left),
            _line_intersection(down, right),
            _line_intersection(up, right),
            _line_intersection(up, left),
        ]
        pts = [p for p in corners if p is not None]
        if len(pts) < 3:
            return None

        # Order points around centroid
        cx = sum(p[0] for p in pts) / float(len(pts))
        cy = sum(p[1] for p in pts) / float(len(pts))
        pts = sorted(pts, key=lambda p: math.atan2(p[1] - cy, p[0] - cx))

        g = ROOT.TGraph(len(pts) + 1)
        g.SetName(gname)
        for i, (qx, wy) in enumerate(pts):
            g.SetPoint(i, qx, wy)
        g.SetPoint(len(pts), pts[0][0], pts[0][1])  # close
        g.SetFillStyle(0)
        g.SetLineStyle(1)
        g.SetLineWidth(3)
        return g

    def _draw_diamond_outline(root_path, eps_label, phi_label, color):
        if root_path is None:
            print("WARNING: No {} epsilon file found for phi = {}".format(eps_label, phi_label))
            return

        h = _build_q2w_hist_for_file(
            root_path,
            "Q2vsW_{}_{}_{}_overlay_fit".format(eps_label.lower(), phi_label, FilenameOverride)
        )
        if not h:
            print("WARNING: Could not build hist for {} epsilon, phi = {}".format(eps_label, phi_label))
            return

        try:
            fits = diamond_fit(h, Q2Val, fitrange=10, threshold=0.0, auto_control=True)
        except Exception as e:
            print("WARNING: diamond_fit failed for {} epsilon, phi = {} ({})".format(eps_label, phi_label, e))
            return

        g = _diamond_graph_from_fits(
            fits,
            "gDiamond_{}_{}_{}".format(eps_label.lower(), phi_label, FilenameOverride)
        )
        if not g:
            print("WARNING: Could not form diamond corners for {} epsilon, phi = {}".format(eps_label, phi_label))
            return

        g.SetLineColor(color)
        g.Draw("L SAME")
        leg.AddEntry(g, "{} #epsilon, {}".format(eps_label, phi_label), "l")

        # keep alive through canvas Print()
        overlay_keepalive.append(g)

    # Define the three phi settings
    phi_defs = [
        ("Center", ["Center", "center"]),
        ("Left",   ["Left", "left"]),
        ("Right",  ["Right", "right"]),
    ]

    # Unique color per (epsilon, phi)
    color_map = {
        ("High", "Center"): ROOT.kBlue,
        ("High", "Left"):   ROOT.kCyan + 1,
        ("High", "Right"):  ROOT.kGreen + 2,
        ("Mid",  "Center"): ROOT.kViolet + 1,
        ("Mid",  "Left"):   ROOT.kAzure + 7,
        ("Mid",  "Right"):  ROOT.kTeal + 3,
        ("Low",  "Center"): ROOT.kRed,
        ("Low",  "Left"):   ROOT.kMagenta + 1,
        ("Low",  "Right"):  ROOT.kOrange + 7,
    }

    def _find_center_fallback_overlay(eps_tag):
        cand = None
        cand_len = 10**9
        for f in glob.glob(OUTPATH + '/*' + ParticleType + '*' + FilenameOverride + '*.root'):
            fl = f.lower()
            if (eps_tag in fl) and ("left" not in fl) and ("right" not in fl):
                if len(f) < cand_len:
                    cand = f
                    cand_len = len(f)
        return cand

    for phi_label, phi_tokens in phi_defs:
        for eps_label, eps_tag in [("High", "high"), ("Mid", "mid"), ("Low", "low")]:
            eps_file = _find_shortest_root(phi_tokens, eps_tag)
            if eps_file is None and phi_label == "Center":
                eps_file = _find_center_fallback_overlay(eps_tag)
            _draw_diamond_outline(eps_file, eps_label, phi_label, color_map[(eps_label, phi_label)])

    # Draw the applied cut polygon as a thick black outline.
    if cut_poly is not None:
        gcut = _tgraph_from_poly(cut_poly, "gCommonCut_{}".format(FilenameOverride))
        if gcut:
            gcut.SetLineColor(ROOT.kBlack)
            gcut.SetLineWidth(5)
            gcut.SetLineStyle(1)
            gcut.Draw("L SAME")
            if common_poly is not None:
                leg.AddEntry(gcut, "Common Cut (overlap)", "l")
            else:
                leg.AddEntry(gcut, "Applied Cut Polygon", "l")
            overlay_keepalive.append(gcut)

    leg.Draw()

    # Close the multipage PDF here
    c1_overlay_allphi.Print(Analysis_Distributions + close_pdf)
    return paramDict
