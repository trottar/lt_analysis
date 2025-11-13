#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-12 14:32:45 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

##################################################################################################################################################

# Import relevant packages
import os
import re
import csv
import math
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ROOT import TFile, TMath, gROOT

import ROOT
from contextlib import contextmanager

@contextmanager
def root_quiet(level=ROOT.kWarning):
    prev = int(ROOT.gErrorIgnoreLevel)
    ROOT.gErrorIgnoreLevel = level
    try:
        yield
    finally:
        ROOT.gErrorIgnoreLevel = prev

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

# -------------------- Configuration --------------------

EPSSET_LST = ["high", "low"]
PHI_SETTING_LST = ["Center", "Left", "Right"]
KIN_VARS = [
    "hsdelta","hsxptar","hsyptar",
    "ssdelta","ssxptar","ssyptar",
    "Q2","W","t","epsilon","MM","th_q","ph_q"
]

# Behavior toggles
SHAPE_ONLY = True         # normalize SIMC to DATA (shape-only comparison)
INCLUDE_OVERFLOW = False
USE_BIN_WIDTH = False
SAVE_PLOTS = True         # diagnostics enabled

# Hard-wired decision thresholds
PVAL_MIN = 1e-10
USE_CHI2_NDF_GATE = True
CHI2_NDF_MIN, CHI2_NDF_MAX = 0.6, 30.0   # acceptably close to 1
USE_HELLINGER_GATE = True
HELLINGER_MAX = 0.15                     # 0 identical → 1 very different
REQUIRE_KS = False        # KS not required for CONTINUE

# -------------------- Helpers --------------------

FNAME_RE = re.compile(
    r"^(?P<pt>[^_]+)_FullAnalysis_Q(?P<Q>[0-9p]+)W(?P<W>[0-9p]+)_(?P<eps>[A-Za-z0-9]+)e\.root$"
)

def parse_tokens(path: Path):
    m = FNAME_RE.match(path.name)
    if not m:
        return None
    pt = m.group("pt")
    Qtok = m.group("Q")
    Wtok = m.group("W")
    eps = m.group("eps")
    return pt, Qtok, Wtok, eps

def list_root_files_for_QW_eps(iter_dir: Path, Q2tok: str, Wtok: str, eps: str,
                               particle: Optional[str] = None) -> List[Path]:
    root_dir = iter_dir / "root"
    if not root_dir.exists():
        return []
    head = f"{particle}_" if particle else "*_"
    patt = f"{head}FullAnalysis_Q{Q2tok}W{Wtok}_{eps}e.root"
    return sorted(root_dir.glob(patt))

def get_hist(file: TFile, path: str):
    obj = file.Get(path)
    if not obj:
        return None
    h = obj.Clone()
    h.SetDirectory(0)
    return h

def dir_exists(file: TFile, path: str) -> bool:
    return bool(file.Get(path))

def same_binning(h1, h2) -> bool:
    if h1.GetNbinsX() != h2.GetNbinsX():
        return False
    x1 = h1.GetXaxis(); x2 = h2.GetXaxis()
    if abs(x1.GetXmin() - x2.GetXmin()) > 1e-9: return False
    if abs(x1.GetXmax() - x2.GetXmax()) > 1e-9: return False
    for i in range(1, h1.GetNbinsX()+1):
        if abs(x1.GetBinLowEdge(i) - x2.GetBinLowEdge(i)) > 1e-9: return False
        if abs(x1.GetBinUpEdge(i)  - x2.GetBinUpEdge(i))  > 1e-9: return False
    return True

def _bin_arrays(h, include_overflow=False):
    nb = h.GetNbinsX()
    idxs = range(0, nb+2) if include_overflow else range(1, nb+1)
    contents, errors, widths = [], [], []
    ax = h.GetXaxis()
    for i in idxs:
        contents.append(float(h.GetBinContent(i)))
        errors.append(float(h.GetBinError(i)))
        widths.append(float(ax.GetBinWidth(i)) if 1 <= i <= nb else 0.0)
    return contents, errors, widths

def _repair_errors(contents, errors):
    """If bin errors are zero/NaN, fall back to Poisson-like sqrt(|content|)."""
    import math
    out = []
    for c, e in zip(contents, errors):
        if e is None or not math.isfinite(e) or e <= 0.0:
            # allow negatives from BG-subtraction
            out.append(math.sqrt(abs(c)) if abs(c) > 0.0 else 0.0)
        else:
            out.append(float(e))
    return out

def _safe_norm(arr):
    import math
    # Clip negatives/NaNs to 0 and renormalize over positives only
    cleaned = [(a if (a is not None and math.isfinite(a) and a > 0.0) else 0.0) for a in arr]
    s = float(sum(cleaned))
    if s <= 0.0:
        return None, 0.0
    return [a/s for a in cleaned], s

def _hellinger(p, q):
    import math
    if p is None or q is None:
        return float("nan")
    # p,q are >=0 and sum to 1 (by _safe_norm). Use math.sqrt to avoid complex.
    s = 0.0
    for pi, qi in zip(p, q):
        # both >=0 from _safe_norm
        s += math.sqrt(pi) * math.sqrt(qi)
    # numeric guard
    s = max(0.0, min(1.0, s))
    return math.sqrt(1.0 - s)

def _js_divergence(p, q, eps=1e-12):
    import math
    if p is None or q is None:
        return float("nan")
    m = [(pi+qi)/2.0 for pi, qi in zip(p, q)]
    def _kl(a, b):
        s = 0.0
        for ai, bi in zip(a, b):
            ai = max(ai, eps); bi = max(bi, eps)
            s += ai * math.log(ai/bi)
        return s
    return 0.5*_kl(p, m) + 0.5*_kl(q, m)

def _wasserstein1_from_bins(p, widths):
    if p is None:
        return lambda q: float("nan")
    pc = []
    acc = 0.0
    for pi in p:
        acc += pi
        pc.append(acc)
    def w1(q):
        if q is None:
            return float("nan")
        acc2 = 0.0
        w = 0.0
        for qi, wi, pci in zip(q, widths, pc):
            acc2 += qi
            w += abs(pci - acc2) * wi
        return w
    return w1

def _poisson_llr(data, mc, shape_only=True):
    n_tot = sum(data); mu_tot = sum(mc)
    alpha = 1.0
    if shape_only and mu_tot > 0:
        alpha = n_tot / mu_tot
    llr = 0.0
    for n, m in zip(data, mc):
        mu = max(alpha*m, 1e-12)
        if n > 0:
            llr += 2.0 * (mu - n + n*math.log(n/mu))
        else:
            llr += 2.0 * mu
    return llr, alpha

def _chi2_with_mc_errors(data, edata, mc, emc, shape_only=True):
    n_tot = sum(data); mu_tot = sum(mc)
    alpha = 1.0
    if shape_only and mu_tot > 0:
        alpha = n_tot / mu_tot
    chi2 = 0.0
    ndof = 0
    for n, en, m, em in zip(data, edata, mc, emc):
        mu = alpha*m
        var = en*en + (alpha*em)**2
        if var <= 0.0:
            continue
        chi2 += (n - mu)**2 / var
        ndof += 1
    if shape_only and ndof > 0:
        ndof -= 1  # fitted normalization
    return chi2, max(ndof, 1), alpha

def compare_th1d(h_data, h_mc,
                 include_overflow=False, use_bin_width=False, shape_only=True) -> Dict[str, float]:
    cD, eD, w = _bin_arrays(h_data, include_overflow=include_overflow)
    cM, eM, _ = _bin_arrays(h_mc,   include_overflow=include_overflow)

    with root_quiet():
        ks_p = h_data.KolmogorovTest(h_mc, "N")
    with root_quiet():
        chi2_p_root = h_data.Chi2Test(h_mc, "UU NORM P")

    # Treat densities as counts if requested
    if use_bin_width:
        cD = [ci*wi for ci, wi in zip(cD, w)]
        cM = [ci*wi for ci, wi in zip(cM, w)]
        eD = [ei*wi for ei, wi in zip(eD, w)]
        eM = [ei*wi for ei, wi in zip(eM, w)]

    # Repair zero/NaN errors so our custom chi2 has usable variances
    eD = _repair_errors(cD, eD)
    eM = _repair_errors(cM, eM)

    # ROOT built-ins: choose options that silence the "NORM needs UU" warning.
    # Use 'UU NORM P' for shape-only; 'WW P' otherwise. Wrap in try/except.
    try:
        ks_p = h_data.KolmogorovTest(h_mc, "N")  # normalized KS (shape)
    except Exception:
        ks_p = float("nan")
    try:
        if shape_only:
            chi2_p_root = h_data.Chi2Test(h_mc, "UU NORM P")
        else:
            chi2_p_root = h_data.Chi2Test(h_mc, "WW P")
    except Exception:
        chi2_p_root = float("nan")

    # Our chi2 with combined uncertainties
    chi2, ndof, alpha_chi2 = _chi2_with_mc_errors(cD, eD, cM, eM, shape_only=shape_only)
    chi2_p_ours = TMath.Prob(chi2, int(ndof)) if ndof > 0 else float("nan")

    # Poisson LLR (always well-defined with small floor in _poisson_llr)
    llr, alpha_llr = _poisson_llr(cD, cM, shape_only=shape_only)

    # Shape distances (robust to negatives via _safe_norm clipping)
    pD, _ = _safe_norm(cD)
    pM, _ = _safe_norm(cM)
    hell = _hellinger(pD, pM)
    jsd  = _js_divergence(pD, pM)
    w1   = _wasserstein1_from_bins(pD, w if not use_bin_width else [1.0]*len(w))(pM) if pD is not None else float("nan")

    return {
        "root_KS_p": ks_p,
        "root_Chi2_p": chi2_p_root,
        "chi2_ndf": (chi2/max(ndof,1)) if ndof > 0 else float("nan"),
        "chi2_p": chi2_p_ours,
        "alpha_norm_used_chi2": alpha_chi2,
        "poisson_llr": llr,
        "alpha_norm_used_llr": alpha_llr,
        "hellinger": hell,
        "js_divergence": jsd,
        "wasserstein_1": w1,
        "shape_only": bool(shape_only),
    }

# -------------------- Plotting (diagnostics) --------------------

def save_overlay_plot(h_data, h_simc, title, outpath, metrics: Dict[str, float], shape_only=True):
    from ROOT import TCanvas, TLegend, kBlue, kRed, gStyle
    gStyle.SetOptStat(0)
    can = TCanvas("c", title, 800, 600)
    hD = h_data
    hM = h_simc.Clone()

    if shape_only:
        aD = hD.Integral(0, hD.GetNbinsX()+1)
        aM = hM.Integral(0, hM.GetNbinsX()+1)
        if aD > 0 and aM > 0:
            hM.Scale(aD/aM)

    hD.SetLineColor(kBlue); hD.SetMarkerColor(kBlue); hD.SetMarkerStyle(20)
    hM.SetLineColor(kRed);  hM.SetMarkerColor(kRed);  hM.SetMarkerStyle(24)

    hmax = max(hD.GetMaximum(), hM.GetMaximum())
    hD.SetMaximum(1.2*hmax if hmax > 0 else 1.0)

    hD.Draw("E1")
    hM.Draw("HIST SAME")

    leg = TLegend(0.12, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.AddEntry(hD, "DATA", "lep")
    leg.AddEntry(hM, "SIMC", "l")
    leg.AddEntry(0, f"KS p={metrics['root_KS_p']:.3g}, Chi2 p={metrics['chi2_p']:.3g}", "")
    leg.AddEntry(0, f"Hellinger={metrics['hellinger']:.3g}, W1={metrics['wasserstein_1']:.3g}", "")
    leg.Draw()

    outpath.parent.mkdir(parents=True, exist_ok=True)
    can.SaveAs(str(outpath))
    can.Close()

# -------------------- Internal directory logic --------------------

# --- replace the old resolvers with these ---

def resolve_hist_paths(file: TFile, phi: Optional[str], kin_var: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Look for histograms at:
      {phi}/data/H_{kin_var}_DATA      and {phi}/simc/H_{kin_var}_SIMC
    and also try a no-phi fallback:
      data/H_{kin_var}_DATA            and simc/H_{kin_var}_SIMC
    """
    candidates = []
    if phi:
        candidates.append((f"{phi}/data/H_{kin_var}_DATA", f"{phi}/simc/H_{kin_var}_SIMC"))
    # no-phi fallback
    candidates.append((f"data/H_{kin_var}_DATA", f"simc/H_{kin_var}_SIMC"))

    for dpath, spath in candidates:
        d = file.Get(dpath)
        s = file.Get(spath)
        if d and s:
            return dpath, spath
    return None, None

def directory_has_phi(file: TFile, phi: str) -> bool:
    """
    Consider a phi 'present' if the directory {phi} exists OR at least one of the
    expected subdirectories {phi}/data or {phi}/simc exists.
    """
    if file.Get(phi):
        return True
    if file.Get(f"{phi}/data") or file.Get(f"{phi}/simc"):
        return True
    # final lightweight probe
    if file.Get(f"{phi}/data/H_Q2_DATA") or file.Get(f"{phi}/simc/H_Q2_SIMC"):
        return True
    return False

# --- replace your existing process_root_file with this version ---

def process_root_file(rpath: Path,
                      eps: str,
                      phi: Optional[str],
                      Q2tok: str,
                      Wtok: str,
                      particle_hint: Optional[str],
                      kin_vars: List[str],
                      shape_only: bool,
                      include_overflow: bool,
                      use_bin_width: bool,
                      save_plots: bool,
                      plots_dir: Path) -> List[Dict[str, object]]:

    rows: List[Dict[str, object]] = []
    f = TFile.Open(str(rpath), "READ")
    if not f or f.IsZombie():
        return rows

    toks = parse_tokens(rpath) or ("", Q2tok, Wtok, eps)
    particle_parsed, Qtok_file, Wtok_file, _ = toks
    particle_use = particle_hint or particle_parsed
    Qtok_use = Qtok_file or Q2tok
    Wtok_use = Wtok_file or Wtok

    # If a phi was requested but not present, record one summary row and return
    if phi and not directory_has_phi(f, phi):
        rows.append({
            "file": str(rpath), "particle": particle_use, "Q2": Qtok_use, "W": Wtok_use,
            "eps": eps, "phi": phi, "kin_var": "", "status": "phi-missing"
        })
        f.Close()
        return rows

    for var in kin_vars:
        dpath, spath = resolve_hist_paths(f, phi, var)
        if not dpath or not spath:
            rows.append({
                "file": str(rpath), "particle": particle_use, "Q2": Qtok_use, "W": Wtok_use,
                "eps": eps, "phi": phi or "", "kin_var": var, "status": "missing-paths"
            })
            continue

        hD = get_hist(f, dpath)
        hS = get_hist(f, spath)
        if not hD or not hS:
            rows.append({
                "file": str(rpath), "particle": particle_use, "Q2": Qtok_use, "W": Wtok_use,
                "eps": eps, "phi": phi or "", "kin_var": var, "status": "missing"
            })
            continue

        if not same_binning(hD, hS):
            rows.append({
                "file": str(rpath), "particle": particle_use, "Q2": Qtok_use, "W": Wtok_use,
                "eps": eps, "phi": phi or "", "kin_var": var, "status": "binning-mismatch"
            })
            continue

        metrics = compare_th1d(
            hD, hS,
            include_overflow=include_overflow,
            use_bin_width=use_bin_width,
            shape_only=shape_only
        )

        print(f"[Chi2] var={var} | PT={particle_use} Q2={Q2tok} W={Wtok} eps={eps} phi={phi or 'NoPhi'} | p={metrics['chi2_p']:.3g} chi2/ndf={metrics['chi2_ndf']:.3g} KS={metrics['root_KS_p']:.3g}")


        row = {
            "file": str(rpath), "particle": particle_use, "Q2": Qtok_use, "W": Wtok_use,
            "eps": eps, "phi": phi or "", "kin_var": var, "status": "ok"
        }
        row.update(metrics)
        rows.append(row)

        if save_plots:
            outname = f"{particle_use or 'PT'}_Q{Q2tok}_W{Wtok}_{eps}_{phi or 'NoPhi'}_{var}.pdf"
            save_overlay_plot(hD, hS, f"{particle_use} Q{Q2tok} W{Wtok} {eps} {phi or ''} :: {var}",
                              plots_dir / outname, metrics, shape_only=shape_only)

    f.Close()
    return rows

# -------------------- Public entry point --------------------

def check_kinematics(inpDict: Dict[str, str], iter_dir: str, iter_num: int) -> Dict[str, object]:
    """
    Inputs:
      inpDict["Q2"]          -> string token like '2p0'
      inpDict["W"]           -> string token like '3p0'
      inpDict["ParticleType"] (optional) -> e.g., 'Proton', 'Kaon', etc.
      iter_dir               -> base iteration directory (contains 'root/')
      iter_num               -> integer for namespacing

    Returns summary including boolean CONTINUE.
    """
    gROOT.SetBatch(True)

    Q2tok: str = str(inpDict["Q2"])
    Wtok:  str = str(inpDict["W"])
    PT:    str = str(inpDict["ParticleType"])

    iter_dir_p = Path(iter_dir)
    monitor_dir = iter_dir_p / "monitor" / "kinematics"
    monitor_dir.mkdir(parents=True, exist_ok=True)

    csv_path = monitor_dir / f"metrics_{PT or 'ALL'}_iter{iter_num}_Q{Q2tok}_W{Wtok}.csv"
    plots_dir = monitor_dir / f"iter{iter_num}" / f"{PT or 'ALL'}_Q{Q2tok}_W{Wtok}"
    if SAVE_PLOTS:
        plots_dir.mkdir(parents=True, exist_ok=True)

    all_rows: List[Dict[str, object]] = []
    eps_missing_files = set()

    for eps in EPSSET_LST:
        rfiles = list_root_files_for_QW_eps(iter_dir_p, Q2tok, Wtok, eps, particle=(PT or None))
        if not rfiles:
            eps_missing_files.add(eps)
            all_rows.append({
                "file": "", "particle": PT, "Q2": Q2tok, "W": Wtok,
                "eps": eps, "phi": "", "kin_var": "", "status": "no-file-for-eps"
            })
            continue

        for rpath in rfiles:
            # Phi-split directories
            for phi in PHI_SETTING_LST:
                rows = process_root_file(
                    rpath=rpath, eps=eps, phi=phi,
                    Q2tok=Q2tok, Wtok=Wtok, particle_hint=PT or None,
                    kin_vars=KIN_VARS,
                    shape_only=SHAPE_ONLY,
                    include_overflow=INCLUDE_OVERFLOW,
                    use_bin_width=USE_BIN_WIDTH,
                    save_plots=SAVE_PLOTS,
                    plots_dir=plots_dir
                )
                all_rows.extend(rows)

            # No-phi subdivision layout
            rows_nophi = process_root_file(
                rpath=rpath, eps=eps, phi=None,
                Q2tok=Q2tok, Wtok=Wtok, particle_hint=PT or None,
                kin_vars=KIN_VARS,
                shape_only=SHAPE_ONLY,
                include_overflow=INCLUDE_OVERFLOW,
                use_bin_width=USE_BIN_WIDTH,
                save_plots=SAVE_PLOTS,
                plots_dir=plots_dir
            )
            if any(r.get("status") == "ok" for r in rows_nophi):
                all_rows.extend(rows_nophi)

    # Write CSV
    fieldnames = [
        "file","particle","Q2","W","eps","phi","kin_var","status",
        "shape_only","root_KS_p","root_Chi2_p","chi2_ndf","chi2_p",
        "alpha_norm_used_chi2","poisson_llr","alpha_norm_used_llr",
        "hellinger","js_divergence","wasserstein_1"
    ]
    with csv_path.open("w", newline="") as fcsv:
        w = csv.DictWriter(fcsv, fieldnames=fieldnames)
        w.writeheader()
        for r in all_rows:
            for k in fieldnames:
                r.setdefault(k, "")
            w.writerow(r)

    # Summary counts
    n_rows = len(all_rows)
    n_ok = sum(1 for r in all_rows if r.get("status") == "ok")
    n_missing = sum(1 for r in all_rows if r.get("status") in ("missing","missing-paths","no-file-for-eps","phi-missing"))
    n_bm = sum(1 for r in all_rows if r.get("status") == "binning-mismatch")

    # CONTINUE decision logic (hard-wired thresholds)
    if eps_missing_files:
        CONTINUE = False
    else:
        def pass_row(r):
            if r.get("status") != "ok":
                return False
            try:
                chi2p   = float(r.get("chi2_p"))
            except Exception:
                chi2p = float("nan")
            try:
                chi2ndf = float(r.get("chi2_ndf"))
            except Exception:
                chi2ndf = float("nan")
            try:
                hell    = float(r.get("hellinger"))
            except Exception:
                hell = float("nan")

            # Hard test
            hard_ok = (chi2p >= PVAL_MIN)

            # Shape-band test
            band_ok = (USE_CHI2_NDF_GATE and not math.isnan(chi2ndf)
                    and CHI2_NDF_MIN <= chi2ndf <= CHI2_NDF_MAX)
            hell_ok = (USE_HELLINGER_GATE and not math.isnan(hell)
                    and hell <= HELLINGER_MAX)

            return hard_ok or (band_ok and hell_ok)


        ok_rows = [r for r in all_rows if r.get("status") == "ok"]
        all_ok_pass = all(pass_row(r) for r in ok_rows) if ok_rows else False
        no_structural_fail = (n_bm == 0) and all(r.get("status") != "missing-paths" for r in all_rows)
        CONTINUE = all_ok_pass and no_structural_fail    

        print(f"\n\nOVERALL: {'PASS' if CONTINUE else 'FAIL'}")   

    return CONTINUE