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
ROOT.gErrorIgnoreLevel = ROOT.kWarning  # hides all "Info" prints (only shows Warning+)

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

# thresholds
F_SYS_YIELD = 0.05            # 5% floor for yields
USE_PVAL = False              # p-value not used for gating
USE_KS   = False              # KS off (binned + ties -> ~0)
CHI2_NDF_MIN, CHI2_NDF_MAX = 0.0, 50.0
HELLINGER_MAX = 0.20
W1NORM_MAX = 0.02             # 2% of x-range
PULL95_MAX = 3.0              # 95th percentile |pull| ≤ 3σ
MEAN_REL_MAX = 0.05           # |Δmean| ≤ 5% of x-range
RMS_REL_MAX  = 0.10           # |ΔRMS|  ≤ 10% of x-range
PVAL_MIN = 1e-10   # only used if USE_PVAL is True

# -------------------- Helpers --------------------

def _x_range(h):
    ax = h.GetXaxis()
    return float(ax.GetXmax() - ax.GetXmin())

def _quantiles(vals, ps):
    if not vals:
        return [float("nan") for _ in ps]
    s = sorted(vals)
    n = len(s)
    out = []
    for p in ps:
        i = min(max(int(round(p*(n-1))), 0), n-1)
        out.append(float(s[i]))
    return out

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

    # Treat densities as counts if requested
    if use_bin_width:
        cD = [ci*wi for ci, wi in zip(cD, w)]
        cM = [ci*wi for ci, wi in zip(cM, w)]
        eD = [ei*wi for ei, wi in zip(eD, w)]
        eM = [ei*wi for ei, wi in zip(eM, w)]

    # Repair zero/NaN errors so our custom chi2 has usable variances
    eD = _repair_errors(cD, eD)
    eM = _repair_errors(cM, eM)

    eD = [math.hypot(ed, F_SYS_YIELD * abs(cd)) for ed, cd in zip(eD, cD)]
    eM = [math.hypot(em, F_SYS_YIELD * abs(cm)) for em, cm in zip(eM, cM)]    

    # ROOT built-ins: choose options that silence the "NORM needs UU" warning.
    # Use 'UU NORM P' for shape-only; 'WW P' otherwise. Wrap in try/except.
    try:
        ks_p = h_data.KolmogorovTest(h_mc, "N")  # normalized KS (shape)
    except Exception:
        ks_p = float("nan")
    try:
        if shape_only:
            chi2_p_root = h_data.Chi2Test(h_mc, "UU NORM")
        else:
            chi2_p_root = h_data.Chi2Test(h_mc, "WW")
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

    # Per-bin pulls (using the same variance model)
    var = [ed*ed + (alpha_chi2*em)**2 for ed, em in zip(eD, eM)]
    pulls = []
    for d, m, v in zip(cD, cM, var):
        if v > 0:
            pulls.append((d - alpha_chi2*m) / math.sqrt(v))
    abs_pulls = [abs(x) for x in pulls]
    pull_p50, pull_p95 = _quantiles(abs_pulls, [0.50, 0.95])

    # Shape metrics (you already compute Hellinger; add W1 normalized)
    xr = max(_x_range(h_data), 1e-12)
    w1_norm = w1 / xr if (isinstance(w1, float) and math.isfinite(w1)) else float("nan")

    # Mean / RMS deltas (scale-free, relative to x-range)
    meanD, meanM = h_data.GetMean(), h_mc.GetMean()
    rmsD,  rmsM  = h_data.GetRMS(),  h_mc.GetRMS()
    mean_rel = abs(meanD - meanM) / xr
    rms_rel  = abs(rmsD  - rmsM ) / xr

    return {
        "root_KS_p": ks_p,
        "root_Chi2_p": chi2_p_root,
        "chi2_ndf": (chi2/max(ndof,1)) if ndof > 0 else float("nan"),
        "chi2_p": chi2_p_ours,
        "alpha_norm_used_chi2": alpha_chi2,
        "poisson_llr": llr,
        "alpha_norm_used_llr": alpha_llr,
        "hellinger": hell,
        "wasserstein_1": w1,
        "w1_norm": w1_norm,
        "pull_p50": pull_p50,
        "pull_p95": pull_p95,
        "mean_rel": mean_rel,
        "rms_rel": rms_rel,
        "shape_only": bool(shape_only),
    }

# -------------------- Plotting (diagnostics) --------------------

# REPLACE your entire save_overlay_plot() with this version (no KS, no p; shows chi2/ndf + shape metrics)

def save_overlay_plot(h_data, h_simc, title, outpath, metrics, shape_only=True):
    from ROOT import TCanvas, TLegend, kBlue, kRed, gStyle
    gStyle.SetOptStat(0)

    can = TCanvas("c", title, 800, 600)
    hD = h_data
    hM = h_simc.Clone()

    # shape-only normalization (match areas)
    if shape_only:
        aD = hD.Integral(0, hD.GetNbinsX()+1)
        aM = hM.Integral(0, hM.GetNbinsX()+1)
        if aD > 0 and aM > 0:
            hM.Scale(aD/aM)

    # styles (smaller markers/legend)
    hD.SetLineColor(kBlue); hD.SetMarkerColor(kBlue); hD.SetMarkerStyle(20); hD.SetMarkerSize(0.6); hD.SetLineWidth(1)
    hM.SetLineColor(kRed);  hM.SetMarkerColor(kRed);  hM.SetMarkerStyle(24); hM.SetMarkerSize(0.6); hM.SetLineWidth(1)

    hmax = max(hD.GetMaximum(), hM.GetMaximum())
    hD.SetMaximum(1.2*hmax if hmax > 0 else 1.0)

    hD.Draw("E1")
    hM.Draw("HIST SAME")

    # derive W1% on the fly in case w1_norm isn't in metrics
    try:
        xr = float(hD.GetXaxis().GetXmax() - hD.GetXaxis().GetXmin())
        w1 = float(metrics.get("wasserstein_1", float("nan")))
        w1pct = (w1/xr*100.0) if (xr > 0 and math.isfinite(w1)) else float("nan")
    except Exception:
        w1pct = float("nan")

    # legend (compact, no KS/p-values)
    leg = TLegend(0.60, 0.75, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.028)
    leg.SetMargin(0.12)

    leg.AddEntry(hD, "DATA", "lep")
    leg.AddEntry(hM, "SIMC", "l")
    leg.AddEntry(0, f"χ²/ndf={metrics.get('chi2_ndf', float('nan')):.2f}", "")
    leg.AddEntry(0, f"Hellinger={metrics.get('hellinger', float('nan')):.3f}, W1%={w1pct:.2f}", "")

    # show optional diagnostics if present
    if 'pull_p95' in metrics:
        leg.AddEntry(0, f"P95|pull|={metrics['pull_p95']:.2f}", "")
    if 'mean_rel' in metrics and 'rms_rel' in metrics:
        leg.AddEntry(0, f"dμ%={metrics['mean_rel']*100:.2f}, dRMS%={metrics['rms_rel']*100:.2f}", "")

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

        print(
            f"{var} | Q2={Q2tok} W={Wtok} eps={eps} phi={phi or 'NoPhi'} | "
            f"chi2/ndf={metrics['chi2_ndf']:.3g}  H={metrics['hellinger']:.3g}  "
            f"W1%={(metrics.get('w1_norm', metrics['wasserstein_1']/(hD.GetXaxis().GetXmax()-hD.GetXaxis().GetXmin()))*100):.2f}"
        )

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
    plots_dir = monitor_dir / "plots"
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
    base = [
        "file","particle","Q2","W","eps","phi","kin_var","status",
        "shape_only","chi2","ndof","chi2_ndf","chi2_p","root_Chi2_p","root_KS_p",
        "hellinger","js_divergence","wasserstein_1","poisson_llr",
        "alpha_norm_used_chi2","alpha_norm_used_llr",
    ]
    # union of keys across rows, preserving order (base first)
    extra_keys = []
    seen = set(base)
    for r in all_rows:
        for k in r.keys():
            if k not in seen:
                seen.add(k)
                extra_keys.append(k)
    fieldnames = base + extra_keys

    with csv_path.open("w", newline="") as fcsv:
        w = csv.DictWriter(fcsv, fieldnames=fieldnames)
        w.writeheader()
        for r in all_rows:
            for k in fieldnames:
                r.setdefault(k, "")
            r = {k: r.get(k, "") for k in fieldnames}
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
                chi2ndf = float(r["chi2_ndf"])
                hell    = float(r["hellinger"])
                w1n     = float(r["w1_norm"])
                p95     = float(r["pull_p95"])
                mrel    = float(r["mean_rel"])
                rrel    = float(r["rms_rel"])
            except Exception:
                return False

            band_ok   = CHI2_NDF_MIN <= chi2ndf <= CHI2_NDF_MAX
            shape_ok  = (hell <= HELLINGER_MAX) and (w1n <= W1NORM_MAX)
            pulls_ok  = (p95 <= PULL95_MAX)
            moments_ok= (mrel <= MEAN_REL_MAX) and (rrel <= RMS_REL_MAX)

            return band_ok and shape_ok and pulls_ok and moments_ok

        # Meaning of checks:
        #  - EPSILON SETS   : required epsilon settings have files (e.g., 'high', 'low')
        #  - FILE STRUCTURE : no structural issues (no binning mismatches, no missing histogram paths)
        #  - METRIC ROWS    : every 'ok' comparison row passes your thresholds

        eps_ok = (len(eps_missing_files) == 0)

        missing_paths = sum(1 for r in all_rows if r.get("status") == "missing-paths")
        bin_mismatch  = n_bm
        struct_ok = (bin_mismatch == 0) and (missing_paths == 0)

        ok_rows  = [r for r in all_rows if r.get("status") == "ok"]
        ok_pass  = sum(1 for r in ok_rows if pass_row(r))
        rows_pass = bool(ok_rows) and (ok_pass == len(ok_rows))

        CONTINUE = eps_ok and struct_ok and rows_pass

        print(f"\n\nEPSILON SETS : {'PASS' if eps_ok else 'FAIL'} | Missing={','.join(sorted(eps_missing_files)) if eps_missing_files else 'none'}")
        print(f"FILE STRUCTURE: {'PASS' if struct_ok else 'FAIL'} | Bin Mismatch={bin_mismatch}, Missing Paths={missing_paths}")
        print(f"METRIC ROWS  : {'PASS' if rows_pass else 'FAIL'} | Pass={ok_pass}/{len(ok_rows)}")

        # ---- sub-breakdown of metrics (quick-glance) ----
        def _finite_float(x):
            try:
                v = float(x)
                return (not math.isnan(v)), v
            except Exception:
                return False, float('nan')

        cndf_pass = hell_pass = w1_pass = p95_pass = mean_pass = rms_pass = 0
        fail_summaries = []
        fail_counts = {"chi2/ndf": 0, "H": 0, "W1": 0, "pull95": 0, "dμ": 0, "dRMS": 0}

        for r in ok_rows:
            reasons = []

            ok_cndf, chi2ndf = _finite_float(r.get("chi2_ndf"))
            if ok_cndf and (CHI2_NDF_MIN <= chi2ndf <= CHI2_NDF_MAX):
                cndf_pass += 1
            else:
                reasons.append("chi2/ndf")

            ok_hell, hell = _finite_float(r.get("hellinger"))
            if ok_hell and (hell <= HELLINGER_MAX):
                hell_pass += 1
            else:
                reasons.append("H")

            ok_w1, w1n = _finite_float(r.get("w1_norm"))
            if ok_w1 and (w1n <= W1NORM_MAX):
                w1_pass += 1
            else:
                reasons.append("W1")

            ok_p95, p95 = _finite_float(r.get("pull_p95"))
            if ok_p95 and (p95 <= PULL95_MAX):
                p95_pass += 1
            else:
                reasons.append("pull95")

            ok_m, mrel = _finite_float(r.get("mean_rel"))
            if ok_m and (mrel <= MEAN_REL_MAX):
                mean_pass += 1
            else:
                reasons.append("dμ")

            ok_r, rrel = _finite_float(r.get("rms_rel"))
            if ok_r and (rrel <= RMS_REL_MAX):
                rms_pass += 1
            else:
                reasons.append("dRMS")

            if reasons:
                for k in reasons:
                    fail_counts[k] += 1
                fail_summaries.append(
                    f"{r.get('kin_var','?')}[{r.get('eps','?')}/{r.get('phi','?')}]:" + ",".join(reasons)
                )

        print("=METRIC BREAKDOWN:")
        print(f"  chi2/ndf in [{CHI2_NDF_MIN},{CHI2_NDF_MAX}]: {cndf_pass}/{len(ok_rows)}")
        print(f"  Hellinger ≤ {HELLINGER_MAX}: {hell_pass}/{len(ok_rows)}")
        print(f"  W1% ≤ {W1NORM_MAX*100:.2f}: {w1_pass}/{len(ok_rows)}")
        print(f"  P95|pull| ≤ {PULL95_MAX:.2f}: {p95_pass}/{len(ok_rows)}")
        print(f"  |Δμ|% ≤ {MEAN_REL_MAX*100:.2f}: {mean_pass}/{len(ok_rows)}")
        print(f"  |ΔRMS|% ≤ {RMS_REL_MAX*100:.2f}: {rms_pass}/{len(ok_rows)}")

        # show the AND-gate result explicitly (this matches METRIC ROWS Pass=X/Y)
        print(f"  COMBINED (chi2+H+W1+pull95+dμ+dRMS): {ok_pass}/{len(ok_rows)}")

        # per-reason fail counts (helps see why the AND drops)
        if any(fail_counts.values()):
            print("  fail counts:", " ".join(f"{k}={v}" for k, v in fail_counts.items()))

        # chunked fail list, 5 per line
        if fail_summaries:
            n_per_line = 5
            for i in range(0, len(fail_summaries), n_per_line):
                line = "; ".join(fail_summaries[i:i+n_per_line])
                if i == 0:
                    print(f"  fails: {line}")
                else:
                    print(f"         {line}")
        else:
            print("  fails: none")

        print(f"OVERALL: {'PASS' if CONTINUE else 'FAIL'}")

    return CONTINUE