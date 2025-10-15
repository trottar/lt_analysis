#!/usr/bin/env python3
"""
scale_q2_dep.py — per-bin AND per-setting Q^2 dependence with a single multi-page PDF.

Builds on your robust script:
- Flexible SF detection (avoids mapping 'T' to 't')
- Effective Q^2 per file (unique -> value; else weighted mean using dQ2; else median)
- Merge on inferred keys, fallback to 'tbin'

Adds:
- One global Q^2 dependence per setting (L, T, LT, TT)
  * POWERLAW: σ_s(Q^2, tbin) = A_{s,tbin} * (Q^2)^{n_s} with shared n_s
  * DIPOLE:   σ_s(Q^2, tbin) = A_{s,tbin} * (1 + Q^2/Λ_s^2)^(-p_s)
               (fix Λ to fit p, or fix p to fit Λ)
- All plots consolidated into a single multi-page PDF (per-bin pages + one global page per setting)
- JSON of global parameters in the same directory as the PDF

Usage:
    python scale_q2_dep.py \
      --file1 KaonLT_Q3p0W2p32.csv \
      --file2 KaonLT_Q4p4W2p74.csv \
      --outdir q2_out \
      --q2_ref 4.3 \
      --global_model powerlaw
"""

from __future__ import annotations
import argparse, json, math, re
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ------------------- config (as in your base) -------------------
SF_CANONICAL = ["L", "T", "LT", "TT"]
SF_ALIASES = {
    "L":  ["sigL", "sigmaL", "sigma_L", "SL", "sL", "L"],
    "T":  ["sigT", "sigmaT", "sigma_T", "ST", "sT", "T"],  # 'T' last to avoid grabbing 't'
    "LT": ["sigLT", "sigmaLT", "sigma_LT", "S_LT", "sLT", "LT"],
    "TT": ["sigTT", "sigmaTT", "sigma_TT", "S_TT", "sTT", "TT"],
}
ERR_SUFFIXES = ["_err", "_stat", "_sys", "_unc", "_uncorr", "_corr"]
NON_BIN_TOKENS = set(["q2","Q2","w","W","epsilon","eps","Epsilon","dQ2","dW"])

def find_sf_columns(df: pd.DataFrame) -> Dict[str, str]:
    cols = list(df.columns)
    cols_lower = {c.lower(): c for c in cols}
    mapping: Dict[str, str] = {}
    for sf in SF_CANONICAL:
        chosen = None
        # prefer exact one-letter (case-sensitive)
        for alias in SF_ALIASES[sf]:
            if len(alias) == 1 and alias in cols:
                chosen = alias
                break
        if chosen is None:
            for alias in SF_ALIASES[sf]:
                al = alias.lower()
                if al in cols_lower:
                    candidate = cols_lower[al]
                    if sf == "T" and candidate.lower() == "t":
                        continue
                    chosen = candidate
                    break
        if chosen is not None:
            mapping[sf] = chosen
    return mapping

def related_error_columns(df: pd.DataFrame, base_col: str) -> List[str]:
    out = []
    for suf in ERR_SUFFIXES:
        cand = f"{base_col}{suf}"
        if cand in df.columns:
            out.append(cand)
    dbase = "d" + base_col
    if dbase in df.columns:
        out.append(dbase)
    return out

def effective_q2(df: pd.DataFrame) -> float | None:
    q2_col = None
    for c in df.columns:
        if c.lower() == "q2":
            q2_col = c
            break
    if q2_col is None:
        return None
    q2_vals = pd.to_numeric(df[q2_col], errors="coerce").to_numpy()
    q2_vals = q2_vals[np.isfinite(q2_vals)]
    if q2_vals.size == 0: return None
    uniq = np.unique(np.round(q2_vals, 6))
    if uniq.size == 1: return float(uniq[0])
    if "dQ2" in df.columns:
        dq = pd.to_numeric(df["dQ2"], errors="coerce").to_numpy()
        if dq.size == q2_vals.size:
            mask = np.isfinite(dq) & (dq > 0)
            w = np.zeros_like(dq, dtype=float)
            w[mask] = 1.0 / (dq[mask]**2)
            if np.any(w > 0):
                return float(np.average(q2_vals, weights=w))
    return float(np.median(q2_vals))

def pick_binning_keys(df: pd.DataFrame, sf_cols: Dict[str, str]) -> List[str]:
    sf_and_err = set(sf_cols.values())
    for base in sf_cols.values():
        sf_and_err.update(related_error_columns(df, base))
    keys = []
    for c in df.columns:
        if c in sf_and_err: continue
        if c in NON_BIN_TOKENS or c.lower() in [t.lower() for t in NON_BIN_TOKENS]:
            continue
        keys.append(c)
    return keys

def compute_n_and_A(s1: float, s2: float, q2_1: float, q2_2: float) -> Tuple[float, float]:
    if not (np.isfinite(s1) and np.isfinite(s2) and s1 > 0 and s2 > 0):
        return (np.nan, np.nan)
    if not (q2_1 > 0 and q2_2 > 0 and q2_1 != q2_2):
        return (np.nan, np.nan)
    denom = math.log(q2_2 / q2_1)
    n = math.log(s2 / s1) / denom
    A1 = s1 / (q2_1 ** n)
    A2 = s2 / (q2_2 ** n)
    A = math.sqrt(A1 * A2) if A1 > 0 and A2 > 0 else np.nan
    return (n, A)

def ln_err_from_abs(y: float, dy: float | None) -> float | None:
    if dy is None or not np.isfinite(dy) or dy <= 0 or y <= 0:
        return None
    return abs(dy / y)

def safe_slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)

# --------------- plotting: return figures (for single PDF) ---------------
def make_per_bin_fig(sf_name: str, bin_label: str,
                     q2_1: float, q2_2: float, s1: float, s2: float,
                     n: float, A: float, q2_ref: float,
                     yerr1: float | None = None, yerr2: float | None = None):
    if not (np.isfinite(n) and np.isfinite(A)):
        return None
    q2_min = min(q2_1, q2_2, q2_ref)
    q2_max = max(q2_1, q2_2, q2_ref)
    xs = np.linspace(0.9*q2_min, 1.1*q2_max, 200)
    xs = xs[xs > 0]
    ys = A * (xs ** n)
    y_ref = A * (q2_ref ** n)

    fig = plt.figure()
    plt.plot(xs, ys, label="fit: σ(Q²)=A·Q^n")
    if yerr1 is not None and np.isfinite(yerr1) and yerr1 >= 0:
        plt.errorbar([q2_1], [s1], yerr=[yerr1], fmt="o", label=f"data @ Q²={q2_1:.3g}")
    else:
        plt.plot([q2_1], [s1], "o", label=f"data @ Q²={q2_1:.3g}")
    if yerr2 is not None and np.isfinite(yerr2) and yerr2 >= 0:
        plt.errorbar([q2_2], [s2], yerr=[yerr2], fmt="s", label=f"data @ Q²={q2_2:.3g}")
    else:
        plt.plot([q2_2], [s2], "s", label=f"data @ Q²={q2_2:.3g}")
    plt.axvline(q2_ref, linestyle="--", label=f"Q²_ref={q2_ref:g}")
    plt.title(f"{sf_name} — {bin_label}")
    plt.xlabel("Q² (GeV²)"); plt.ylabel(f"σ_{sf_name} (arb.)")
    txt = f"n = {n:.4f}\nA = {A:.4g}\nσ(Q²_ref) = {y_ref:.4g}"
    plt.annotate(txt, xy=(0.02, 0.98), xycoords="axes fraction", va="top", ha="left")
    plt.legend(); plt.tight_layout()
    return fig

def make_global_powerlaw_fig(setting_name: str, q2_1: float, q2_2: float, q2_ref: float,
                             y_points: List[Tuple[str, Tuple[float,float,float|None], Tuple[float,float,float|None]]],
                             n_s: float, dn_s: float):
    xs = np.linspace(0.9*min(q2_1,q2_2,q2_ref), 1.1*max(q2_1,q2_2,q2_ref), 200)
    fig = plt.figure()
    # Overlay each bin's line with shared slope n_s and its own intercept through (Q1,s1)
    for (label,(Q1,s1,dy1),(Q2,s2,dy2)) in y_points:
        if not (s1>0 and s2>0): continue
        Aj = s1 / (Q1**n_s)
        plt.plot(xs, Aj*(xs**n_s), alpha=0.55)
        if dy1 is not None: plt.errorbar([Q1],[s1],yerr=[dy1],fmt="o", label=f"{label}@Q²={Q1:g}")
        else:               plt.plot([Q1],[s1],"o", label=f"{label}@Q²={Q1:g}")
        if dy2 is not None: plt.errorbar([Q2],[s2],yerr=[dy2],fmt="s")
        else:               plt.plot([Q2],[s2],"s")
    plt.axvline(q2_ref, linestyle="--", label=f"Q²_ref={q2_ref:g}")
    plt.title(f"{setting_name} — global power law: n={n_s:.3f}±{dn_s:.3f}")
    plt.xlabel("Q² (GeV²)"); plt.ylabel(f"σ_{setting_name} (arb.)")
    plt.legend(ncol=2, fontsize=8); plt.tight_layout()
    return fig

# -------------------- global per-setting fits --------------------
def combine_bin_slopes(bin_slopes: List[Tuple[float,float]]) -> Tuple[float,float]:
    """bin_slopes: (n_j, dn_j). Returns (n, dn) with inverse-variance weighting."""
    ws, xs = [], []
    for n_j, dn_j in bin_slopes:
        if np.isfinite(n_j) and np.isfinite(dn_j) and dn_j > 0:
            ws.append(1.0 / (dn_j**2))
            xs.append(n_j)
    if not ws:
        return (np.nan, np.nan)
    wsum = float(np.sum(ws))
    n = float(np.sum(np.array(ws) * np.array(xs)) / wsum)
    dn = math.sqrt(1.0 / wsum)
    return (n, dn)

def fit_global_powerlaw_for_setting(merged: pd.DataFrame, sf_col1: str, sf_col2: str,
                                    q2_1: float, q2_2: float, q2_ref: float,
                                    setting_name: str,
                                    collect_figs: List):
    # Build per-bin slope list and plotting payload
    bin_slopes = []
    y_points = []
    for idx, row in merged.iterrows():
        s1 = float(row[sf_col1 + "_f1"]); s2 = float(row[sf_col2 + "_f2"])
        if not (s1>0 and s2>0): 
            continue
        n_j = (math.log(s2) - math.log(s1)) / (math.log(q2_2) - math.log(q2_1))
        dy1 = row.get("d"+sf_col1 + "_f1", None)
        dy2 = row.get("d"+sf_col2 + "_f2", None)
        e1 = ln_err_from_abs(s1, float(dy1) if dy1 is not None else None)
        e2 = ln_err_from_abs(s2, float(dy2) if dy2 is not None else None)
        var_dln = 0.0
        if e1: var_dln += e1**2
        if e2: var_dln += e2**2
        dn_j = math.sqrt(var_dln) / abs(math.log(q2_2) - math.log(q2_1)) if var_dln>0 else np.nan
        bin_slopes.append((n_j, dn_j))
        label = f"tbin={row['tbin']}" if "tbin" in row else f"row{idx}"
        y_points.append((label, (q2_1, s1, float(dy1) if dy1 is not None else None),
                                (q2_2, s2, float(dy2) if dy2 is not None else None)))
    n_s, dn_s = combine_bin_slopes(bin_slopes)
    fig = make_global_powerlaw_fig(setting_name, q2_1, q2_2, q2_ref, y_points, n_s, dn_s)
    collect_figs.append(fig)
    return n_s, dn_s

def fit_global_dipole_one_free(merged: pd.DataFrame, sf_col1: str, sf_col2: str,
                               q2_1: float, q2_2: float, q2_ref: float,
                               setting_name: str,
                               fix_Lambda: float | None, fix_p: float | None,
                               collect_figs: List):
    assert (fix_Lambda is None) ^ (fix_p is None), "Fix exactly one of (Lambda, p)."
    if fix_Lambda is not None:
        F1 = 1.0 + q2_1/(fix_Lambda**2)
        F2 = 1.0 + q2_2/(fix_Lambda**2)
        dlnF = math.log(F2) - math.log(F1)
        p_vals, wp = [], []
        y_points = []
        for idx, row in merged.iterrows():
            s1 = float(row[sf_col1 + "_f1"]); s2 = float(row[sf_col2 + "_f2"])
            if s1 <= 0 or s2 <= 0: continue
            pj = - (math.log(s2) - math.log(s1)) / dlnF
            dy1 = row.get("d"+sf_col1 + "_f1", None)
            dy2 = row.get("d"+sf_col2 + "_f2", None)
            e1 = ln_err_from_abs(s1, float(dy1) if dy1 is not None else None)
            e2 = ln_err_from_abs(s2, float(dy2) if dy2 is not None else None)
            var = (e1 or 0)**2 + (e2 or 0)**2
            if var > 0:
                w = 1.0 / (var / (dlnF**2))
                p_vals.append(pj); wp.append(w)
            label = f"tbin={row['tbin']}" if "tbin" in row else f"row{idx}"
            y_points.append((label, (q2_1, s1, float(dy1) if dy1 is not None else None),
                                    (q2_2, s2, float(dy2) if dy2 is not None else None)))
        p = float(np.average(p_vals, weights=wp)) if wp else np.nan
        dp = math.sqrt(1.0/np.sum(wp)) if wp else np.nan
        # Plot overlay with shared p and fixed Λ
        xs = np.linspace(0.9*min(q2_1,q2_2,q2_ref), 1.1*max(q2_1,q2_2,q2_ref), 200)
        fig = plt.figure()
        for (label,(Q1,s1,dy1),(Q2,s2,dy2)) in y_points:
            Aj = s1 / ((1.0 + Q1/(fix_Lambda**2))**(-p))
            ys = Aj * ((1.0 + xs/(fix_Lambda**2))**(-p))
            plt.plot(xs, ys, alpha=0.55)
            if dy1 is not None: plt.errorbar([Q1],[s1],yerr=[dy1],fmt="o", label=f"{label}@Q²={Q1:g}")
            else:               plt.plot([Q1],[s1],"o", label=f"{label}@Q²={Q1:g}")
            if dy2 is not None: plt.errorbar([Q2],[s2],yerr=[dy2],fmt="s")
            else:               plt.plot([Q2],[s2],"s")
        plt.axvline(q2_ref, linestyle="--", label=f"Q²_ref={q2_ref:g}")
        plt.title(f"{setting_name} — dipole (Λ={fix_Lambda:.3g} GeV): p={p:.3f}±{dp:.3f}")
        plt.xlabel("Q² (GeV²)"); plt.ylabel(f"σ_{setting_name} (arb.)")
        plt.legend(ncol=2, fontsize=8); plt.tight_layout()
        collect_figs.append(fig)
        return (p, dp, fix_Lambda, None)

    # fix_p is not None, fit Λ by 1D scan (weighted)
    def pred_delta_ln(L):
        return (-fix_p) * (math.log(1+q2_2/L**2) - math.log(1+q2_1/L**2))
    def objective(L):
        if L <= 0: return np.inf
        pred = pred_delta_ln(L)
        num = den = 0.0
        for _, row in merged.iterrows():
            s1 = float(row[sf_col1 + "_f1"]); s2 = float(row[sf_col2 + "_f2"])
            if s1 <= 0 or s2 <= 0: continue
            y = (math.log(s2) - math.log(s1))
            dy1 = row.get("d"+sf_col1 + "_f1", None)
            dy2 = row.get("d"+sf_col2 + "_f2", None)
            e1 = ln_err_from_abs(s1, float(dy1) if dy1 is not None else None)
            e2 = ln_err_from_abs(s2, float(dy2) if dy2 is not None else None)
            var = (e1 or 0)**2 + (e2 or 0)**2
            w = 0.0 if var==0 else 1.0/var
            num += w * (y - pred)**2
            den += w
        return num/den if den>0 else np.inf
    L_grid = np.geomspace(0.2, 3.0, 200)  # GeV (tune if you like)
    vals = [objective(L) for L in L_grid]
    i = int(np.argmin(vals))
    Lhat = float(L_grid[i])
    # crude error from curvature in log L
    if 1 < i < len(L_grid)-2:
        h = math.log(L_grid[i+1]) - math.log(L_grid[i])
        curv = (vals[i+1] - 2*vals[i] + vals[i-1]) / (h*h)
        dL = None if curv<=0 else Lhat * math.sqrt(2.0/curv)
    else:
        dL = None

    # Plot overlay with fixed p and fitted Λ
    xs = np.linspace(0.9*min(q2_1,q2_2,q2_ref), 1.1*max(q2_1,q2_2,q2_ref), 200)
    fig = plt.figure()
    for idx, row in merged.iterrows():
        s1 = float(row[sf_col1 + "_f1"]); s2 = float(row[sf_col2 + "_f2"])
        if s1 <= 0 or s2 <= 0: continue
        Q1, Q2 = q2_1, q2_2
        Aj = s1 / ((1.0 + Q1/(Lhat**2))**(-fix_p))
        ys = Aj * ((1.0 + xs/(Lhat**2))**(-fix_p))
        label = f"tbin={row['tbin']}" if "tbin" in row else f"row{idx}"
        dy1 = row.get("d"+sf_col1 + "_f1", None)
        dy2 = row.get("d"+sf_col2 + "_f2", None)
        plt.plot(xs, ys, alpha=0.55)
        if dy1 is not None: plt.errorbar([Q1],[s1],yerr=[dy1],fmt="o", label=f"{label}@Q²={Q1:g}")
        else:               plt.plot([Q1],[s1],"o", label=f"{label}@Q²={Q1:g}")
        if dy2 is not None: plt.errorbar([Q2],[s2],yerr=[dy2],fmt="s")
        else:               plt.plot([Q2],[s2],"s")
    plt.axvline(q2_ref, linestyle="--", label=f"Q²_ref={q2_ref:g}")
    ttl = f"{setting_name} — dipole (p={fix_p:.3g} fixed): Λ={Lhat:.3g}"
    ttl += f" ± {dL:.3g}" if dL is not None else ""
    plt.title(ttl)
    plt.xlabel("Q² (GeV²)"); plt.ylabel(f"σ_{setting_name} (arb.)")
    plt.legend(ncol=2, fontsize=8); plt.tight_layout()
    collect_figs.append(fig)
    return (fix_p, None, Lhat, dL)

# -------------------- main --------------------
def main():
    ap = argparse.ArgumentParser(description="Q^2 dependence: per-bin and per-setting (single PDF).")
    ap.add_argument("--file1", type=Path, required=True)
    ap.add_argument("--file2", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--q2_ref", type=float, default=4.3)
    ap.add_argument("--left_on", type=str, default="")
    ap.add_argument("--right_on", type=str, default="")
    ap.add_argument("--global_model", choices=["powerlaw","dipole"], default="powerlaw")
    ap.add_argument("--dipole_fix_Lambda", type=float, default=None, help="GeV; if set, fit p")
    ap.add_argument("--dipole_fix_p", type=float, default=None, help="dimensionless; if set, fit Λ")
    ap.add_argument("--no_per_bin_plots", action="store_true")
    args = ap.parse_args()

    # Load
    df1 = pd.read_csv(args.file1); df2 = pd.read_csv(args.file2)
    q2_1 = effective_q2(df1); q2_2 = effective_q2(df2)
    if q2_1 is None or q2_2 is None:
        raise ValueError("Could not determine effective Q2 from one or both files (need a 'Q2' column).")

    sf1 = find_sf_columns(df1); sf2 = find_sf_columns(df2)
    common_sfs = sorted(set(sf1.keys()) & set(sf2.keys()))
    if "T" in common_sfs and (sf1["T"].lower() == "t" or sf2["T"].lower() == "t"):
        common_sfs.remove("T")
    if not common_sfs:
        raise ValueError(f"No common structure functions detected. file1: {sf1}  file2: {sf2}")

    # Merge
    if args.left_on.strip() and args.right_on.strip():
        left_keys = [s.strip() for s in args.left_on.split(",") if s.strip()]
        right_keys = [s.strip() for s in args.right_on.split(",") if s.strip()]
    else:
        left_keys  = pick_binning_keys(df1, sf1)
        right_keys = pick_binning_keys(df2, sf2)

    def cols_to_keep(df, sf_cols_map, keys):
        keep = set(keys)
        for base in sf_cols_map.values():
            keep.add(base)
            for err in related_error_columns(df, base):
                keep.add(err)
        for c in df.columns:
            if c.lower() == "q2": keep.add(c)
        if "tbin" in df.columns: keep.add("tbin")
        return [c for c in df.columns if c in keep]

    df1_s = df1[cols_to_keep(df1, sf1, left_keys)].copy()
    df2_s = df2[cols_to_keep(df2, sf2, right_keys)].copy()

    merged = pd.merge(df1_s, df2_s, left_on=left_keys, right_on=right_keys, how="inner", suffixes=("_f1", "_f2"))
    if merged.empty and ("tbin" in df1.columns and "tbin" in df2.columns):
        merged = pd.merge(df1_s, df2_s, on="tbin", how="inner", suffixes=("_f1", "_f2"))
    if merged.empty:
        raise ValueError("No overlapping bins after merge. Provide explicit --left_on/--right_on.")

    # Prepare outputs
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    pdf_path = outdir / "plots.pdf"
    json_path = outdir / "global_q2_params.json"

    perbin_params = {sf: [] for sf in common_sfs}
    global_params = {}

    # Helper
    def bin_label_from_row(row, keys):
        parts = []
        for k in keys:
            if k in row: parts.append(f"{k}={row[k]}")
        if not parts and "tbin" in row: parts.append(f"tbin={row['tbin']}")
        return ", ".join(parts) if parts else "bin"

    figs: List = []
    # Per-bin pages
    if not args.no_per_bin_plots:
        for _, row in merged.iterrows():
            bin_label = bin_label_from_row(row, left_keys)
            for sf in common_sfs:
                base1 = sf1[sf]; base2 = sf2[sf]
                s1 = float(row[base1 + "_f1"]); s2 = float(row[base2 + "_f2"])
                n, A = compute_n_and_A(s1, s2, q2_1, q2_2)
                if not (np.isfinite(n) and np.isfinite(A)): continue
                # pick errors if present
                dy1 = row.get("d"+base1 + "_f1", None)
                dy2 = row.get("d"+base2 + "_f2", None)
                yerr1 = float(dy1) if dy1 is not None and np.isfinite(dy1) else None
                yerr2 = float(dy2) if dy2 is not None and np.isfinite(dy2) else None
                perbin_params[sf].append({"bin": bin_label, "A": A, "n": n})
                fig = make_per_bin_fig(sf, bin_label, q2_1, q2_2, s1, s2, n, A, args.q2_ref, yerr1, yerr2)
                if fig is not None: figs.append(fig)

    # Global pages
    for sf in common_sfs:
        base1 = sf1[sf]; base2 = sf2[sf]
        if args.global_model == "powerlaw":
            n_s, dn_s = fit_global_powerlaw_for_setting(
                merged, base1, base2, q2_1, q2_2, args.q2_ref, sf, figs
            )
            global_params[sf] = {"model": "powerlaw", "n": n_s, "dn": dn_s}
        else:
            if (args.dipole_fix_Lambda is None) == (args.dipole_fix_p is None):
                raise ValueError("For dipole model, set exactly one of --dipole_fix_Lambda or --dipole_fix_p.")
            p, dp, Lam, dLam = fit_global_dipole_one_free(
                merged, base1, base2, q2_1, q2_2, args.q2_ref, sf,
                args.dipole_fix_Lambda, args.dipole_fix_p, figs
            )
            global_params[sf] = {"model": "dipole", "p": p, "dp": dp, "Lambda": Lam, "dLambda": dLam}

    # Write JSON and single PDF
    with open(json_path, "w") as f:
        json.dump(global_params, f, indent=2)
    with PdfPages(pdf_path) as pdf:
        for fig in figs:
            pdf.savefig(fig)
            plt.close(fig)

    print("Saved JSON:", json_path.resolve())
    print("Saved PDF:", pdf_path.resolve())

if __name__ == "__main__":
    main()
