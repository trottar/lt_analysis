#!/usr/bin/env python3
"""
q2_fit_plot.py (robust)

Fit and visualize Q^2 dependence for structure functions (L, T, LT, TT) using two datasets,
assuming a power-law:  sigma(Q^2) = A * (Q^2)^n

Improvements:
- No requirement for exactly constant Q2; compute an effective Q2 per file:
    * If Q2 has one unique value -> use it.
    * Else if dQ2 exists -> inverse-variance weighted mean of Q2 using dQ2.
    * Else -> median(Q2).
- Recognizes uncertainty columns of the form 'dsigL', 'dsigT', 'dsigLT', 'dsigTT', in addition to *_err/_stat/_sys.
- Flexible SF aliases ('L', 'sigL', 'sigma_L', etc.) while avoiding mapping 'T' to the Mandelstam 't'.
- Merges bins by inferred keys (excluding SFs, their errors, and obvious kinematics). Fallback: merge on 'tbin'.

Outputs:
- One PNG per bin per SF with data points at (Q2_file1, s1) and (Q2_file2, s2), fitted curve,
  a vertical line at Q2_ref, and annotated (A, n, sigma@Q2_ref).
- q2_fit_params.json summarizing A and n for each bin and SF.

Usage:
    python q2_fit_plot.py \
        --file1 KaonLT_Q3p0W2p32.csv \
        --file2 KaonLT_Q4p4W2p74.csv \
        --outdir q2_fits_out \
        --q2_ref 4.3

If automatic bin-key inference fails (e.g., column name mismatches), pass explicit keys:
    --left_on "tbin" --right_on "tbin"
"""

from __future__ import annotations
import argparse, json, math, re
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np, pandas as pd, matplotlib.pyplot as plt

SF_CANONICAL = ["L", "T", "LT", "TT"]
SF_ALIASES = {
    "L":  ["sigL", "sigmaL", "sigma_L", "SL", "sL", "L"],
    "T":  ["sigT", "sigmaT", "sigma_T", "ST", "sT", "T"],  # 'T' last so we don't grab 't'
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
        # Prefer exact case-sensitive single-letter alias
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
    if q2_vals.size == 0:
        return None
    uniq = np.unique(np.round(q2_vals, 6))
    if uniq.size == 1:
        return float(uniq[0])
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
        if c in sf_and_err:
            continue
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

def safe_slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)

def plot_bin_q2_dependence(outdir: Path, sf_name: str, bin_label: str,
                           q2_1: float, q2_2: float, s1: float, s2: float,
                           n: float, A: float, q2_ref: float,
                           yerr1: float | None = None, yerr2: float | None = None):
    if not (np.isfinite(n) and np.isfinite(A)):
        return
    q2_min = min(q2_1, q2_2, q2_ref)
    q2_max = max(q2_1, q2_2, q2_ref)
    xs = np.linspace(0.9*q2_min, 1.1*q2_max, 200)
    xs = xs[xs > 0]
    ys = A * (xs ** n)
    y_ref = A * (q2_ref ** n)
    plt.figure()
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
    outpath = outdir / sf_name / f"{safe_slug(sf_name)}__{safe_slug(bin_label)}.png"
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=200); plt.close()

def main():
    ap = argparse.ArgumentParser(description="Fit Q^2 dependence for L,T,LT,TT from two CSVs and plot.")
    ap.add_argument("--file1", type=Path, required=True)
    ap.add_argument("--file2", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, default="")
    ap.add_argument("--q2_ref", type=float, default=4.3)
    ap.add_argument("--left_on", type=str, default="")
    ap.add_argument("--right_on", type=str, default="")
    args = ap.parse_args()

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

    if args.left_on.strip() and args.right_on.strip():
        left_keys = [s.strip() for s in args.left_on.split(",") if s.strip()]
        right_keys = [s.strip() for s in args.right_on.split(",") if s.strip()]
    else:
        left_keys = pick_binning_keys(df1, sf1)
        right_keys = pick_binning_keys(df2, sf2)

    def cols_to_keep(df, sf_cols_map, keys):
        keep = set(keys)
        for base in sf_cols_map.values():
            keep.add(base)
            for err in related_error_columns(df, base):
                keep.add(err)
        for c in df.columns:
            if c.lower() == "q2":
                keep.add(c)
        return [c for c in df.columns if c in keep]

    df1_s = df1[cols_to_keep(df1, sf1, left_keys)].copy()
    df2_s = df2[cols_to_keep(df2, sf2, right_keys)].copy()

    merged = pd.merge(df1_s, df2_s, left_on=left_keys, right_on=right_keys, how="inner", suffixes=("_f1", "_f2"))
    if merged.empty and ("tbin" in df1.columns and "tbin" in df2.columns):
        merged = pd.merge(df1_s, df2_s, on="tbin", how="inner", suffixes=("_f1", "_f2"))
    if merged.empty:
        raise ValueError("No overlapping bins after merge. Provide explicit --left_on/--right_on with matching keys.")

    params_summary = {sf: [] for sf in common_sfs}

    def bin_label_from_row(row, keys):
        parts = []
        for k in keys:
            if k in row:
                v = row[k]
                parts.append(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}")
        if not parts and "tbin" in row:
            parts.append(f"tbin={row['tbin']}")
        return ", ".join(parts) if parts else "bin"

    for _, row in merged.iterrows():
        bin_label = bin_label_from_row(row, left_keys)
        for sf in common_sfs:
            base1 = sf1[sf]; base2 = sf2[sf]
            s1 = float(row[base1 + "_f1"]); s2 = float(row[base2 + "_f2"])

            def grab_err(df_src, base, suffix):
                dname = "d" + base + suffix
                if dname in merged.columns:
                    return float(row[dname])
                ecols = related_error_columns(df_src, base)
                for nm in ecols:
                    if nm.endswith("_err") and (nm + suffix) in merged.columns:
                        return float(row[nm + suffix])
                stat = sys = None
                if (base + "_stat" + suffix) in merged.columns: stat = float(row[base + "_stat" + suffix])
                if (base + "_sys" + suffix) in merged.columns:  sys  = float(row[base + "_sys" + suffix])
                if stat is not None or sys is not None:
                    stat = 0.0 if stat is None else stat
                    sys  = 0.0 if sys  is None else sys
                    return math.sqrt(stat**2 + sys**2)
                return None

            yerr1 = grab_err(df1, base1, "_f1")
            yerr2 = grab_err(df2, base2, "_f2")

            n, A = compute_n_and_A(s1, s2, q2_1, q2_2)
            if not (np.isfinite(n) and np.isfinite(A)):
                continue

            params_summary[sf].append({
                "bin": bin_label,
                "Q2_file1_eff": q2_1,
                "Q2_file2_eff": q2_2,
                "A": A,
                "n": n,
                "sigma_at_Q2ref": float(A * (args.q2_ref ** n)),
            })

            plot_bin_q2_dependence(args.outdir, sf, bin_label,
                                   q2_1, q2_2, s1, s2, n, A, args.q2_ref,
                                   yerr1=yerr1, yerr2=yerr2)

    args.outdir.mkdir(parents=True, exist_ok=True)
    with open(args.outdir / "q2_fit_params.json", "w") as f:
        json.dump(params_summary, f, indent=2)

    print("Effective Q2:", q2_1, q2_2)
    print("Fitted SFs:", ", ".join(common_sfs))
    print("Plots ->", args.outdir.resolve())
    print("Params ->", (args.outdir / "q2_fit_params.json").resolve())

if __name__ == "__main__":
    main()
