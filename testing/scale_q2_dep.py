#!/usr/bin/env python3
"""
q2_fit_plot.py

Fit and visualize Q^2 dependence for structure functions (L, T, LT, TT) using two datasets,
assuming a power-law form:  sigma(Q^2) = A * (Q^2)^n

- Each input CSV must contain a constant 'Q2' column (one unique value).
- Flexible column names are supported for SFs (e.g., 'L', 'sigma_L', 'sigLT', etc.).
- Uncertainties (if present) are carried to plots as y-error bars (suffixes like _err, _stat, _sys).
- Per-bin fits use the two Q^2 points to extract n and A:
      n = ln(s2/s1) / ln(Q2_2/Q2_1)
      A = geometric mean of A1 and A2, where A1 = s1/Q2_1^n and A2 = s2/Q2_2^n
- Produces one PNG per bin per SF with the data points, the fitted curve, and the value at Q2_ref.
- Writes a parameters JSON for all bins and all SFs.

Usage:
    python q2_fit_plot.py \
        --file1 KaonLT_Q3p0W2p32.csv \
        --file2 KaonLT_Q4p4W2p74.csv \
        --outdir q2_fits_out \
        --q2_ref 4.3

If automatic bin-key inference fails (e.g., column name mismatches), pass explicit keys:
    --left_on "t,phi" --right_on "t,phi"
"""

from __future__ import annotations
import argparse
import json
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ----- Config -----
SF_CANONICAL = ["L", "T", "LT", "TT"]
SF_ALIASES = {
    "L":  ["L", "sigL", "sigmaL", "sigma_L", "SL", "sL"],
    "T":  ["T", "sigT", "sigmaT", "sigma_T", "ST", "sT"],
    "LT": ["LT", "sigLT", "sigmaLT", "sigma_LT", "S_LT", "sLT"],
    "TT": ["TT", "sigTT", "sigmaTT", "sigma_TT", "S_TT", "sTT"],
}
ERR_SUFFIXES = ["_err", "_stat", "_sys", "_unc", "_uncorr", "_corr"]
NON_BIN_TOKENS = set(["q2", "Q2", "w", "W", "epsilon", "eps", "Epsilon"])


# ----- Utilities -----
def find_sf_columns(df: pd.DataFrame) -> Dict[str, str]:
    cols_lower = {c.lower(): c for c in df.columns}
    mapping: Dict[str, str] = {}
    for sf in SF_CANONICAL:
        for alias in SF_ALIASES[sf]:
            if alias.lower() in cols_lower:
                mapping[sf] = cols_lower[alias.lower()]
                break
    return mapping

def related_error_columns(df: pd.DataFrame, base_col: str) -> List[str]:
    out = []
    for suf in ERR_SUFFIXES:
        cand = f"{base_col}{suf}"
        if cand in df.columns:
            out.append(cand)
    return out

def detect_constant_q2(df: pd.DataFrame) -> float | None:
    for c in df.columns:
        if c.lower() == "q2":
            uniq = df[c].dropna().unique()
            if len(uniq) == 1:
                return float(uniq[0])
    return None

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
    """
    Return (n, A) for sigma = A * Q^n using two points.
    Returns (nan, nan) if invalid.
    """
    if not (np.isfinite(s1) and np.isfinite(s2) and s1 > 0 and s2 > 0):
        return (np.nan, np.nan)
    if not (q2_1 > 0 and q2_2 > 0 and q2_1 != q2_2):
        return (np.nan, np.nan)
    denom = math.log(q2_2 / q2_1)
    if denom == 0:
        return (np.nan, np.nan)
    n = math.log(s2 / s1) / denom
    # A from each point (they should match ideally)
    A1 = s1 / (q2_1 ** n)
    A2 = s2 / (q2_2 ** n)
    if A1 > 0 and A2 > 0:
        A = math.sqrt(A1 * A2)  # geometric mean (robust)
    else:
        A = np.nan
    return (n, A)

def safe_slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


# ----- Plotting -----
def plot_bin_q2_dependence(
    outdir: Path,
    sf_name: str,
    bin_label: str,
    q2_1: float,
    q2_2: float,
    s1: float,
    s2: float,
    n: float,
    A: float,
    q2_ref: float,
    yerr1: float | None = None,
    yerr2: float | None = None,
):
    """
    Make a per-bin plot for a given structure function:
      - data points at (Q2_1, s1), (Q2_2, s2)
      - fitted curve sigma(Q2) = A * Q^n
      - vertical line at Q2_ref, annotated with predicted sigma
    """
    if not (np.isfinite(n) and np.isfinite(A)):
        return  # skip bad fits

    q2_min = min(q2_1, q2_2, q2_ref)
    q2_max = max(q2_1, q2_2, q2_ref)
    x = np.linspace(0.9*q2_min, 1.1*q2_max, 200)
    x = x[x > 0]
    y = A * (x ** n)
    y_ref = A * (q2_ref ** n)

    plt.figure()
    # fitted curve
    plt.plot(x, y, label=f"fit: σ(Q²) = A·Q²^n")
    # data points
    if yerr1 is not None and np.isfinite(yerr1) and yerr1 >= 0:
        plt.errorbar([q2_1], [s1], yerr=[yerr1], fmt="o", label=f"data @ Q²={q2_1:g}")
    else:
        plt.plot([q2_1], [s1], "o", label=f"data @ Q²={q2_1:g}")
    if yerr2 is not None and np.isfinite(yerr2) and yerr2 >= 0:
        plt.errorbar([q2_2], [s2], yerr=[yerr2], fmt="s", label=f"data @ Q²={q2_2:g}")
    else:
        plt.plot([q2_2], [s2], "s", label=f"data @ Q²={q2_2:g}")

    # Q2_ref indicator
    plt.axvline(q2_ref, linestyle="--", label=f"Q²_ref = {q2_ref:g}")
    plt.title(f"{sf_name} — {bin_label}")
    plt.xlabel("Q² (GeV²)")
    plt.ylabel(f"σ_{sf_name} (arb. units)")
    # parameterization text
    txt = f"Fit params:\n n = {n:.4f}\n A = {A:.4g}\n σ(Q²_ref) = {y_ref:.4g}"
    plt.annotate(txt, xy=(0.02, 0.98), xycoords="axes fraction", va="top", ha="left")
    plt.legend()
    fname = f"{safe_slug(sf_name)}__{safe_slug(bin_label)}.png"
    outpath = outdir / sf_name / fname
    outpath.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


# ----- Main -----
def main():
    ap = argparse.ArgumentParser(description="Fit Q^2 dependence for L,T,LT,TT from two CSVs and plot.")
    ap.add_argument("--file1", type=Path, required=True, help="CSV for dataset 1 (e.g., Q²≈3)")
    ap.add_argument("--file2", type=Path, required=True, help="CSV for dataset 2 (e.g., Q²≈4.4)")
    ap.add_argument("--outdir", type=Path, default="", help="Output directory for plots & params.json")
    ap.add_argument("--q2_ref", type=float, default=4.3, help="Reference Q² to annotate in plots")
    ap.add_argument("--left_on", type=str, default="", help="Comma-separated bin keys in file1 (optional)")
    ap.add_argument("--right_on", type=str, default="", help="Comma-separated bin keys in file2 (optional)")
    args = ap.parse_args()

    df1 = pd.read_csv(args.file1)
    df2 = pd.read_csv(args.file2)

    # Validate constant Q2 in each file
    q2_1 = detect_constant_q2(df1)
    q2_2 = detect_constant_q2(df2)
    if q2_1 is None or q2_2 is None:
        raise ValueError("Each CSV must contain a constant 'Q2' column with a single unique value.")

    # Identify SF columns
    sf1 = find_sf_columns(df1)
    sf2 = find_sf_columns(df2)
    common_sfs = sorted(set(sf1.keys()) & set(sf2.keys()))
    if not common_sfs:
        raise ValueError(f"No common structure functions detected. file1: {sf1}  file2: {sf2}")

    # Choose merge keys
    if args.left_on.strip() and args.right_on.strip():
        left_keys = [s.strip() for s in args.left_on.split(",") if s.strip()]
        right_keys = [s.strip() for s in args.right_on.split(",") if s.strip()]
    else:
        left_keys = pick_binning_keys(df1, sf1)
        right_keys = pick_binning_keys(df2, sf2)

    # Slim to relevant columns
    def cols_to_keep(df: pd.DataFrame, sf_cols_map: Dict[str, str], keys: List[str]) -> List[str]:
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
    if merged.empty:
        raise ValueError(
            "No overlapping bins after merge. Provide explicit --left_on/--right_on with matching keys."
        )

    # Prepare parameters store
    params_summary = {sf: [] for sf in common_sfs}

    # Loop over bins (rows) and SFs
    # Build a human-readable bin label from key columns
    def bin_label_from_row(row: pd.Series, keys: List[str]) -> str:
        parts = []
        for k in keys:
            val = row[k]
            if isinstance(val, float):
                parts.append(f"{k}={val:.4g}")
            else:
                parts.append(f"{k}={val}")
        return ", ".join(parts) if parts else "bin"

    for idx, row in merged.iterrows():
        bin_label = bin_label_from_row(row, left_keys)  # left/right keys correspond after merge

        for sf in common_sfs:
            c1 = sf1[sf] + "_f1"
            c2 = sf2[sf] + "_f2"
            s1 = float(row[c1])
            s2 = float(row[c2])

            # Optional y-errors: prefer a single combined *_err, else quadrature of stat/sys if present
            def grab_err(prefix: str) -> float | None:
                base = sf1[sf] if prefix == "_f1" else sf2[sf]
                ecols = related_error_columns(df1 if prefix == "_f1" else df2, base)
                # pick a single '_err' if present
                for nm in ecols:
                    if nm.endswith("_err"):
                        colname = nm + prefix
                        if colname in merged.columns:
                            return float(row[colname])
                # else combine stat & sys if present
                stat = sys = None
                if (base + "_stat" + prefix) in merged.columns:
                    stat = float(row[base + "_stat" + prefix])
                if (base + "_sys" + prefix) in merged.columns:
                    sys = float(row[base + "_sys" + prefix])
                if stat is not None or sys is not None:
                    stat = 0.0 if stat is None else stat
                    sys = 0.0 if sys is None else sys
                    return math.sqrt(stat**2 + sys**2)
                return None

            yerr1 = grab_err("_f1")
            yerr2 = grab_err("_f2")

            n, A = compute_n_and_A(s1, s2, q2_1, q2_2)
            if not (np.isfinite(n) and np.isfinite(A)):
                continue

            # Save parameters
            params_summary[sf].append({
                "bin": bin_label,
                "Q2_1": q2_1,
                "Q2_2": q2_2,
                "A": A,
                "n": n,
                "sigma_at_Q2ref": float(A * (args.q2_ref ** n)),
            })

            # Plot
            plot_bin_q2_dependence(
                outdir=args.outdir,
                sf_name=sf,
                bin_label=bin_label,
                q2_1=q2_1, q2_2=q2_2,
                s1=s1, s2=s2,
                n=n, A=A,
                q2_ref=args.q2_ref,
                yerr1=yerr1, yerr2=yerr2
            )

    # Write parameters JSON
    args.outdir.mkdir(parents=True, exist_ok=True)
    with open(args.outdir / "q2_fit_params.json", "w") as f:
        json.dump(params_summary, f, indent=2)

    # Console summary
    print("Fitted structure functions:", ", ".join(common_sfs))
    print("Plots written under:", args.outdir.resolve())
    print("Parameters JSON:", (args.outdir / "q2_fit_params.json").resolve())


if __name__ == "__main__":
    main()