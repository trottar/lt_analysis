#!/usr/bin/env python3
"""
q2_only_global.py — One Q² exponent per setting (L, T, LT, TT), collapse both datasets to Q²_ref,
and plot σ_scaled(t) overlays to verify alignment. Saves ONE PDF and ONE JSON in outdir.

Usage:
  python q2_only_global.py \
    --file1 KaonLT_Q3p0W2p32.csv \
    --file2 KaonLT_Q4p4W2p74.csv \
    --outdir q2_out \
    --q2_ref 4.3

What it does
------------
- Determine effective Q² for each file.
- For each setting s∈{L,T,LT,TT}, compute per-bin slopes
    n_{s,j} = [ln σ2 - ln σ1] / [ln Q²2 - ln Q²1]
  with uncertainty from error propagation (dsigma/ sigma).
- Combine across bins: n_s = weighted mean(n_{s,j}) with inverse-variance weights.
- Scale both datasets to Q²_ref with that single n_s.
- Plot σ_scaled(t) overlays (dataset 1 vs dataset 2) — one page per setting — in a single PDF.
- Save JSON { setting: {n, dn, Q2_file1_eff, Q2_file2_eff, Q2_ref} }.

Notes
-----
- This script *does not* fit t-dependence; that’s handled elsewhere.
- It expects SF columns named like sigL, sigT, sigLT, sigTT (aliases supported).
- Errors: prefers dsig* columns if present.
"""

from __future__ import annotations
import argparse, math, json
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# ---------- Config ----------
SF_ALIASES = {
    "L":  ["sigL", "sigmaL", "sigma_L", "SL", "sL", "L"],
    "T":  ["sigT", "sigmaT", "sigma_T", "ST", "sT", "T"],
    "LT": ["sigLT", "sigmaLT", "sigma_LT", "S_LT", "sLT", "LT"],
    "TT": ["sigTT", "sigmaTT", "sigma_TT", "S_TT", "sTT", "TT"],
}
NON_BIN_TOKENS = set(["q2","Q2","w","W","epsilon","eps","Epsilon","dQ2","dW","chi2"])

# ---------- Helpers ----------
def effective_q2(df: pd.DataFrame) -> Optional[float]:
    q2_col = next((c for c in df.columns if c.lower()=="q2"), None)
    if q2_col is None:
        return None
    q2 = pd.to_numeric(df[q2_col], errors="coerce").to_numpy()
    q2 = q2[np.isfinite(q2)]
    if q2.size == 0:
        return None
    if np.unique(np.round(q2,6)).size == 1:
        return float(q2[0])
    if "dQ2" in df.columns:
        dq = pd.to_numeric(df["dQ2"], errors="coerce").to_numpy()
        if dq.size == q2.size:
            mask = np.isfinite(dq) & (dq>0)
            w = np.zeros_like(dq, dtype=float)
            w[mask] = 1.0/(dq[mask]**2)
            if np.any(w>0):
                return float(np.average(q2, weights=w))
    return float(np.median(q2))

def find_sf_columns(df: pd.DataFrame) -> Dict[str,str]:
    cols = list(df.columns)
    lowers = {c.lower():c for c in cols}
    out = {}
    for key, aliases in SF_ALIASES.items():
        # prefer specific long names over single letters
        chosen = None
        for a in aliases:
            if a in cols:
                chosen = a; break
        if chosen is None:
            for a in aliases:
                al = a.lower()
                if al in lowers:
                    cand = lowers[al]
                    if key=="T" and cand.lower()=="t":  # avoid mapping to Mandelstam t
                        continue
                    chosen = cand; break
        if chosen: out[key] = chosen
    return out

def dsig_column_name(base: str, suffix: str) -> Optional[str]:
    # prefer dsig* pattern
    cand = "d"+base+suffix
    return cand

def pick_merge_strategy(df1: pd.DataFrame, df2: pd.DataFrame) -> Tuple[str, float]:
    """
    Return merge key and tolerance. Prefer 'tbin' if in both.
    If 't' is in both, we later use a nearest-neighbor match with tolerance.
    """
    if "tbin" in df1.columns and "tbin" in df2.columns:
        return ("tbin", 0.0)
    if "t" in df1.columns and "t" in df2.columns:
        return ("t", 1e-6)  # tolerance used later
    # fallback: intersection of numeric columns minus obvious kinematics / SFs
    return ("tbin", 0.0)  # safest default; user can override upstream if needed

def nearest_merge_on_t(df1: pd.DataFrame, df2: pd.DataFrame, tol: float) -> pd.DataFrame:
    # simple nearest merge on 't' within tol
    a = df1.copy(); b = df2.copy()
    a = a.sort_values("t").reset_index(drop=True)
    b = b.sort_values("t").reset_index(drop=True)
    idx2 = 0; rows = []
    for i, r1 in a.iterrows():
        t1 = float(r1["t"])
        # advance idx2 until b.t is closest
        while idx2+1 < len(b) and abs(float(b.loc[idx2+1,"t"])-t1) <= abs(float(b.loc[idx2,"t"])-t1):
            idx2 += 1
        t2 = float(b.loc[idx2,"t"])
        if abs(t2 - t1) <= tol:
            rows.append((r1, b.loc[idx2]))
    if not rows:
        return pd.DataFrame()
    left = pd.DataFrame([r[0] for r in rows]).reset_index(drop=True)
    right = pd.DataFrame([r[1] for r in rows]).reset_index(drop=True)
    merged = pd.concat([left.add_suffix("_f1"), right.add_suffix("_f2")], axis=1)
    return merged

def ln_err_from_abs(y: float, dy: Optional[float]) -> Optional[float]:
    if dy is None or not np.isfinite(dy) or dy<=0 or y<=0:
        return None
    return abs(dy/y)

def combine_bin_slopes(slopes: List[Tuple[float,float]]) -> Tuple[float,float]:
    ws, xs = [], []
    for n, dn in slopes:
        if np.isfinite(n) and np.isfinite(dn) and dn>0:
            ws.append(1.0/(dn*dn)); xs.append(n)
    if not ws:
        return (np.nan, np.nan)
    wsum = float(np.sum(ws))
    nbar = float(np.sum(np.array(ws)*np.array(xs))/wsum)
    dnbar = math.sqrt(1.0/wsum)
    return (nbar, dnbar)

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(description="One global Q² exponent per setting; collapse to Q²_ref.")
    ap.add_argument("--file1", type=Path, required=True)
    ap.add_argument("--file2", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, default="")
    ap.add_argument("--q2_ref", type=float, default=4.3)
    args = ap.parse_args()

    df1 = pd.read_csv(args.file1)
    df2 = pd.read_csv(args.file2)

    q2_1 = effective_q2(df1); q2_2 = effective_q2(df2)
    if q2_1 is None or q2_2 is None:
        raise ValueError("Could not determine effective Q2 from one or both files (need a 'Q2' column).")

    sf1 = find_sf_columns(df1); sf2 = find_sf_columns(df2)
    settings = [s for s in ["L","T","LT","TT"] if (s in sf1 and s in sf2)]
    if not settings:
        raise ValueError(f"No common structure functions detected. file1: {sf1}  file2: {sf2}")

    # Merge strategy
    key, tol = pick_merge_strategy(df1, df2)
    if key == "tbin":
        keep1 = [*df1.columns]; keep2 = [*df2.columns]
        # inner merge on tbin
        merged = pd.merge(df1, df2, on="tbin", how="inner", suffixes=("_f1","_f2"))
    else:
        merged = nearest_merge_on_t(df1, df2, tol)
    if merged.empty:
        raise ValueError("No overlapping bins. Ensure common 'tbin' or compatible 't' ranges.")

    # Prepare outputs
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    pdf_path = outdir / "plots.pdf"
    json_path = outdir / "global_q2_params.json"

    # Compute global n_s per setting
    results = {}
    with PdfPages(pdf_path) as pdf:
        for s in settings:
            base1 = sf1[s]; base2 = sf2[s]
            ds1_col = "d"+base1+"_f1" if ("d"+base1+"_f1" in merged.columns) else None
            ds2_col = "d"+base2+"_f2" if ("d"+base2+"_f2" in merged.columns) else None

            # per-bin slopes
            slopes: List[Tuple[float,float]] = []
            for _, row in merged.iterrows():
                sig1 = float(row[base1+"_f1"]); sig2 = float(row[base2+"_f2"])
                if not (sig1>0 and sig2>0):
                    continue
                n_j = (math.log(sig2) - math.log(sig1)) / (math.log(q2_2) - math.log(q2_1))
                # error on ln sigma from absolute errors
                dy1 = float(row[ds1_col]) if (ds1_col and pd.notna(row[ds1_col])) else None
                dy2 = float(row[ds2_col]) if (ds2_col and pd.notna(row[ds2_col])) else None
                e1 = ln_err_from_abs(sig1, dy1); e2 = ln_err_from_abs(sig2, dy2)
                var_dln = (e1 or 0.0)**2 + (e2 or 0.0)**2
                dn_j = math.sqrt(var_dln) / abs(math.log(q2_2) - math.log(q2_1)) if var_dln>0 else np.nan
                slopes.append((n_j, dn_j))

            n_s, dn_s = combine_bin_slopes(slopes)
            results[s] = {
                "n": n_s, "dn": dn_s,
                "Q2_file1_eff": q2_1, "Q2_file2_eff": q2_2, "Q2_ref": args.q2_ref
            }

            # Scale both datasets to Q2_ref using single n_s
            # Build scaled tables for plotting vs t/tbin
            def scale_df(df: pd.DataFrame, base: str, q2_src: float, tag: str):
                out = pd.DataFrame()
                if key=="tbin" and "tbin" in df.columns: out["tbin"] = df["tbin"]
                if "t" in df.columns: out["t"] = df["t"]
                if base in df.columns:
                    out["sigma_scaled"] = df[base] * (args.q2_ref / q2_src) ** n_s
                else:
                    raise ValueError(f"Missing column {base} in {tag}")
                # carry error if present
                dcol = "d"+base
                if dcol in df.columns:
                    out["dsigma_scaled"] = df[dcol] * (args.q2_ref / q2_src) ** n_s
                return out

            df1_scaled = scale_df(df1, base1, q2_1, "file1")
            df2_scaled = scale_df(df2, base2, q2_2, "file2")

            # Merge the scaled ones for overlay by tbin (or t)
            if key=="tbin" and "tbin" in df1_scaled.columns and "tbin" in df2_scaled.columns:
                overlay = pd.merge(df1_scaled.add_suffix("_f1"),
                                   df2_scaled.add_suffix("_f2"),
                                   left_on="tbin_f1", right_on="tbin_f2", how="inner")
                # choose plotting x = t if available in either
                x_label = "tbin"
                x_f1 = overlay["tbin_f1"].to_numpy()
                y1 = overlay["sigma_scaled_f1"].to_numpy()
                y2 = overlay["sigma_scaled_f2"].to_numpy()
                dy1 = overlay.get("dsigma_scaled_f1", pd.Series([np.nan]*len(y1))).to_numpy()
                dy2 = overlay.get("dsigma_scaled_f2", pd.Series([np.nan]*len(y2))).to_numpy()
                # try to pick a representative t array for x if present
                if "t_f1" in merged.columns:
                    # map tbin to t via merged
                    map_t = merged[["tbin", "t_f1"]].drop_duplicates().set_index("tbin")["t_f1"]
                    if set(overlay["tbin_f1"]).issubset(set(map_t.index)):
                        x_vals = np.array([map_t[tb] for tb in x_f1])
                        x_label = "t (GeV²)"
                    else:
                        x_vals = x_f1
                elif "t_f2" in merged.columns:
                    map_t = merged[["tbin", "t_f2"]].drop_duplicates().set_index("tbin")["t_f2"]
                    if set(overlay["tbin_f1"]).issubset(set(map_t.index)):
                        x_vals = np.array([map_t[tb] for tb in x_f1])
                        x_label = "t (GeV²)"
                    else:
                        x_vals = x_f1
                else:
                    x_vals = x_f1
            else:
                # fall back to nearest on t
                # reindex on t
                a = df1_scaled.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)
                b = df2_scaled.dropna(subset=["t"]).sort_values("t").reset_index(drop=True)
                m = min(len(a), len(b))
                a = a.iloc[:m]; b = b.iloc[:m]
                x_vals = a["t"].to_numpy()
                y1 = a["sigma_scaled"].to_numpy()
                y2 = b["sigma_scaled"].to_numpy()
                dy1 = a.get("dsigma_scaled", pd.Series([np.nan]*len(y1))).to_numpy()
                dy2 = b.get("dsigma_scaled", pd.Series([np.nan]*len(y2))).to_numpy()
                x_label = "t (GeV²)"

            # Plot page for this setting: σ_scaled(t), both datasets
            plt.figure()
            plt.errorbar(x_vals, y1, yerr=dy1 if np.any(np.isfinite(dy1)) else None,
                         fmt="o", label=f"{s} @ Q²={q2_1:.3g} → scaled")
            plt.errorbar(x_vals, y2, yerr=dy2 if np.any(np.isfinite(dy2)) else None,
                         fmt="s", label=f"{s} @ Q²={q2_2:.3g} → scaled")
            plt.title(f"{s}: single global n = {n_s:.3f} ± {dn_s:.3f}  (to Q²_ref={args.q2_ref:g})")
            plt.xlabel(x_label); plt.ylabel(f"σ_{s}(Q²_ref)")
            plt.legend(); plt.tight_layout()
            pdf.savefig(); plt.close()

    # Save JSON
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)

    print("Saved:", pdf_path.resolve())
    print("Saved:", json_path.resolve())

if __name__ == "__main__":
    main()
