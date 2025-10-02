#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Build per-(Q2,W) output files from an input CSV of Fpi2-like bins.

INPUT CSV (one row per t-bin) must contain at least these columns:
    Q2, dQ2, W, dW, t, dt,
    sigL, dsigL, sigT, dsigT, sigLT, dsigLT, sigTT, dsigTT, chi2
Optional:
    tbin  (integer label for the t-bin; if absent it will be 1..N)

OUTPUT (created per unique (Q2,W)):
  - averages/avek.Q<Q2val>W<Wval>.dat
        Columns:  W  dW  Q2  dQ2  t  dt  theta*  t-bin
  - xsects/x_sep.<pol>_Q<Q2valstring>W<Wvalstring>.dat
        Columns:  sigL  dsigL  sigT  dsigT  sigLT  dsigLT  sigTT  dsigTT  chi2  t  W  Q2  theta*

Filename tokens:
  - <Q2val>/<Wval>           : decimal point removed (1.6 → "16", 2.22 → "222")
  - <Q2valstring>/<Wvalstring>: '.' replaced by 'p' (1.6 → "1p6", 2.22 → "2p22")

Particle/channel selection:
  --ptype {pion,pionlt,kaon,kaonlt}
    * pion / pionlt :  γ* p → π+ n   (mπ, m_n)
    * kaon / kaonlt :  γ* p → K+ Λ   (mK, m_Λ)
The “lt” suffix is only for your convenience in naming; physics masses are the same.
If --pol is not given, the filename <pol> tag defaults to --ptype.

theta* (deg) is computed in the γ*–p CM frame using two-body kinematics:
    m_γ*^2 = −Q^2, target = proton, final = meson + baryon (per --ptype).
The third argument may be either t (<0) or −t (>0); sign handled internally.

Usage examples:
  python build_products.py fpi2_table.csv --ptype pionlt
  python build_products.py kaon_table.csv --ptype kaonlt --pol sep
"""

import argparse
import math
import os
from typing import Tuple

import pandas as pd

# --------- Physical constants (GeV) ----------
MP  = 0.9382720813
MN  = 0.9395654133
MPI = 0.13957039
MK  = 0.493677
ML  = 1.115683  # Lambda^0

# --------- Helpers ----------
def kallen(x: float, y: float, z: float) -> float:
    return x*x + y*y + z*z - 2.0*(x*y + x*z + y*z)

def cm_energy_momentum(s: float, m1sq: float, m2sq: float, W: float) -> Tuple[float, float]:
    E = (s + m1sq - m2sq) / (2.0 * W)
    p = (kallen(s, m1sq, m2sq) ** 0.5) / (2.0 * W) if s > 0 else 0.0
    return E, max(0.0, p)

def theta_cm_deg(Q2: float, W: float, t_or_tneg: float, m_mes: float, m_baryon: float) -> float:
    """
    Compute theta* [deg] for γ*(Q2)+p -> meson + baryon from (Q^2, W, t).
    Accepts either t (<0) or -t (>0) as the third argument.
    """
    # Interpret input as either t (<0) or -t (>0)
    t = t_or_tneg if t_or_tneg < 0.0 else -abs(t_or_tneg)
    s = W * W

    m1sq = -Q2            # virtual photon
    m2sq = MP * MP        # proton target
    m3sq = m_mes * m_mes
    m4sq = m_baryon * m_baryon

    E1, p1 = cm_energy_momentum(s, m1sq, m2sq, W)
    E3, p3 = cm_energy_momentum(s, m3sq, m4sq, W)

    # t = m1^2 + m3^2 - 2(E1 E3 - p1 p3 cosθ*)
    den = 2.0 * p1 * p3
    if den == 0.0:
        return 0.0
    num = t - m1sq - m3sq + 2.0 * E1 * E3
    cos_th = max(-1.0, min(1.0, num / den))
    return math.degrees(math.acos(cos_th))

def fmt_number_token(x: float, decimals: int = 3) -> str:
    s = f"{x:.{decimals}f}".rstrip("0").rstrip(".")
    return s

def token_no_dot(x: float) -> str:
    return fmt_number_token(x).replace(".", "")

def token_p_dot(x: float) -> str:
    return fmt_number_token(x).replace(".", "p")

# --------- Core ----------
def build_outputs(df: pd.DataFrame, pol_tag: str, ptype: str) -> None:
    # Channel masses
    ptype_l = ptype.lower()
    if ptype_l.startswith("pion"):
        m_mes, m_baryon = MPI, MN      # π+ n
    elif ptype_l.startswith("kaon"):
        m_mes, m_baryon = MK, ML       # K+ Λ
    else:
        raise ValueError("Unknown --ptype (use: pion, pionlt, kaon, kaonlt)")

    # Required columns (case-insensitive)
    required = {
        "Q2", "dQ2", "W", "dW", "t", "dt",
        "sigL", "dsigL", "sigT", "dsigT", "sigLT", "dsigLT", "sigTT", "dsigTT", "chi2"
    }
    # Build lower->actual name map
    canon_map = {c.lower(): c for c in df.columns}
    missing = [c for c in required if c.lower() not in canon_map]
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    # Helper to access with canonical names regardless of input case
    def C(name): return canon_map[name.lower()]

    # Ensure t-bin labels; if absent create per (Q2,W) group
    tbin_col = canon_map.get("tbin")
    if tbin_col is None:
        df["tbin"] = None
        for (_, _), idx in df.groupby([C("Q2"), C("W")]).groups.items():
            df.loc[idx, "tbin"] = range(1, len(idx) + 1)
        tbin_col = "tbin"

    # Compute theta* for each row
    df["theta_cm_deg"] = [
        theta_cm_deg(
            Q2=row[C("Q2")],
            W=row[C("W")],
            t_or_tneg=row[C("t")],
            m_mes=m_mes,
            m_baryon=m_baryon,
        )
        for _, row in df.iterrows()
    ]

    # Prepare output dirs
    os.makedirs("averages", exist_ok=True)
    os.makedirs("xsects", exist_ok=True)

    # Group by (Q2,W) and write files
    for (Q2val, Wval), sub in df.groupby([C("Q2"), C("W")], sort=False):
        sub = sub.sort_values(by=[C("t")]).reset_index(drop=True)

        q_token_num = token_no_dot(Q2val)
        w_token_num = token_no_dot(Wval)
        q_token_str = token_p_dot(Q2val)
        w_token_str = token_p_dot(Wval)

        # --- averages file ---
        avek_path = os.path.join("averages", f"avek.Q{q_token_num}W{w_token_num}.dat")
        with open(avek_path, "w") as f:
            for _, r in sub.iterrows():
                f.write(
                    f"{r[C('W')]:.6f} {r[C('dW')]:.6f} "
                    f"{r[C('Q2')]:.6f} {r[C('dQ2')]:.6f} "
                    f"{r[C('t')]:.6f} {r[C('dt')]:.6f} "
                    f"{r['theta_cm_deg']:.6f} {int(r[tbin_col])}\n"
                )

        # --- cross sections file ---
        pol = pol_tag if pol_tag else ptype_l   # default <pol> tag from --ptype
        x_path = os.path.join("xsects", f"x_sep.{pol}_Q{q_token_str}W{w_token_str}.dat")
        with open(x_path, "w") as f:
            for _, r in sub.iterrows():
                f.write(
                    f"{r[C('sigL')]:.6f} {r[C('dsigL')]:.6f} "
                    f"{r[C('sigT')]:.6f} {r[C('dsigT')]:.6f} "
                    f"{r[C('sigLT')]:.6f} {r[C('dsigLT')]:.6f} "
                    f"{r[C('sigTT')]:.6f} {r[C('dsigTT')]:.6f} "
                    f"{r[C('chi2')]:.6f} "
                    f"{r[C('t')]:.6f} {r[C('W')]:.6f} {r[C('Q2')]:.6f} {r['theta_cm_deg']:.6f}\n"
                )

        print(f"Wrote: {avek_path}")
        print(f"Wrote: {x_path}")

# --------- CLI ----------
def main():
    ap = argparse.ArgumentParser(
        description="Produce averages/xsects files with computed theta* from a Fpi2-style CSV."
    )
    ap.add_argument("csv", help="Path to input CSV with per-bin kinematics and separated cross sections")
    ap.add_argument("--ptype", required=True, choices=["pion", "pionlt", "kaon", "kaonlt"],
                    help="Particle/channel type: sets masses (π+n or K+Λ). "
                         "Filename <pol> tag defaults to this value unless --pol is given.")
    ap.add_argument("--pol", default="pl", choices=["pl", "mn"]
                    help="Optional tag to use in x_sep.<pol>_Q...W... filename (default: use --ptype)")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    build_outputs(df, pol_tag=args.pol, ptype=args.ptype)

if __name__ == "__main__":
    main()
