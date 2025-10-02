#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Produce ONE file per setting (all t-bins appended). Filenames are defined
EXCLUSIVELY by two REQUIRED CSV columns you provide:

  Q2_token , W_token   (strings like: 16 for 1.6, 222 for 2.22, 44 for 4.4, 274 for 2.74)

Example outputs:
  averages/avek.Q44W274.dat
  xsects/x_sep.pl_Q44W274.dat

INPUT CSV (one row per t-bin) — column names are case-insensitive:
  Required:
    Q2, dQ2, W, dW, t, dt,
    sigL, dsigL, sigT, dsigT, sigLT, dsigLT, sigTT, dsigTT, chi2,
    Q2_token, W_token
  Optional:
    tbin  (if absent, bins will be auto-numbered 1..N per setting)

CLI:
  --ptype {pion,pionlt,kaon,kaonlt}  # sets final-state masses (π⁺n or K⁺Λ) for θ*
  --pol   {pl,mn}                    # tag embedded in x_sep filename

Physics:
  θ* (deg) is computed in the γ*–p CM frame from (Q², W, t). The CSV may
  provide either t (<0) or −t (>0); sign is handled automatically.
"""

import argparse
import math
import os
from typing import Tuple

import pandas as pd

# ---- Physical constants (GeV) ----
MP  = 0.9382720813
MN  = 0.9395654133
MPI = 0.13957039
MK  = 0.493677
ML  = 1.115683  # Lambda^0

# ---- Helpers ----
def kallen(x: float, y: float, z: float) -> float:
    return x*x + y*y + z*z - 2.0*(x*y + x*z + y*z)

def cm_energy_momentum(s: float, m1sq: float, m2sq: float, W: float) -> Tuple[float, float]:
    E = (s + m1sq - m2sq) / (2.0 * W)
    p = (kallen(s, m1sq, m2sq) ** 0.5) / (2.0 * W) if s > 0 else 0.0
    return E, max(0.0, p)

def theta_cm_deg(Q2: float, W: float, t_or_tneg: float, m_mes: float, m_baryon: float) -> float:
    """
    Compute θ* [deg] for γ*(Q2)+p -> meson + baryon from (Q^2, W, t).
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

    den = 2.0 * p1 * p3
    if den == 0.0:
        return 0.0
    num = t - m1sq - m3sq + 2.0 * E1 * E3
    cos_th = max(-1.0, min(1.0, num / den))
    return math.degrees(math.acos(cos_th))

# ---- Core ----
def build_outputs(df: pd.DataFrame, pol_tag: str, ptype: str) -> None:
    # Channel masses
    ptype_l = ptype.lower()
    if ptype_l.startswith("pion"):
        m_mes, m_baryon = MPI, MN      # π+ n
    elif ptype_l.startswith("kaon"):
        m_mes, m_baryon = MK, ML       # K+ Λ
    else:
        raise ValueError("Unknown --ptype (use: pion, pionlt, kaon, kaonlt)")

    # Canonicalize (case-insensitive) column lookup
    canon = {c.lower(): c for c in df.columns}

    # Required columns
    required = {
        "Q2_token","W_token",
        "Q2","dQ2","W","dW","t","dt",
        "sigL","dsigL","sigT","dsigT","sigLT","dsigLT","sigTT","dsigTT","chi2"        
    }
    missing = [k for k in required if k not in canon]
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    def C(name): return canon[name]

    # Ensure tbin labels; if absent, create per setting
    tbin_col = canon.get("tbin")
    if tbin_col is None:
        df["tbin"] = None
        tbin_col = "tbin"

    # Compute θ* for each row
    df["theta_cm_deg"] = [
        theta_cm_deg(
            Q2=row[C("q2")],
            W=row[C("w")],
            t_or_tneg=row[C("t")],
            m_mes=m_mes,
            m_baryon=m_baryon,
        )
        for _, row in df.iterrows()
    ]

    # Prepare output dirs
    os.makedirs("averages", exist_ok=True)
    os.makedirs("xsects", exist_ok=True)

    # Group strictly by the provided tokens → ONE file per setting
    for (qtok, wtok), sub in df.groupby([C("Q2_token"), C("W_token")], sort=False):
        # If tbin absent, enumerate 1..N within this setting
        if canon.get("tbin") is None:
            sub = sub.sort_values(by=[C("t")]).reset_index(drop=True)
            sub["tbin"] = range(1, len(sub) + 1)
        else:
            sub = sub.sort_values(by=[C("t")]).reset_index(drop=True)

        avek_path = os.path.join("averages", f"avek.Q{qtok}W{wtok}.dat")
        x_path    = os.path.join("xsects",  f"x_sep.{pol_tag}_Q{qtok}W{wtok}.dat")

        # --- averages file (all t-bins appended) ---
        with open(avek_path, "w") as f:
            for _, r in sub.iterrows():
                f.write(
                    f"{r[C('w')]:.6f} {r[C('dw')]:.6f} "
                    f"{r[C('q2')]:.6f} {r[C('dq2')]:.6f} "
                    f"{r[C('t')]:.6f} {r[C('dt')]:.6f} "
                    f"{r['theta_cm_deg']:.2f} {int(r['tbin'])}\n"
                )

        # --- cross sections file (all t-bins appended) ---
        with open(x_path, "w") as f:
            for _, r in sub.iterrows():
                f.write(
                    f"{r[C('sigl')]:.6f} {r[C('dsigl')]:.6f} "
                    f"{r[C('sigt')]:.6f} {r[C('dsigt')]:.6f} "
                    f"{r[C('siglt')]:.6f} {r[C('dsiglt')]:.6f} "
                    f"{r[C('sigtt')]:.6f} {r[C('dsigtt')]:.6f} "
                    f"{r[C('chi2')]:.6f} "
                    f"{r[C('t')]:.6f} {r[C('w')]:.6f} {r[C('q2')]:.6f} {r['theta_cm_deg']:.2f}\n"
                )

        print(f"Wrote (all t-bins): {avek_path}")
        print(f"Wrote (all t-bins): {x_path}")

# ---- CLI ----
def main():
    ap = argparse.ArgumentParser(
        description="Produce per-setting files using CSV-provided Q2_token/W_token. One file per setting; all t-bins appended."
    )
    ap.add_argument("csv", help="Path to input CSV with per-bin kinematics and separated cross sections")
    ap.add_argument("--ptype", required=True, choices=["pion", "pionlt", "kaon", "kaonlt"],
                    help="Particle/channel type: sets masses (π+n or K+Λ).")
    ap.add_argument("--pol", required=True, choices=["pl", "mn"],
                    help="Tag used in x_sep.<pol>_Q...W... filename (e.g., pl or mn).")
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    build_outputs(df, pol_tag=args.pol, ptype=args.ptype)

if __name__ == "__main__":
    main()
