#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Make per-setting data files from a CSV of t-bins.

REQUIREMENTS (case-insensitive; leading/trailing spaces in headers are OK):
  Required columns:
    Q2, dQ2, W, dW, t, dt,
    sigL, dsigL, sigT, dsigT, sigLT, dsigLT, sigTT, dsigTT, chi2,
    Q2_token, W_token
  Optional:
    tbin  (if absent, bins are auto-numbered 1..N within each setting)

Filenames are defined exclusively by the CSV tokens:
  averages/avek.Q<Q2_token>W<W_token>.dat
  xsects/x_sep.<pol>_Q<Q2_token>W<W_token>.dat

CLI:
  --ptype {pion,pionlt,kaon,kaonlt}   # sets channel masses for θ*
  --pol   {pl,mn}                     # tag embedded in x_sep filename

Physics:
  θ* (deg) computed in the γ*–p CM frame from (Q², W, t).
  t may be provided as t (<0) or −t (>0); the sign is handled internally.
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

    # Normalize column names: strip spaces, lowercase
    df = df.copy()
    df.columns = [c.strip().lower() for c in df.columns]

    # Required columns (lowercase)
    required = [
        "q2","dq2","w","dw","t","dt",
        "sigl","dsigl","sigt","dsigt","siglt","dsiglt","sigtt","dsigtt","chi2",
        "q2_token","w_token"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    # Ensure tbin labels; if absent create per setting
    if "tbin" not in df.columns:
        df["tbin"] = None

    # Compute θ* for each row
    df["theta_cm_deg"] = [
        theta_cm_deg(
            Q2=row["q2"],
            W=row["w"],
            t_or_tneg=row["t"],
            m_mes=m_mes,
            m_baryon=m_baryon,
        )
        for _, row in df.iterrows()
    ]

    # Prepare output dirs
    os.makedirs("averages", exist_ok=True)
    os.makedirs("xsects", exist_ok=True)

    # Group strictly by provided tokens → ONE file per setting (all t-bins)
    for (qtok, wtok), sub in df.groupby(["q2_token", "w_token"], sort=False):
        # Sort by t; auto-enumerate tbin if needed
        sub = sub.sort_values(by=["t"]).reset_index(drop=True)
        if df["tbin"].isna().any():
            sub["tbin"] = range(1, len(sub) + 1)

        avek_path = os.path.join("averages", f"avek.Q{str(qtok).strip()}W{str(wtok).strip()}.dat")
        x_path    = os.path.join("xsects",  f"x_sep.{pol_tag}_Q{str(qtok).strip()}W{str(wtok).strip()}.dat")

        # --- averages file (all t-bins appended) ---
        with open(avek_path, "w") as f:
            for _, r in sub.iterrows():
                f.write(
                    f"{r['w']:.6f} {r['dw']:.6f} "
                    f"{r['q2']:.6f} {r['dq2']:.6f} "
                    f"{r['t']:.6f} {r['dt']:.6f} "
                    f"{r['theta_cm_deg']:.2f} {int(r['tbin'])}\n"
                )

        # --- cross sections file (all t-bins appended) ---
        with open(x_path, "w") as f:
            for _, r in sub.iterrows():
                f.write(
                    f"{r['sigl']:.6f} {r['dsigl']:.6f} "
                    f"{r['sigt']:.6f} {r['dsigt']:.6f} "
                    f"{r['siglt']:.6f} {r['dsiglt']:.6f} "
                    f"{r['sigtt']:.6f} {r['dsigtt']:.6f} "
                    f"{r['chi2']:.6f} "
                    f"{r['t']:.6f} {r['w']:.6f} {r['q2']:.6f} {r['theta_cm_deg']:.2f}\n"
                )

        print(f"Wrote (all t-bins): {avek_path}")
        print(f"Wrote (all t-bins): {x_path}")

# ---- CLI ----
def main():
    import sys
    ap = argparse.ArgumentParser(
        description="Produce per-setting files using CSV-provided Q2_token/W_token (robust to header spacing/case)."
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
