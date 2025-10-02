#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Produce ONE file per setting (all t-bins appended), where the setting name used
in the filenames is defined by TWO OPTIONAL CSV COLUMNS supplied by you:

    Q2_token , W_token

If present, these string tokens are used verbatim in the filenames:
    averages/avek.Q<Q2_token>W<W_token>.dat
    xsects/x_sep.<pol>_Q<Q2_token>W<W_token>.dat

If Q2_token / W_token are NOT present, the script falls back to the
"decimal-removed" style derived from the numeric values:
    1.6  -> '16'
    2.22 -> '222'

INPUT CSV (one row per t-bin), case-insensitive column names accepted:
    Required:
      Q2, dQ2, W, dW, t, dt,
      sigL, dsigL, sigT, dsigT, sigLT, dsigLT, sigTT, dsigTT, chi2
    Optional:
      tbin        (int; if absent auto-numbered 1..N per setting)
      Q2_token    (str; filename token override for Q^2)
      W_token     (str; filename token override for W)

OUTPUT (per setting):
  - averages/avek.Q<Q2_token_or_auto>W<W_token_or_auto>.dat
        Columns per line:
            W  dW  Q2  dQ2  t  dt  theta*  t-bin
  - xsects/x_sep.<pol>_Q<Q2_token_or_auto>W<W_token_or_auto>.dat
        Columns per line:
            sigL  dsigL  sigT  dsigT  sigLT  dsigLT  sigTT  dsigTT  chi2  t  W  Q2  theta*

CLI:
  --ptype {pion,pionlt,kaon,kaonlt}   # sets final-state masses (π⁺n or K⁺Λ)
  --pol   {pl,mn}                     # tag embedded in x_sep filename

Physics:
  theta* (deg) is computed in the γ*–p CM frame from (Q^2, W, t).
  t may be provided as t (<0) or −t (>0); sign is handled automatically.
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

    den = 2.0 * p1 * p3
    if den == 0.0:
        return 0.0
    num = t - m1sq - m3sq + 2.0 * E1 * E3
    cos_th = max(-1.0, min(1.0, num / den))
    return math.degrees(math.acos(cos_th))

def fmt_number_token(x: float, decimals: int = 3) -> str:
    """Format a float with up to `decimals` decimals, trimming trailing zeros/dot."""
    s = f"{x:.{decimals}f}".rstrip("0").rstrip(".")
    return s

def token_no_dot_from_float(x: float) -> str:
    """Number token with decimal point removed: 1.6 -> '16', 2.22 -> '222'."""
    return fmt_number_token(x).replace(".", "")

def token_no_dot_from_any(x) -> str:
    """
    Make a safe token from an arbitrary CSV field:
      - If it's a string: strip spaces, replace '.' with nothing.
      - If it's numeric: use token_no_dot_from_float.
    """
    if isinstance(x, str):
        return x.strip().replace(".", "")
    try:
        xf = float(x)
        return token_no_dot_from_float(xf)
    except Exception:
        # Fallback: keep as-is without dots/spaces
        return str(x).strip().replace(".", "")

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

    # Canonicalize column lookup (case-insensitive)
    canon = {c.lower(): c for c in df.columns}

    required = {
        "Q2", "dQ2", "W", "dW", "t", "dt",
        "sigL", "dsigL", "sigT", "dsigT", "sigLT", "dsigLT", "sigTT", "dsigTT", "chi2"
    }
    missing = [k for k in required if k.lower() not in canon]
    if missing:
        raise ValueError(f"Input CSV is missing required columns: {missing}")

    def C(name): return canon[name.lower()]

    # Optional token overrides
    qtok_col = canon.get("q2_token")
    wtok_col = canon.get("w_token")

    # t-bin labels; if absent, create per setting
    tbin_col = canon.get("tbin")
    if tbin_col is None:
        df["tbin"] = None
        tbin_col = "tbin"

    # Compute theta*
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

    # Prepare tokens row-wise (use provided tokens if present, else fallback)
    if qtok_col is None:
        df["_Qtok"] = [token_no_dot_from_float(v) for v in df[C("Q2")]]
    else:
        df["_Qtok"] = [token_no_dot_from_any(v) for v in df[qtok_col]]

    if wtok_col is None:
        df["_Wtok"] = [token_no_dot_from_float(v) for v in df[C("W")]]
    else:
        df["_Wtok"] = [token_no_dot_from_any(v) for v in df[wtok_col]]

    # If tbin missing, populate 1..N per setting (as defined by tokens)
    if canon.get("tbin") is None:
        for (qt, wt), idx in df.groupby(["_Qtok", "_Wtok"]).groups.items():
            df.loc[idx, "tbin"] = range(1, len(idx) + 1)

    # Prepare output dirs
    os.makedirs("averages", exist_ok=True)
    os.makedirs("xsects", exist_ok=True)

    # Group by the *tokens* to define a setting; write ONE file per setting
    for (qtok, wtok), sub in df.groupby(["_Qtok", "_Wtok"], sort=False):
        # Keep each t-bin line; sort by t for readability
        sub = sub.sort_values(by=[C("t")]).reset_index(drop=True)

        avek_path = os.path.join("averages", f"avek.Q{qtok}W{wtok}.dat")
        x_path    = os.path.join("xsects",  f"x_sep.{pol_tag}_Q{qtok}W{wtok}.dat")

        # --- averages file (all t-bins appended) ---
        with open(avek_path, "w") as f:
            for _, r in sub.iterrows():
                f.write(
                    f"{r[C('W')]:.6f} {r[C('dW')]:.6f} "
                    f"{r[C('Q2')]:.6f} {r[C('dQ2')]:.6f} "
                    f"{r[C('t')]:.6f} {r[C('dt')]:.6f} "
                    f"{r['theta_cm_deg']:.6f} {int(r[tbin_col])}\n"
                )

        # --- cross sections file (all t-bins appended) ---
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

        print(f"Wrote (all t-bins): {avek_path}")
        print(f"Wrote (all t-bins): {x_path}")

# --------- CLI ----------
def main():
    ap = argparse.ArgumentParser(
        description="Produce per-setting files using CSV-provided Q2_token/W_token (or numeric fallbacks)."
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
