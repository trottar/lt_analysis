#! /usr/bin/python
# -*- coding: utf-8 -*-

#
# Description:
"""
Read two parameter files (new_par.pl_Q??W??.dat), parse par1..par16 + errors + chi2,
parse Q^2 and W from the path, then plot a 4-panel comparison:
  - Set A
  - Set B
  - Avg (weighted)  — original behavior (fixed-effects w = (1/chi2)*(1/err^2))
  - Best (selected) — NEW: per-function selection (FE vs RE) with diagnostics

We quantify per-parameter and per-function agreement (L, T, LT, TT):
  • per-parameter pulls z_k, tau^2 (random-effects)
  • per-function summary: max|z|, mean I^2, model ΔRMS across t
Decision per function:
  - if GOOD  → use Fixed-Effects (FE) for that function’s parameters
  - if not   → use Random-Effects (RE) for that function’s parameters
  - if BAD   → flagged in output as likely requiring functional-form adaptation

All original functionality is preserved:
  • weighted-average output file is still written
  • 4-panel plot is produced

NEW files saved:
  • new_par.<file_str>_best.dat              (per-function FE/RE mixed “best”)
  • new_par.<file_str>_diagnostics.csv       (per-parameter + per-function table)
"""

# ================================================================
# Time-stamp: "2025-10-03"
# ================================================================

import sys, os, re, math, csv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

##################################################################################################################################################

# ---------- ltsep paths ----------
from ltsep import Root
lt = Root(os.path.realpath(__file__), "Plot_LTSep")
OUTPATH   = lt.OUTPATH
LTANAPATH = lt.LTANAPATH
TEMP_CACHEPATH = f"{OUTPATH}/testing_env"

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import load_equations, prepare_equations

##################################################################################################################################################

# ---------- t-range ----------
TMIN = 0.02
TMAX = 0.60

# Output filenames
file_str           = "pl_Q1p6-2p4W2p22"
AVG_OUTFILE        = f"{LTANAPATH}/testing/parameters/new_par.{file_str}.dat"               # original fixed-effects output (preserved)
BEST_OUTFILE       = f"{LTANAPATH}/testing/parameters/new_par.{file_str}_best.dat"          # NEW: per-function FE/RE mix
DIAG_CSV_OUTFILE   = f"{LTANAPATH}/testing/parameters/new_par.{file_str}_diagnostics.csv"   # NEW: diagnostics
AVG_ROW_CHI2       = 3.0  # default chi2 for every parameter row in saved files

# Load equations from model of given setting
equations = load_equations("Q1p6W2p22.model")

# ========================= Tunables: give the model more leeway =========================
TUNABLES = {
    # If parameter tension is tiny, call it FE regardless of shape metric (protects small-signal cases)
    "Z_SMALL":   1.8,   # was 1.5
    "I2_SMALL":  0.15,  # was 0.10

    # FE region (GOOD): be more permissive
    "Z_FE":      2.5,   # was 2.0
    "I2_FE":     0.35,  # was 0.25
    "DRMS_FE":   0.50,  # was 0.35

    # RE region (OK): still acceptable, just inflate errors
    "Z_RE":      4.0,   # was 3.0
    "I2_RE":     0.75,  # was 0.50
    "DRMS_RE":   0.90,  # was 0.60

    # Robust ΔRMS settings
    "TAU_MASK":  1e-3,  # mask t-points where both predictions are below tau * scale
    "EPS_FRAC":  1e-6,  # additive stabilizer fraction of signal scale

    "CHI2_DATA_FE_MAX": 1.25,  # treat function as FE if its data χ² (from rows) ≤ this
}

# ======================= Function → parameter index mapping (explicit) =======================
# Indices are 0-based here. Adjust these lists if you change functional forms.
# Mapped based on use in sigma_* below:
#   L uses par1..par4     → idx 0..3
#   T uses par5..par6     → idx 4..7
#   LT uses par9..par12   → idx 8..11
#   TT uses par13         → idx 12..15
# Any other parameters (e.g., 6–7, 13–16 depending on your model) are treated as "OTHER".
PARAM_GROUPS = {
    "L":  [0,1,2,3],
    "T":  [4,5,6,7],
    "LT": [8,9,10,11],
    "TT": [12,13,14,15],
}
ALL_USED = sorted(set(sum(PARAM_GROUPS.values(), [])))
OTHER_IDX = [i for i in range(16) if i not in ALL_USED]  # keep behavior, but mark as OTHER group
if OTHER_IDX:
    PARAM_GROUPS["OTHER"] = OTHER_IDX  # will follow FE by default unless RE is clearly better

# ---------- Constants ----------
mtar  = 0.93827231
mpipl = 0.139570

# ---------- Physics pieces ----------
wfactor = prepare_equations(equations, 'wfactor')
sigma_L = prepare_equations(equations, 'sig_L')
sigma_T = prepare_equations(equations, 'sig_T')
sigma_LT = prepare_equations(equations, 'sig_LT')
sigma_TT = prepare_equations(equations, 'sig_TT')

# ---------- Group utilities ----------
def function_slices(pars, group_name):
    """Return a view (copy) of parameters for a given function group."""
    idxs = PARAM_GROUPS[group_name]
    return np.array([pars[i] for i in idxs], dtype=float), idxs

def inject_group(pars_base, group_name, group_vals):
    """Return a copy of pars_base with group_name indices replaced by group_vals."""
    out = np.array(pars_base, dtype=float)
    for j, i in enumerate(PARAM_GROUPS[group_name]):
        out[i] = group_vals[j]
    return out

# ---------- Helpers ----------
def parse_Q2_W_from_path(p: Path):
    s = str(p)
    m = re.search(r'Q(\d+)p(\d+)W(\d+)p(\d+)', s)
    if m:
        q2 = float(f"{m.group(1)}.{m.group(2)}")
        w  = float(f"{m.group(3)}.{m.group(4)}")
        return q2, w
    m = re.search(r'_Q(\d+)W(\d+)\b', s)
    if m:
        q2 = float(m.group(1)) / 10.0
        w  = float(m.group(2)) / 100.0
        return q2, w
    raise ValueError(f"Could not parse Q^2 and W from path: {p}")

def read_params_file_rows(par_path: Path):
    """
    Returns:
      vals:      (16,) parameter central values
      errs:      (16,) parameter 1σ errors (second column)
      chi2_g:    scalar (global) = median of row χ² (last column)
      chi2_rows: (16,) per-parameter row χ² (order 1..16)
    """
    values, errors, chi2_rows = {}, {}, {}
    chi2_list = []
    with open(par_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                val   = float(parts[0])
                err   = float(parts[1])
                idx   = int(float(parts[2]))
                chi2v = float(parts[3])   # last column
            except ValueError:
                continue
            values[idx]    = val
            errors[idx]    = err
            chi2_rows[idx] = chi2v
            chi2_list.append(chi2v)

    vals       = np.array([values.get(i, 0.0) for i in range(1, 17)], dtype=float)
    errs       = np.array([errors.get(i, 1.0) for i in range(1, 17)], dtype=float)
    chi2_rowsV = np.array([chi2_rows.get(i, 1.0) for i in range(1, 17)], dtype=float)
    chi2_g     = float(np.median(chi2_list)) if chi2_list else 1.0

    errs   = np.clip(errs, 1e-18, None)
    chi2_g = max(chi2_g, 1e-18)
    return vals, errs, chi2_g, chi2_rowsV

def read_params_file(par_path: Path):
    values, errors, chi2_list = {}, {}, []
    with open(par_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                val   = float(parts[0])
                err   = float(parts[1])
                idx   = int(float(parts[2]))
                chi2v = float(parts[3])
            except ValueError:
                continue
            values[idx] = val
            errors[idx] = err
            chi2_list.append(chi2v)

    vals  = np.array([values.get(i, 0.0) for i in range(1, 17)], dtype=float)
    errs  = np.array([errors.get(i, 1.0) for i in range(1, 17)], dtype=float)
    chi2  = float(np.median(chi2_list)) if chi2_list else 1.0
    errs = np.clip(errs, 1e-18, None)
    chi2 = max(chi2, 1e-18)
    return vals, errs, chi2

def write_params_file(filename, values, errors, chi2_per_row=3.0):
    with open(filename, "w") as f:
        for idx in range(1, 17):
            v = values[idx-1]
            e = errors[idx-1]
            line = f"{v: .8e}   {e: .8e}   {idx:2d}          {chi2_per_row:.1f}\n"
            f.write(line)

def weighted_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B):
    wA = (1.0/chi2A) * (1.0/np.square(errsA))
    wB = (1.0/chi2B) * (1.0/np.square(errsB))
    wsum = np.clip(wA + wB, 1e-18, None)
    pbar = (wA*valsA + wB*valsB) / wsum
    sbar = 1.0 / np.sqrt(wsum)
    return pbar, sbar

# ---------- Plotting ----------
def make_four_panel(datasets, npts=1200, theta_cm=math.pi/2):
    obs = {
        "L":  (sigma_L,  r'$\sigma_L$'),
        "T":  (sigma_T,  r'$\sigma_T$'),
        "LT": (sigma_LT, r'$\sigma_{LT}$'),
        "TT": (sigma_TT, r'$\sigma_{TT}$'),
    }
    t_neg = -np.linspace(TMIN, TMAX, npts)
    tpos  = -t_neg

    fig = plt.figure(figsize=(10, 8))
    for i, (key, (func, ylabel)) in enumerate(obs.items(), start=1):
        ax = fig.add_subplot(2, 2, i)
        for data in datasets:
            y = func(data["pars"], data["qq"], t_neg, theta_cm, data["ww"])
            ax.plot(tpos, y, label=f'{data["name"]}  (Q^2={data["qq"]:.3g}, W={data["ww"]:.3g})')
        ax.set_xlabel(r'$-t\ \mathrm{[GeV^2]}$')
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.grid(True)
        if i == 1:
            ax.legend()
    fig.suptitle("Comparison of $\\sigma_L,\\sigma_T,\\sigma_{LT},\\sigma_{TT}$")
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(f"{TEMP_CACHEPATH}/fits_{file_str}.pdf")
    plt.close(fig)

# ========================= Diagnostics + robust combination =========================
def param_pulls(valsA, errsA, valsB, errsB):
    """Per-parameter pulls: z_k = (A - B) / sqrt(sA^2 + sB^2)."""
    denom = np.sqrt(np.square(errsA) + np.square(errsB))
    denom = np.clip(denom, 1e-18, None)
    return (valsA - valsB) / denom

def dersimonian_laird_tau2(vA, sA, vB, sB, chi2A=None, chi2B=None):
    """
    Method-of-moments tau^2 per parameter. If chi2A/B provided, we honor the script's
    1/chi2 down-weight by inflating within-fit variances before computing tau^2.
    """
    if chi2A is not None and chi2B is not None:
        sA2 = np.square(sA) * chi2A
        sB2 = np.square(sB) * chi2B
    else:
        sA2 = np.square(sA)
        sB2 = np.square(sB)
    wA = 1.0 / np.clip(sA2, 1e-18, None)
    wB = 1.0 / np.clip(sB2, 1e-18, None)
    wsum = np.clip(wA + wB, 1e-18, None)
    m = (wA * vA + wB * vB) / wsum
    Q = wA * np.square(vA - m) + wB * np.square(vB - m)
    denom = wsum - (np.square(wA) + np.square(wB)) / wsum
    tau2 = np.maximum((Q - 1.0) / np.clip(denom, 1e-18, None), 0.0)
    return tau2

def random_effects_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B):
    pulls = param_pulls(valsA, errsA, valsB, errsB)
    tau2 = dersimonian_laird_tau2(valsA, errsA, valsB, errsB, chi2A, chi2B)
    varA = np.square(errsA) * chi2A + tau2
    varB = np.square(errsB) * chi2B + tau2
    wA = 1.0 / np.clip(varA, 1e-18, None)
    wB = 1.0 / np.clip(varB, 1e-18, None)
    wsum = np.clip(wA + wB, 1e-18, None)
    pbar = (wA * valsA + wB * valsB) / wsum
    sbar = np.sqrt(1.0 / wsum)
    return pbar, sbar, pulls, tau2

def global_heterogeneity_stats(pulls):
    """Two-study case per parameter: Q_k = z_k^2 (df=1)."""
    Q_per = np.square(pulls)
    Q_tot = float(np.sum(Q_per))
    df_tot = len(pulls)
    I2_per = np.maximum(0.0, (Q_per - 1.0) / np.clip(Q_per, 1e-18, None))
    I2_mean = float(np.mean(I2_per))
    pval = None
    try:
        from scipy.stats import chi2 as _chi2
        pval = float(_chi2.sf(Q_tot, df_tot))
    except Exception:
        pval = None
    return Q_per, Q_tot, df_tot, I2_per, I2_mean, pval

def model_delta_rms(per_func, parsA, parsB, q2A, q2B, w, npts=400, theta_cm=math.pi/2):
    """
    Robust, scale-aware shape disagreement:
      dRMS = 0.5 * [ NRMSE(q2A) + NRMSE(q2B) ]
    where NRMSE(x) = sqrt(mean((A-B)^2 over valid t)) / (RMS( (|A|+|B|)/2 over valid t ) + eps)

    - Masks t-points where both predictions are below a scale-dependent floor.
    - Uses an additive stabilizer eps to avoid division blowups.
    - Returns np.nan for non-σ groups (e.g., "OTHER").
    """
    funcs = {"L": sigma_L, "T": sigma_T, "LT": sigma_LT, "TT": sigma_TT}
    if per_func not in funcs:
        return float("nan")

    f = funcs[per_func]
    t = -np.linspace(TMIN, TMAX, npts)
    tau = float(TUNABLES["TAU_MASK"])
    eps_frac = float(TUNABLES["EPS_FRAC"])

    def nrmse_at_q2(q2):
        A = f(1.6, 2.22, q2, w, t, theta_cm, parsA)
        B = f(2.45, 2.22, q2, w, t, theta_cm, parsB)
        avg_mag = 0.5 * (np.abs(A) + np.abs(B))
        # robust scale for masking + denominator scaling
        ref = np.percentile(avg_mag, 75)  # upper-mid scale
        mask = avg_mag > (tau * max(ref, 1e-18))
        if not np.any(mask):
            return 0.0  # everything is tiny → no meaningful shape tension

        A_m, B_m = A[mask], B[mask]
        avg_mag_m = avg_mag[mask]

        num = np.sqrt(np.mean((A_m - B_m)**2))
        denom_scale = np.sqrt(np.mean(avg_mag_m**2))
        eps = eps_frac * max(denom_scale, 1.0)  # scale-aware stabilizer
        return float(num / (denom_scale + eps))

    return 0.5 * (nrmse_at_q2(q2A) + nrmse_at_q2(q2B))

def per_function_chi2_from_rows(chi2_rows):
    """
    Build {"L": χ²_L, "T": χ²_T, "LT": χ²_LT, "TT": χ²_TT, "OTHER": χ²_other?}
    using the median χ² of the parameters belonging to each group.
    """
    def med(arr, idxs):
        if not idxs:
            return float("nan")
        return float(np.median(arr[idxs]))
    out = {}
    for g, idxs in PARAM_GROUPS.items():
        out[g] = med(chi2_rows, idxs)
    return out

def per_function_decision(group_name, pulls, I2_per, parsA, parsB, q2A, q2B, w, chi2_data=None):
    """
    Decide FE/RE/BAD for a parameter group with sensible gating:
      - Always compute robust dRMS first for physics groups.
      - Primary decision from (zmax, I2_mean, dRMS) with relaxed thresholds.
      - Data-χ² override can upgrade at most one level (BAD→RE or RE→FE),
        and only when tension/shape are not extreme.
    """

    # ---------- Tunables (use your TUNABLES dict if present) ----------
    ZS, IS = 1.8, 0.15        # small-tension override
    ZF, IF, DF = 2.5, 0.35, 0.50   # FE thresholds
    ZR, IR, DR = 4.0, 0.75, 0.90   # RE thresholds
    CHI2_FE_MAX = 1.25             # data-chi2 "good" threshold

    # "extreme tension" guardrails to block χ² override
    Z_EXTREME   = 6.0   # if zmax exceeds this, do not upgrade
    DRMS_EXTREME= 1.5   # if dRMS exceeds this, do not upgrade (clearly different shapes)

    # ---------- Pull/I^2 for this group ----------
    idxs = PARAM_GROUPS[group_name]
    z = pulls[idxs]
    I2 = I2_per[idxs]
    zmax = float(np.max(np.abs(z))) if len(z) else 0.0
    I2m  = float(np.mean(I2)) if len(I2) else 0.0

    # ---------- Compute dRMS BEFORE any decision (physics groups only) ----------
    dRMS = model_delta_rms(group_name, parsA, parsB, q2A, q2B, w)  # NaN for OTHER

    # ---------- 0) Small-tension override (protect tiny-signal cases) ----------
    # Only applies if both tension metrics are tiny; independent of dRMS.
    if (zmax <= ZS) and (I2m <= IS):
        return "FE", {"zmax": zmax, "I2_mean": I2m, "dRMS": dRMS}

    # ---------- 1) Primary decision from thresholds ----------
    # OTHER group: no sigma model → ignore dRMS (it's NaN), use z/I2 only
    if np.isnan(dRMS):
        if (zmax <= ZF) and (I2m < IF):
            primary = "FE"
        elif (zmax <= ZR) and (I2m < IR):
            primary = "RE"
        else:
            primary = "BAD"
    else:
        if (zmax <= ZF) and (I2m < IF) and (dRMS < DF):
            primary = "FE"
        elif (zmax <= ZR) and (I2m < IR) and (dRMS < DR):
            primary = "RE"
        else:
            primary = "BAD"

    # ---------- 2) Gated χ² override (at most one-level upgrade) ----------
    # Only consider if we actually have per-function chi2.
    if (chi2_data is not None) and (group_name in chi2_data):
        chi2v = chi2_data[group_name]
        if (chi2v == chi2v) and (chi2v <= CHI2_FE_MAX):
            # Allow upgrade only if not extreme tension/shape
            not_extreme = (zmax <= Z_EXTREME) and (np.isnan(dRMS) or dRMS <= DRMS_EXTREME)

            if not_extreme:
                if primary == "BAD":
                    primary = "RE"  # upgrade BAD→RE
                elif primary == "RE":
                    primary = "FE"  # upgrade RE→FE
                # if FE already, keep FE
            # else: keep primary (no upgrade in extreme cases)

    return primary, {"zmax": zmax, "I2_mean": I2m, "dRMS": dRMS}

# ---------- Main (adapted; preserves original outputs and adds per-function selection) ----------
def main(argv):
    # Defaults (preserved)
    default_files = [
        f"{TEMP_CACHEPATH}/pion/Q1p6W2p22/2025October02_H12M16S29/new_par.pl_Q16W222.dat",
        f"{TEMP_CACHEPATH}/pion/Q2p4W2p22/2025October02_H13M02S25/new_par.pl_Q24W222.dat"
    ]
    if len(argv) == 1:
        files = default_files
        print("No input files provided, using defaults:")
        for f in files: print("  ", f)
    elif len(argv) == 3:
        files = argv[1:3]
    else:
        print("Usage: python3 script.py <fileA> <fileB>")
        sys.exit(2)

    A_path = Path(files[0]); B_path = Path(files[1])
    q2A, wA = parse_Q2_W_from_path(A_path)
    q2B, wB = parse_Q2_W_from_path(B_path)
    q2_avg = 0.5*(q2A + q2B)
    w_avg  = 0.5*(wA  + wB )

    # NEW: Read both parameter files WITH per-row chi2
    valsA, errsA, chi2A, chi2A_rows = read_params_file_rows(A_path)
    valsB, errsB, chi2B, chi2B_rows = read_params_file_rows(B_path)

    # Compute per-function chi2 for each file (median within each group)
    chi2A_func = per_function_chi2_from_rows(chi2A_rows)
    chi2B_func = per_function_chi2_from_rows(chi2B_rows)

    # Combine per-function chi2 across datasets; conservative max(A,B)
    chi2_data = {}
    for g in PARAM_GROUPS.keys():
        a = chi2A_func.get(g, float("nan"))
        b = chi2B_func.get(g, float("nan"))
        chi2_data[g] = max(a, b) if (a == a and b == b) else (a if a == a else b)

    print("\nPer-function χ² (max of A,B medians over that group's parameter rows):")
    for g in sorted(chi2_data.keys()):
        print(f"  {g:5s}: {chi2_data[g]:.3f}")

    # ---- Original fixed-effects average (preserved) ----
    fe_vals, fe_errs = weighted_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B)
    write_params_file(AVG_OUTFILE, fe_vals, fe_errs, chi2_per_row=AVG_ROW_CHI2)
    print(f"\n[FE] Wrote original fixed-effects average to: {AVG_OUTFILE}")

    # ---- Random-effects global (for reference) + pulls ----
    re_vals, re_errs, pulls, tau2 = random_effects_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B)
    Q_per, Q_tot, df_tot, I2_per, I2_mean, pval = global_heterogeneity_stats(pulls)
    max_abs_z = float(np.max(np.abs(pulls))) if len(pulls) else 0.0
    print("\n=== Global diagnostics (all parameters) ===")
    if pval is None:
        print(f"Q_total = {Q_tot:.2f} (df={df_tot}), mean I^2 = {I2_mean:.3f}, max|z| = {max_abs_z:.2f}")
    else:
        print(f"Q_total = {Q_tot:.2f} (df={df_tot}), p = {pval:.3f}, mean I^2 = {I2_mean:.3f}, max|z| = {max_abs_z:.2f}")

    # ---- Per-function decisions (explicit) ----
    # We build a per-function "best" parameter vector, mixing FE/RE per group.
    best_vals = np.array(fe_vals, dtype=float)
    best_errs = np.array(fe_errs, dtype=float)
    perfunc_summary = []
    for func in ["L","T","LT","TT"] + (["OTHER"] if "OTHER" in PARAM_GROUPS else []):
        choice, mets = per_function_decision(func, pulls, I2_per,
                                            valsA, valsB, q2A, q2B, w_avg,
                                            chi2_data=chi2_data)
        idxs = PARAM_GROUPS[func]
        if choice == "FE":
            # keep FE (already in best_*)
            pass
        elif choice == "RE":
            best_vals[idxs] = re_vals[idxs]
            best_errs[idxs] = re_errs[idxs]
        elif choice == "BAD":
            # choose RE to be conservative, but flag as BAD (functional form likely needs change)
            best_vals[idxs] = re_vals[idxs]
            best_errs[idxs] = re_errs[idxs]
        perfunc_summary.append((func, choice, mets["zmax"], mets["I2_mean"], mets["dRMS"]))

    # ---- Save the selected per-function mixed parameters ----
    write_params_file(BEST_OUTFILE, best_vals, best_errs, chi2_per_row=AVG_ROW_CHI2)
    print(f"\n[BEST] Wrote per-function selected parameters to: {BEST_OUTFILE}")

    # ---- Print explicit per-function decisions ----
    print("\n=== Per-function decisions (explicit) ===")
    print(" func  | decision | max|z|  |  I2_mean  |  ΔRMS(A vs B) ")
    print("-------+----------+--------+-----------+---------------")
    for func, choice, zmax, I2m, dR in perfunc_summary:
        dR_str = f"{dR: .3f}" if (dR == dR) else "   NaN"  # NaN-safe formatting
        print(f" {func:5s} | {choice:8s} | {zmax:6.2f} | {I2m:9.3f} | {dR_str:>13s}")    
    # Guidance:
    if any(c=="BAD" for _,c,_,_,_ in perfunc_summary):
        print("\nNOTE: One or more functions are marked BAD — parameters disagree beyond tolerance and model shapes differ (ΔRMS high).")
        print("      This strongly suggests the **functional form** for those functions requires Q^2 evolution/adaptation.")

    # ---- Save diagnostics CSV ----
    try:
        with open(DIAG_CSV_OUTFILE, "w", newline="") as fcsv:
            w = csv.writer(fcsv)
            w.writerow(["Per-parameter diagnostics"])
            w.writerow(["par_index","pull_z","Q=z^2","tau2","I2"])
            for k in range(16):
                w.writerow([k+1, pulls[k], Q_per[k], tau2[k], I2_per[k]])
            w.writerow([])
            w.writerow(["Per-function summary"])
            w.writerow(["function","decision","max|z|","I2_mean","deltaRMS"])
            for row in perfunc_summary:
                w.writerow(list(row))
            w.writerow([])
            w.writerow(["Global"])
            w.writerow(["Q_total", Q_tot, "df_total", df_tot, "p_value", (float('nan') if pval is None else pval), "I2_mean", I2_mean, "max|z|", max_abs_z])
        print(f"Saved diagnostics to: {DIAG_CSV_OUTFILE}")
    except Exception as e:
        print("Warning: could not write diagnostics CSV:", e)

    # ---- Plot, now including the per-function 'Best' (mixed FE/RE) set ----
    datasets = [
        {"name": "Set A",           "pars": valsA,     "qq": q2A,     "ww": wA},
        {"name": "Set B",           "pars": valsB,     "qq": q2B,     "ww": wB},
        {"name": "Avg (weighted)",  "pars": fe_vals,   "qq": q2_avg,  "ww": w_avg},  # preserved
        {"name": "Best (per-func)", "pars": best_vals, "qq": q2_avg,  "ww": w_avg},
    ]
    make_four_panel(datasets)

if __name__ == "__main__":
    main(sys.argv)
