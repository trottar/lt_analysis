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
  - Best (selected) — per-function selection (FE vs RE) with diagnostics

This version FIXES all calls to the dynamically generated functions produced by
`prepare_equations(...)`. Every call now uses a consistent dispatcher that matches
the signature:

    sig_*_optimized(q2_set, w_set, qq, ww, tt, theta_cm, parA, parB, parC, parD)

We also use per-row χ² (last column) to compute per-function χ² and include
robust ΔRMS shape diagnostics.
"""

# ================================================================
# Time-stamp: "2025-10-03"
# ================================================================

import sys, os, re, math, csv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------- ltsep paths ----------
from ltsep import Root
lt = Root(os.path.realpath(__file__), "Plot_LTSep")
OUTPATH   = lt.OUTPATH
LTANAPATH = lt.LTANAPATH
TEMP_CACHEPATH = f"{OUTPATH}/testing_env"

# ---------- import utility ----------
sys.path.append("utility")
from utility import load_equations, prepare_equations

# ---------- t-range ----------
TMIN = 0.025
TMAX = 0.600

# Output filenames
file_str           = "pl_Q1p6-2p4W2p22_trotta"
AVG_OUTFILE        = f"{LTANAPATH}/testing/parameters/new_par.{file_str}.dat"               # original fixed-effects output (preserved)
BEST_OUTFILE       = f"{LTANAPATH}/testing/parameters/new_par.{file_str}_best.dat"          # per-function FE/RE mix
DIAG_CSV_OUTFILE   = f"{LTANAPATH}/testing/parameters/new_par.{file_str}_diagnostics.csv"   # diagnostics
AVG_ROW_CHI2       = 3.0  # default chi2 for every parameter row in saved files

# Load equations from model
equations = load_equations("Q1p6W2p22.model")

# ========================= Tunables =========================
TUNABLES = {
    # Small-tension override → FE
    "Z_SMALL":   1.8,
    "I2_SMALL":  0.15,

    # FE region (GOOD)
    "Z_FE":      2.5,
    "I2_FE":     0.35,
    "DRMS_FE":   0.50,

    # RE region (OK)
    "Z_RE":      4.0,
    "I2_RE":     0.75,
    "DRMS_RE":   0.90,

    # Robust ΔRMS settings
    "TAU_MASK":  1e-3,   # mask t-points where both predictions are tiny
    "EPS_FRAC":  1e-6,   # additive stabilizer fraction of signal scale

    # Data-χ² (from rows) override: allow at most one-level upgrade if χ² is good
    "CHI2_DATA_FE_MAX": 1.25,
}

# ======================= Function → parameter index mapping (0-based) =======================
PARAM_GROUPS = {
    "L":  [0,1,2,3],        # par1..par4
    "T":  [4,5,6,7],        # par5..par8
    "LT": [8,9,10,11],      # par9..par12
    "TT": [12,13,14,15],    # par13..par16
}
ALL_USED = sorted(set(sum(PARAM_GROUPS.values(), [])))
OTHER_IDX = [i for i in range(16) if i not in ALL_USED]
if OTHER_IDX:
    PARAM_GROUPS["OTHER"] = OTHER_IDX

# ---------- Constants ----------
mtar  = 0.93827231
mpipl = 0.139570

# ---------- Build dynamic optimized functions ----------
wfactor_opt = prepare_equations(equations, 'wfactor')
sigL_opt    = prepare_equations(equations, 'sig_L')
sigT_opt    = prepare_equations(equations, 'sig_T')
sigLT_opt   = prepare_equations(equations, 'sig_LT')
sigTT_opt   = prepare_equations(equations, 'sig_TT')

OPT_FUNCS = {
    "L":  sigL_opt,
    "T":  sigT_opt,
    "LT": sigLT_opt,
    "TT": sigTT_opt,
}

# ---------- Consistent dispatcher (MATCHES prepare_equations signature) ----------
def _slice_pars(group_name, pars):
    """Return the 4 parameters (in order) for a given group from the 16-vector."""
    idxs = PARAM_GROUPS[group_name]
    return (pars[idxs[0]], pars[idxs[1]], pars[idxs[2]], pars[idxs[3]])

def eval_sigma(group_name, pars, q2_set, w_set, qq, ww, tt, theta_cm):
    """
    Always call the prepared functions via this dispatcher:

        sig_*_optimized(q2_set, w_set, qq, ww, tt, theta_cm, p1, p2, p3, p4)

    Returns a 1D numpy array with the SAME length as `tt` (broadcasting scalars if needed).
    """
    if group_name not in OPT_FUNCS:
        raise KeyError(f"Missing optimized function for {group_name}")
    f = OPT_FUNCS[group_name]
    p1, p2, p3, p4 = _slice_pars(group_name, pars)

    # Call the optimized function
    y = f(q2_set, w_set, qq, ww, tt, theta_cm, p1, p2, p3, p4)

    # --- Shape normalization (critical for plotting) ---
    tt_arr = np.asarray(tt)
    y_arr  = np.asarray(y)

    # If scalar or 0-d, broadcast to tt shape
    if y_arr.ndim == 0:
        y_arr = np.full_like(tt_arr, float(y_arr), dtype=float)

    # If (1,) or otherwise not matching tt length, try broadcasting
    elif y_arr.shape[0] != tt_arr.shape[0]:
        try:
            y_arr = np.broadcast_to(y_arr, tt_arr.shape)
        except ValueError:
            # Fallback: force to same length with uniform fill of first element
            y_arr = np.full_like(tt_arr, float(np.ravel(y_arr)[0]), dtype=float)

    # Ensure 1D vector
    return np.ravel(y_arr)

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

# ---------- Plotting (uses eval_sigma) ----------
def make_four_panel(datasets, npts=1200, theta_cm=math.pi/2):
    obs = ["L", "T", "LT", "TT"]
    t_neg = -np.linspace(TMIN, TMAX, npts)
    tpos  = -t_neg

    fig = plt.figure(figsize=(10, 8))
    for i, key in enumerate(obs, start=1):
        ax = fig.add_subplot(2, 2, i)
        for data in datasets:
            y = eval_sigma(
                key,
                pars=data["pars"],
                q2_set=data["qq"], w_set=data["ww"],
                qq=data["qq"],    ww=data["ww"],
                tt=t_neg, theta_cm=theta_cm
            )
            label = f'{data["name"]}  (Q^2={data["qq"]:.3g}, W={data["ww"]:.3g})'
            ax.plot(tpos, y, label=label)
        ax.set_xlabel(r'$-t\ \mathrm{[GeV^2]}$')
        ax.set_ylabel(rf'$\sigma_{{{key}}}$')
        ax.set_title(rf'$\sigma_{{{key}}}$')
        ax.grid(True)
        if i == 1:
            ax.legend()
    fig.suptitle("Comparison of $\\sigma_L,\\sigma_T,\\sigma_{LT},\\sigma_{TT}$")
    fig.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(f"{TEMP_CACHEPATH}/fits_{file_str}.pdf")
    plt.close(fig)

# ========================= Diagnostics + robust combination =========================
def param_pulls(valsA, errsA, valsB, errsB):
    denom = np.sqrt(np.square(errsA) + np.square(errsB))
    denom = np.clip(denom, 1e-18, None)
    return (valsA - valsB) / denom

def dersimonian_laird_tau2(vA, sA, vB, sB, chi2A=None, chi2B=None):
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

# ---------- Robust ΔRMS using the dispatcher ----------
def model_delta_rms(group_name, parsA, parsB, q2A, q2B, w, npts=400, theta_cm=math.pi/2):
    """
    Robust, scale-aware shape disagreement for a given sigma group (L/T/LT/TT):
      dRMS = 0.5 * [ NRMSE(q2A) + NRMSE(q2B) ]
    Using the SAME calling convention as plotting via eval_sigma(...).
    """
    if group_name not in {"L","T","LT","TT"}:
        return float("nan")

    t = -np.linspace(TMIN, TMAX, npts)
    tau     = float(TUNABLES["TAU_MASK"])
    eps_frac= float(TUNABLES["EPS_FRAC"])

    def nrmse_at(q2_set):
        A = eval_sigma(group_name, parsA, q2_set, w, q2_set, w, t, theta_cm)
        B = eval_sigma(group_name, parsB, q2_set, w, q2_set, w, t, theta_cm)
        avg_mag = 0.5*(np.abs(A)+np.abs(B))
        ref = np.percentile(avg_mag, 75)
        mask = avg_mag > (tau * max(ref, 1e-18))
        if not np.any(mask):
            return 0.0
        A_m, B_m   = A[mask], B[mask]
        avg_mag_m  = avg_mag[mask]
        num  = np.sqrt(np.mean((A_m - B_m)**2))
        denS = np.sqrt(np.mean(avg_mag_m**2))
        eps  = eps_frac * max(denS, 1.0)
        return float(num / (denS + eps))

    return 0.5*(nrmse_at(q2A) + nrmse_at(q2B))

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
      - Compute robust dRMS first (for L/T/LT/TT).
      - Primary decision from (zmax, I2_mean, dRMS) with relaxed thresholds.
      - Data-χ² override can upgrade at most one level (BAD→RE or RE→FE),
        and only when tension/shape are not extreme.
    """
    ZS, IS = TUNABLES["Z_SMALL"], TUNABLES["I2_SMALL"]
    ZF, IF, DF = TUNABLES["Z_FE"], TUNABLES["I2_FE"], TUNABLES["DRMS_FE"]
    ZR, IR, DR = TUNABLES["Z_RE"], TUNABLES["I2_RE"], TUNABLES["DRMS_RE"]
    CHI2_FE_MAX = TUNABLES["CHI2_DATA_FE_MAX"]

    Z_EXTREME    = 6.0
    DRMS_EXTREME = 1.5

    idxs = PARAM_GROUPS[group_name]
    z = pulls[idxs]
    I2 = I2_per[idxs]
    zmax = float(np.max(np.abs(z))) if len(z) else 0.0
    I2m  = float(np.mean(I2)) if len(I2) else 0.0

    dRMS = model_delta_rms(group_name, parsA, parsB, q2A, q2B, w)  # NaN for OTHER

    # 0) Small-tension override
    if (zmax <= ZS) and (I2m <= IS):
        return "FE", {"zmax": zmax, "I2_mean": I2m, "dRMS": dRMS}

    # 1) Primary decision
    if np.isnan(dRMS):  # OTHER group
        if (zmax <= ZF) and (I2m < IF): primary = "FE"
        elif (zmax <= ZR) and (I2m < IR): primary = "RE"
        else: primary = "BAD"
    else:
        if (zmax <= ZF) and (I2m < IF) and (dRMS < DF): primary = "FE"
        elif (zmax <= ZR) and (I2m < IR) and (dRMS < DR): primary = "RE"
        else: primary = "BAD"

    # 2) Gated χ² override (upgrade at most one level)
    if (chi2_data is not None) and (group_name in chi2_data):
        chi2v = chi2_data[group_name]
        if (chi2v == chi2v) and (chi2v <= CHI2_FE_MAX):
            not_extreme = (zmax <= Z_EXTREME) and (np.isnan(dRMS) or dRMS <= DRMS_EXTREME)
            if not_extreme:
                if primary == "BAD": primary = "RE"
                elif primary == "RE": primary = "FE"

    return primary, {"zmax": zmax, "I2_mean": I2m, "dRMS": dRMS}

# ---------- Main ----------
def main(argv):
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

    # Read both parameter files WITH per-row chi2
    valsA, errsA, chi2A, chi2A_rows = read_params_file_rows(A_path)
    valsB, errsB, chi2B, chi2B_rows = read_params_file_rows(B_path)

    # Per-function χ² from rows (median over that group's indices)
    def per_function_chi2_from_rows_local(chi2_rows):
        out = {}
        for g, idxs in PARAM_GROUPS.items():
            if not idxs: out[g] = float("nan")
            else: out[g] = float(np.median(chi2_rows[idxs]))
        return out

    chi2A_func = per_function_chi2_from_rows_local(chi2A_rows)
    chi2B_func = per_function_chi2_from_rows_local(chi2B_rows)

    # Combine per-function χ² conservatively (max)
    chi2_data = {}
    for g in PARAM_GROUPS.keys():
        a, b = chi2A_func.get(g, float("nan")), chi2B_func.get(g, float("nan"))
        chi2_data[g] = max(a, b) if (a == a and b == b) else (a if a == a else b)

    print("\nPer-function χ² (max of A,B medians over that group's rows):")
    for g in sorted(chi2_data.keys()):
        print(f"  {g:5s}: {chi2_data[g]:.3f}")

    # Fixed-effects average (preserved)
    fe_vals, fe_errs = weighted_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B)
    write_params_file(AVG_OUTFILE, fe_vals, fe_errs, chi2_per_row=AVG_ROW_CHI2)
    print(f"\n[FE] Wrote original fixed-effects average to: {AVG_OUTFILE}")

    # Random-effects + pulls
    re_vals, re_errs, pulls, tau2 = random_effects_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B)
    Q_per, Q_tot, df_tot, I2_per, I2_mean, pval = global_heterogeneity_stats(pulls)
    max_abs_z = float(np.max(np.abs(pulls))) if len(pulls) else 0.0
    print("\n=== Global diagnostics (all parameters) ===")
    if pval is None:
        print(f"Q_total = {Q_tot:.2f} (df={df_tot}), mean I^2 = {I2_mean:.3f}, max|z| = {max_abs_z:.2f}")
    else:
        print(f"Q_total = {Q_tot:.2f} (df={df_tot}), p = {pval:.3f}, mean I^2 = {I2_mean:.3f}, max|z| = {max_abs_z:.2f}")

    # Per-function decisions
    best_vals = np.array(fe_vals, dtype=float)
    best_errs = np.array(fe_errs, dtype=float)
    perfunc_summary = []
    order = ["L","T","LT","TT"] + (["OTHER"] if "OTHER" in PARAM_GROUPS else [])
    for func in order:
        choice, mets = per_function_decision(func, pulls, I2_per,
                                             valsA, valsB, q2A, q2B, w_avg,
                                             chi2_data=chi2_data)
        idxs = PARAM_GROUPS[func]
        if choice == "FE":
            pass
        elif choice in ("RE","BAD"):
            best_vals[idxs] = re_vals[idxs]
            best_errs[idxs] = re_errs[idxs]
        perfunc_summary.append((func, choice, mets["zmax"], mets["I2_mean"], mets["dRMS"]))

    write_params_file(BEST_OUTFILE, best_vals, best_errs, chi2_per_row=AVG_ROW_CHI2)
    print(f"\n[BEST] Wrote per-function selected parameters to: {BEST_OUTFILE}")

    # Print decisions
    print("\n=== Per-function decisions (explicit) ===")
    print(" func  | decision | max|z|  |  I2_mean  |  ΔRMS(A vs B) ")
    print("-------+----------+--------+-----------+---------------")
    for func, choice, zmax, I2m, dR in perfunc_summary:
        dR_str = f"{dR: .3f}" if (dR == dR) else "   NaN"
        print(f" {func:5s} | {choice:8s} | {zmax:6.2f} | {I2m:9.3f} | {dR_str:>13s}")

    if any(c=="BAD" for _,c,_,_,_ in perfunc_summary):
        print("\nNOTE: One or more functions are marked BAD — parameters disagree beyond tolerance and/or shapes differ.")
        print("      This suggests the functional form for those functions may need Q^2 evolution/adaptation.")

    # Save diagnostics CSV
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

    # Plot (now uses eval_sigma and thus the correct optimized signature)
    datasets = [
        {"name": "Set A",           "pars": valsA,     "qq": q2A,     "ww": wA},
        {"name": "Set B",           "pars": valsB,     "qq": q2B,     "ww": wB},
        {"name": "Avg (weighted)",  "pars": fe_vals,   "qq": q2_avg,  "ww": w_avg},
        {"name": "Best (per-func)", "pars": best_vals, "qq": q2_avg,  "ww": w_avg},
    ]
    make_four_panel(datasets)

if __name__ == "__main__":
    main(sys.argv)
