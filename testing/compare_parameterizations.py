#! /usr/bin/python

#
# Description:
"""
Read two parameter files (new_par.pl_Q??W??.dat), parse par1..par16 + errors + chi2,
parse Q^2 and W from the path, then plot a 4-panel comparison:
  - Set A
  - Set B
  - Weighted average of A and B (using errors + chi2):
      w_k^(f) = (1/chi2_f) * (1/err_{k,f}^2)
      pbar_k  = sum_f w_k^(f) p_{k,f} / sum_f w_k^(f)
      sbar_k  = 1 / sqrt(sum_f w_k^(f))
Also writes the averaged parameters to:
  new_par.pl_Q1p6-2p4W2p22.dat
with the same spacing/columns as the inputs and a default chi2 of 3.0 per row.

Defaults (used if no args):
  <OUTPATH>/testing_env/pion/Q1p6W2p22/2025October02_H12M16S29/new_par.pl_Q16W222.dat
  <OUTPATH>/testing_env/pion/Q2p4W2p22/2025October02_H13M02S25/new_par.pl_Q24W222.dat
"""
# ================================================================
# Time-stamp: "2025-10-01 08:44:00 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
###############################################################################################################################################

import sys,os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

###############################################################################################################################################
# ltsep package import and pathing definitions

from ltsep import Root
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

USER=lt.USER
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
CACHEPATH=lt.CACHEPATH

TEMP_CACHEPATH=f"{OUTPATH}/testing_env"

###############################################################################################################################################
# t-range
TMIN = 0.02
TMAX = 0.60

file_str = "pl_Q1p6-2p4W2p22"

# Output filename for averaged parameters
AVG_OUTFILE = f"{LTANAPATH}/testing/parameters/new_par.{file_str}.dat"
AVG_ROW_CHI2 = 3.0  # default chi2 for every parameter row in the saved file

###############################################################################################################################################
# ---------- Constants ----------
mtar  = 0.93827231
mpipl = 0.139570

# ---------- Physics pieces ----------
def wfactor(ww):
    return 1.0 / ((ww**2) - (mtar**2))**2.0

def sigma_L(pars, qq, t, theta_cm, ww):
    par1, par2, par3, par4 = pars[0], pars[1], pars[2], pars[3]
    return (wfactor(ww)
            * (par1 * qq / (1 + par2*qq + 0.05*(qq**2))**2)
            * np.exp((par3 - par4*np.log(qq)) * -np.abs(t)))

def sigma_T(pars, qq, t, theta_cm, ww):
    par5, par6 = pars[4], pars[5]
    return wfactor(ww) * ((par5 / qq) + (par6 / (qq**2))) * np.ones_like(t)

def sigma_LT(pars, qq, t, theta_cm, ww):
    par9, par10, par11, par12 = pars[8], pars[9], pars[10], pars[11]
    return (wfactor(ww)
            * (np.exp(par9 + (par10 / np.sqrt(qq)) * -np.abs(t)) + par11 - (par12 / (qq**2)))
            * np.sin(theta_cm))

def sigma_TT(pars, qq, t, theta_cm, ww):
    par13 = float(pars[12])
    abs_t = np.abs(t)
    denom = (-abs_t + mpipl**2)**2
    denom = np.clip(denom, 1e-18, None)
    return (wfactor(ww)
            * (par13 / (qq**2))
            * (-abs_t / denom)
            * (np.sin(theta_cm)**2))

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

def read_pars_file(par_path: Path):
    values = {}
    errors = {}
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

def weighted_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B):
    wA = (1.0/chi2A) * (1.0/np.square(errsA))
    wB = (1.0/chi2B) * (1.0/np.square(errsB))
    wsum = np.clip(wA + wB, 1e-18, None)
    pbar = (wA*valsA + wB*valsB) / wsum
    sbar = 1.0 / np.sqrt(wsum)
    return pbar, sbar

def write_params_file(filename, values, errors, chi2_per_row=3.0):
    """
    Write lines with the same spacing/columns as your inputs:
      <value>   <error>   <index>          <chi2>
    Using scientific notation and aligned columns.
    """
    with open(filename, "w") as f:
        for idx in range(1, 17):
            v = values[idx-1]
            e = errors[idx-1]
            line = f"{v: .8e}   {e: .8e}   {idx:2d}          {chi2_per_row:.1f}\n"
            f.write(line)

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

# ---------- Main ----------
def main(argv):
    default_files = [
        f"{TEMP_CACHEPATH}/pion/Q1p6W2p22/2025October02_H12M16S29/new_par.pl_Q16W222.dat",
        f"{TEMP_CACHEPATH}/pion/Q2p4W2p22/2025October02_H13M02S25/new_par.pl_Q24W222.dat"
    ]

    if len(argv) == 1:
        files = default_files
        print("No input files provided, using defaults:")
        for f in files:
            print("  ", f)
    elif len(argv) == 3:
        files = argv[1:3]
    else:
        print("Usage: python3 script.py <fileA> <fileB>")
        sys.exit(2)

    A_path = Path(files[0]); B_path = Path(files[1])
    q2A, wA = parse_Q2_W_from_path(A_path)
    q2B, wB = parse_Q2_W_from_path(B_path)

    valsA, errsA, chi2A = read_pars_file(A_path)
    valsB, errsB, chi2B = read_pars_file(B_path)

    avg_vals, avg_errs = weighted_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B)
    q2_avg = 0.5*(q2A + q2B)
    w_avg  = 0.5*(wA  + wB )

    print("\nWeighted-average parameters (w = (1/chi2)*(1/err^2)):")
    for i, (v, s) in enumerate(zip(avg_vals, avg_errs), start=1):
        print(f"  par{i:02d} = {v:+.8e}   +/- {s:.3e}")

    write_params_file(AVG_OUTFILE, avg_vals, avg_errs, chi2_per_row=AVG_ROW_CHI2)
    print(f"\nWrote averaged parameters to: {AVG_OUTFILE}")
    print(f"  (Q^2_avg={q2_avg:.3g}, W_avg={w_avg:.3g}, chi2 per row written = {AVG_ROW_CHI2:.1f})")

    datasets = [
        {"name": "Set A",           "pars": valsA,    "qq": q2A,    "ww": wA},
        {"name": "Set B",           "pars": valsB,    "qq": q2B,    "ww": wB},
        {"name": "Avg (weighted)",  "pars": avg_vals, "qq": q2_avg, "ww": w_avg},
    ]

    make_four_panel(datasets)

if __name__ == "__main__":
    main(sys.argv)

# ========================= ADDED UTILITIES =========================
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

def combine_params(valsA, errsA, chi2A, valsB, errsB, chi2B, mode="auto", pull_thresh=2.0, verbose=True):
    fe_vals, fe_errs = weighted_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B)
    if mode == "fe":
        return fe_vals, fe_errs, None
    re_vals, re_errs, pulls, tau2 = random_effects_average_two_fits(valsA, errsA, chi2A, valsB, errsB, chi2B)
    if mode == "re":
        return re_vals, re_errs, {"pulls":pulls, "tau2":tau2}
    use_re = np.abs(pulls) > float(pull_thresh)
    avg_vals = np.where(use_re, re_vals, fe_vals)
    avg_errs = np.where(use_re, re_errs, fe_errs)
    return avg_vals, avg_errs, {"pulls":pulls, "tau2":tau2, "use_re":use_re}
# ========================= END ADDED UTILITIES =========================
