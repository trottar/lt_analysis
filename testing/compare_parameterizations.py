#! /usr/bin/python

#
# Description:
"""
Read two parameter files (new_par.pl_Q??W??.dat), parse par1..par16 + errors + chi2,
parse Q^2 and W from the path, then plot:
  4-panel comparison (sigma_L, sigma_T, sigma_LT, sigma_TT) for:
    - Set A
    - Set B
    - Weighted average of A and B (using errors + chi2)
Default files (used if no arguments are provided):
  Q1p6W2p22/2025October02_H12M16S29/new_par.pl_Q16W222.dat
  Q2p4W2p22/2025October02_H13M02S25/new_par.pl_Q24W222.dat
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
TMIN = 0.020
TMAX = 0.600

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
    plt.show()

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

    datasets = [
        {"name": "Set A",           "pars": valsA,    "qq": q2A,    "ww": wA},
        {"name": "Set B",           "pars": valsB,    "qq": q2B,    "ww": wB},
        {"name": "Avg (weighted)",  "pars": avg_vals, "qq": q2_avg, "ww": w_avg},
    ]

    make_four_panel(datasets)

if __name__ == "__main__":
    main(sys.argv)
