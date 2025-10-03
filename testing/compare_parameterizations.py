#! /usr/bin/python

#
# Description:
"""
Read two parameter files (new_par.pl_Q??W??.dat), parse par1..par16, parse Q^2 and W
from the pathname (e.g., Q1p6W2p22/...), then plot:
  1) 4-panel comparison (sigma_L, sigma_T, sigma_LT, sigma_TT)
  2) zoomed sigma_TT with the -t = m_pi^2 singularity marked

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
# Copyright (c) trottar
#

###############################################################################################################################################
# ltsep package import and pathing definitions

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
CACHEPATH=lt.CACHEPATH

TEMP_CACHEPATH=f"{OUTPATH}/testing_env"

###############################################################################################################################################

import sys
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---------- Constants ----------
mtar  = 0.93827231
mpipl = 0.139570

# ---------- Physics pieces (from your .model) ----------
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
                idx   = int(float(parts[2]))
            except ValueError:
                continue
            values[idx] = val
    return np.array([values.get(i, 0.0) for i in range(1, 17)], dtype=float)

# ---------- Plotting ----------
def make_four_panel(datasets, tmin=0.001, tmax=0.60, npts=1200, theta_cm=math.pi/2):
    obs = {
        "L":  (sigma_L,  r'$\sigma_L$'),
        "T":  (sigma_T,  r'$\sigma_T$'),
        "LT": (sigma_LT, r'$\sigma_{LT}$'),
        "TT": (sigma_TT, r'$\sigma_{TT}$'),
    }
    t_neg = -np.linspace(tmin, tmax, npts)
    tpos  = -t_neg

    fig = plt.figure(figsize=(10, 8))
    for i, (key, (func, ylabel)) in enumerate(obs.items(), start=1):
        ax = fig.add_subplot(2, 2, i)
        for data in datasets:
            y = func(data["pars"], data["qq"], t_neg, theta_cm, data["ww"])
            ax.plot(tpos, y, label=f'{data["name"]}  (Q^2={data["qq"]}, W={data["ww"]})')
        ax.set_xlabel(r'$-t\ \mathrm{[GeV^2]}$')
        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.grid(True)
        if i == 1:
            ax.legend()
    fig.suptitle("Comparison of $\\sigma_L,\\sigma_T,\\sigma_{LT},\\sigma_{TT}$")
    fig.tight_layout(rect=[0,0,1,0.96])
    plt.show()

def make_tt_zoom(datasets, tmin=0.001, tmax=0.60, npts=1200, theta_cm=math.pi/2, exclude=0.004):
    t_neg = -np.linspace(tmin, tmax, npts)
    tpos  = -t_neg
    t_sing = mpipl**2

    all_vals = []
    for data in datasets:
        y = sigma_TT(data["pars"], data["qq"], t_neg, theta_cm, data["ww"])
        mask = (np.abs(tpos - t_sing) > exclude)
        all_vals.append(y[mask])
    finite_vals = np.concatenate(all_vals)
    ymax = np.nanpercentile(np.abs(finite_vals), 99.5)

    plt.figure(figsize=(7,5))
    for data in datasets:
        y = sigma_TT(data["pars"], data["qq"], t_neg, theta_cm, data["ww"])
        plt.plot(tpos, y, label=f'{data["name"]}  (Q^2={data["qq"]}, W={data["ww"]})')
    plt.axvline(t_sing, linestyle=':', linewidth=2, label=r'Singularity at $-t=m_\pi^2$')
    plt.ylim(-ymax, ymax)
    plt.xlim(tmin, tmax)
    plt.xlabel(r'$-t\ \mathrm{[GeV^2]}$')
    plt.ylabel(r'$\sigma_{TT}$ (arb.)')
    plt.title(r'$\sigma_{TT}$ zoomed (singularity marked)$\,$')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

# ---------- Main ----------
def main(argv):
    default_files = [
        "{TEMP_CACHEPATH}/pion/Q1p6W2p22/2025October02_H12M16S29/new_par.pl_Q16W222.dat",
        "{TEMP_CACHEPATH}/pion/Q2p4W2p22/2025October02_H13M02S25/new_par.pl_Q24W222.dat"
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
        print("If no arguments given, uses default files.")
        sys.exit(2)

    datasets = []
    for i, f in enumerate(files):
        path = Path(f)
        q2, w = parse_Q2_W_from_path(path)
        pars = read_pars_file(path)
        datasets.append({
            "name": f"Set {chr(ord('A')+i)}",
            "pars": pars,
            "qq": q2,
            "ww": w
        })

    make_four_panel(datasets)
    make_tt_zoom(datasets)

if __name__ == "__main__":
    main(sys.argv)
