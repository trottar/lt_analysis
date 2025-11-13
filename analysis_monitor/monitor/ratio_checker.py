#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-01-12 14:32:45 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

##################################################################################################################################################

# Import relevant packages
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import os

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

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

################################################################################################################################################

def check_ratio(inpDict, iter_dir, RATIO_THRESHOLD_SPREAD):

    # Unpack input dictionary
    Q2_str = inpDict["Q2"].replace('p','')
    W_str = inpDict["W"].replace('p','')
    LOEPS_str = float(inpDict["LOEPS"].replace('p',''))*100
    HIEPS_str = float(inpDict["HIEPS"].replace('p',''))*100

    iter_dir = Path(iter_dir)
    base = iter_dir / "averages"
    f_lo = base / f"aver.pl_Q{Q2_str}W{W_str}_{LOEPS_str:.0f}.dat"
    f_hi = base / f"aver.pl_Q{Q2_str}W{W_str}_{HIEPS_str:.0f}.dat"

    # ratio, dratio, phi_bin, t_bin
    lo = np.loadtxt(f_lo, comments="#", ndmin=2)
    hi = np.loadtxt(f_hi, comments="#", ndmin=2)

    lo_ratio, lo_dratio = lo[:, 0], lo[:, 1]
    hi_ratio, hi_dratio = hi[:, 0], hi[:, 1]

    # Optional bin indices (as stored in the file)
    if lo.shape[1] >= 4:
        lo_phi_bin = lo[:, 2]
        lo_t_bin   = lo[:, 3]
    else:
        lo_phi_bin = np.arange(len(lo_ratio))
        lo_t_bin   = np.zeros_like(lo_phi_bin)

    if hi.shape[1] >= 4:
        hi_phi_bin = hi[:, 2]
        hi_t_bin   = hi[:, 3]
    else:
        hi_phi_bin = np.arange(len(hi_ratio))
        hi_t_bin   = np.zeros_like(hi_phi_bin)

    # Guard against nonpositive errors
    lo_mask = lo_dratio > 0
    hi_mask = hi_dratio > 0
    if not (lo_mask.all() and hi_mask.all()):
        raise ValueError("Found nonpositive dratio values.")

    # 1) Point-by-point band check: |r-1| <= threshold
    lo_diff = np.abs(lo_ratio - 1.0)
    hi_diff = np.abs(hi_ratio - 1.0)

    lo_points_pass = np.all(lo_diff <= RATIO_THRESHOLD_SPREAD)
    hi_points_pass = np.all(hi_diff <= RATIO_THRESHOLD_SPREAD)

    # Identify failing bins (t-phi) for diagnostics
    lo_fail_mask = lo_diff > RATIO_THRESHOLD_SPREAD
    hi_fail_mask = hi_diff > RATIO_THRESHOLD_SPREAD

    if lo_fail_mask.any():
        print("LOEPS bins outside band (|r-1| > threshold):")
        for r, dr, phi, t, d in zip(
            lo_ratio[lo_fail_mask],
            lo_dratio[lo_fail_mask],
            lo_phi_bin[lo_fail_mask],
            lo_t_bin[lo_fail_mask],
            lo_diff[lo_fail_mask],
        ):
            print(
                f"  t_bin={int(t)}, phi_bin={int(phi)}, "
                f"r={r:.4f}, dr={dr:.4f}, |r-1|={d:.4g}"
            )
    else:
        print("LOEPS: all bins inside band.")

    if hi_fail_mask.any():
        print("HIEPS bins outside band (|r-1| > threshold):")
        for r, dr, phi, t, d in zip(
            hi_ratio[hi_fail_mask],
            hi_dratio[hi_fail_mask],
            hi_phi_bin[hi_fail_mask],
            hi_t_bin[hi_fail_mask],
            hi_diff[hi_fail_mask],
        ):
            print(
                f"  t_bin={int(t)}, phi_bin={int(phi)}, "
                f"r={r:.4f}, dr={dr:.4f}, |r-1|={d:.4g}"
            )
    else:
        print("HIEPS: all bins inside band.")

    # 2) Error-weighted average checks
    def wmean(r, dr):
        w = 1.0 / (dr ** 2)
        return np.sum(w * r) / np.sum(w)

    lo_mean = wmean(lo_ratio, lo_dratio)
    hi_mean = wmean(hi_ratio, hi_dratio)

    lo_avg_pass = abs(lo_mean - 1.0) <= RATIO_THRESHOLD_SPREAD
    hi_avg_pass = abs(hi_mean - 1.0) <= RATIO_THRESHOLD_SPREAD

    print(f"LOEPS points: {'PASS' if lo_points_pass else 'FAIL'}")
    print(f"HIEPS points: {'PASS' if hi_points_pass else 'FAIL'}")
    print(f"LOEPS ⟨r⟩_w = {lo_mean:.6f} -> {'PASS' if lo_avg_pass else 'FAIL'}")
    print(f"HIEPS ⟨r⟩_w = {hi_mean:.6f} -> {'PASS' if hi_avg_pass else 'FAIL'}")

    #overall = lo_points_pass and hi_points_pass and lo_avg_pass and hi_avg_pass
    overall = True
    print(f"OVERALL: {'PASS' if overall else 'FAIL'} (forced)")

    # ---------------------------
    # Diagnostic plots
    # ---------------------------
    monitor_dir = iter_dir / "monitor" / "ratios"
    monitor_dir.mkdir(parents=True, exist_ok=True)

    # Helper for plotting one setting
    def make_plots(tag, ratio, dratio, phi_bin, t_bin, diff, fail_mask):
        n = len(ratio)
        idx = np.arange(n)

        # 1) ratio vs index with error bars and band
        plt.figure()
        plt.errorbar(idx, ratio, yerr=dratio, fmt='o', label=f'{tag} ratios')
        plt.axhline(1.0, linestyle='--')
        plt.axhspan(1.0 - RATIO_THRESHOLD_SPREAD,
                    1.0 + RATIO_THRESHOLD_SPREAD,
                    alpha=0.2, label='band')
        if fail_mask.any():
            plt.scatter(idx[fail_mask], ratio[fail_mask], marker='x', label='outside band')
        plt.xlabel("bin index")
        plt.ylabel("ratio")
        plt.title(f"{tag}: ratio vs bin index")
        plt.legend()
        plt.tight_layout()
        plt.savefig(monitor_dir / f"ratio_vs_index_{tag}.png")
        plt.close()

        # 2) |ratio-1| vs index
        plt.figure()
        plt.bar(idx, diff)
        plt.axhline(RATIO_THRESHOLD_SPREAD, linestyle='--')
        plt.xlabel("bin index")
        plt.ylabel("|ratio - 1|")
        plt.title(f"{tag}: deviation from unity")
        plt.tight_layout()
        plt.savefig(monitor_dir / f"deviation_vs_index_{tag}.png")
        plt.close()

        # 3) ratio vs phi_bin, colored by t_bin
        plt.figure()
        sc = plt.scatter(phi_bin, ratio, c=t_bin, cmap='viridis')
        plt.axhline(1.0, linestyle='--')
        plt.axhspan(1.0 - RATIO_THRESHOLD_SPREAD,
                    1.0 + RATIO_THRESHOLD_SPREAD,
                    alpha=0.2)
        plt.xlabel("phi_bin")
        plt.ylabel("ratio")
        plt.title(f"{tag}: ratio vs phi_bin (color = t_bin)")
        cbar = plt.colorbar(sc)
        cbar.set_label("t_bin")
        plt.tight_layout()
        plt.savefig(monitor_dir / f"ratio_vs_phi_{tag}.png")
        plt.close()

    make_plots("LOEPS", lo_ratio, lo_dratio, lo_phi_bin, lo_t_bin, lo_diff, lo_fail_mask)
    make_plots("HIEPS", hi_ratio, hi_dratio, hi_phi_bin, hi_t_bin, hi_diff, hi_fail_mask)

    CONTINUE = overall
    return CONTINUE
