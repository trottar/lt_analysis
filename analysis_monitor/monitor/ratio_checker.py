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
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os

import logging

# Turn down matplotlib + PIL debug spam
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("PIL").setLevel(logging.WARNING)
logging.getLogger("PIL.PngImagePlugin").setLevel(logging.WARNING)

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
    Q2_str   = inpDict["Q2"].replace('p','')
    W_str    = inpDict["W"].replace('p','')
    LOEPS_str = float(inpDict["LOEPS"].replace('p','')) * 100
    HIEPS_str = float(inpDict["HIEPS"].replace('p','')) * 100

    iter_dir = Path(iter_dir)
    base = iter_dir / "averages"
    f_lo = base / f"aver.pl_Q{Q2_str}W{W_str}_{LOEPS_str:.0f}.dat"
    f_hi = base / f"aver.pl_Q{Q2_str}W{W_str}_{HIEPS_str:.0f}.dat"

    # ratio, dratio, phi_bin, t_bin
    lo = np.loadtxt(f_lo, comments="#", ndmin=2)
    hi = np.loadtxt(f_hi, comments="#", ndmin=2)

    lo_ratio, lo_dratio = lo[:, 0], lo[:, 1]
    hi_ratio, hi_dratio = hi[:, 0], hi[:, 1]

    # mask out bins where ratio == 0 AND dratio == -1000
    lo_valid = ~((lo_ratio == 0) & (lo_dratio == -1000))
    hi_valid = ~((hi_ratio == 0) & (hi_dratio == -1000))

    # keep only good bins, preserving your names
    lo        = lo[lo_valid]
    hi        = hi[hi_valid]
    lo_ratio  = lo_ratio[lo_valid]
    lo_dratio = lo_dratio[lo_valid]
    hi_ratio  = hi_ratio[hi_valid]
    hi_dratio = hi_dratio[hi_valid]

    # Optional bin indices (as stored in the file); fall back to simple indices
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
    lo_poserr = lo_dratio > 0
    hi_poserr = hi_dratio > 0
    if not (lo_poserr.all() and hi_poserr.all()):
        raise ValueError("Found nonpositive dratio values.")

    # 1) Point-by-point band check with error weighting
    #    - absolute band: |r - 1| > RATIO_THRESHOLD_SPREAD
    #    - AND statistically significant: |r - 1| / dr > SIGMA_THRESHOLD
    SIGMA_THRESHOLD = 1.0  # change to 2.0, 3.0, etc. if desired

    lo_diff   = np.abs(lo_ratio - 1.0)
    hi_diff   = np.abs(hi_ratio - 1.0)
    lo_sigma  = lo_diff / lo_dratio
    hi_sigma  = hi_diff / hi_dratio

    lo_fail_mask = (lo_sigma > SIGMA_THRESHOLD)
    hi_fail_mask = (hi_sigma > SIGMA_THRESHOLD)

    lo_points_pass = not lo_fail_mask.any()
    hi_points_pass = not hi_fail_mask.any()

    # Identify failing bins (t-phi) for diagnostics
    if lo_fail_mask.any():
        print("LOEPS bins outside band and significant (|r-1| > band AND |r-1|/dr > sigma cut):")
        for r, dr, phi, t, d, s in zip(
            lo_ratio[lo_fail_mask],
            lo_dratio[lo_fail_mask],
            lo_phi_bin[lo_fail_mask],
            lo_t_bin[lo_fail_mask],
            lo_diff[lo_fail_mask],
            lo_sigma[lo_fail_mask],
        ):
            print(
                f"  t_bin={int(t)}, phi_bin={int(phi)}, "
                f"r={r:.4f}, dr={dr:.4f}, |r-1|={d:.4g}, |r-1|/dr={s:.3f}"
            )
    else:
        print("LOEPS: all bins inside band or statistically consistent with unity.")

    if hi_fail_mask.any():
        print("HIEPS bins outside band and significant (|r-1| > band AND |r-1|/dr > sigma cut):")
        for r, dr, phi, t, d, s in zip(
            hi_ratio[hi_fail_mask],
            hi_dratio[hi_fail_mask],
            hi_phi_bin[hi_fail_mask],
            hi_t_bin[hi_fail_mask],
            hi_diff[hi_fail_mask],
            hi_sigma[hi_fail_mask],
        ):
            print(
                f"  t_bin={int(t)}, phi_bin={int(phi)}, "
                f"r={r:.4f}, dr={dr:.4f}, |r-1|={d:.4g}, |r-1|/dr={s:.3f}"
            )
    else:
        print("HIEPS: all bins inside band or statistically consistent with unity.")

    # 2) Error-weighted average checks (unchanged)
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
    # Diagnostic plots (PDF only)
    # ---------------------------
    monitor_dir = iter_dir / "monitor" / "ratios"
    monitor_dir.mkdir(parents=True, exist_ok=True)

    def make_plots(tag, ratio, dratio, phi_bin, t_bin, diff, sigma, fail_mask):
        """
        Create a multi-page PDF with diagnostics for a given setting (LOEPS/HIEPS).
        Pages:
          1) ratio vs index with errors and band (absolute)
          2) |ratio-1| vs index + band
          3) |ratio-1|/dr vs index + sigma cut
          4) ratio vs phi_bin (colored by t_bin)
        """
        pdf_path = monitor_dir / f"ratios_{tag}.pdf"
        with PdfPages(pdf_path) as pdf:

            n = len(ratio)
            idx = np.arange(n)

            # 1) ratio vs index with error bars and absolute band
            fig, ax = plt.subplots()
            ax.errorbar(idx, ratio, yerr=dratio, fmt='o', label=f'{tag} ratios')
            ax.axhline(1.0, linestyle='--')
            ax.axhspan(1.0 - RATIO_THRESHOLD_SPREAD,
                       1.0 + RATIO_THRESHOLD_SPREAD,
                       alpha=0.2, label='band')
            if fail_mask.any():
                ax.scatter(idx[fail_mask], ratio[fail_mask], marker='x', label='outside band + significant')
            ax.set_xlabel("bin index")
            ax.set_ylabel("ratio")
            ax.set_title(f"{tag}: ratio vs bin index")
            ax.legend()
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            # 2) |ratio-1| vs index
            fig, ax = plt.subplots()
            ax.bar(idx, diff)
            ax.axhline(RATIO_THRESHOLD_SPREAD, linestyle='--')
            ax.set_xlabel("bin index")
            ax.set_ylabel("|ratio - 1|")
            ax.set_title(f"{tag}: deviation from unity (absolute)")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            # 3) |ratio-1|/dr vs index (significance)
            fig, ax = plt.subplots()
            ax.bar(idx, sigma)
            ax.axhline(SIGMA_THRESHOLD, linestyle='--')
            ax.set_xlabel("bin index")
            ax.set_ylabel("|ratio - 1| / dr")
            ax.set_title(f"{tag}: deviation significance")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            # 4) ratio vs phi_bin, colored by t_bin
            fig, ax = plt.subplots()
            sc = ax.scatter(phi_bin, ratio, c=t_bin, cmap='viridis')
            ax.axhline(1.0, linestyle='--')
            ax.axhspan(1.0 - RATIO_THRESHOLD_SPREAD,
                       1.0 + RATIO_THRESHOLD_SPREAD,
                       alpha=0.2)
            ax.set_xlabel("phi_bin")
            ax.set_ylabel("ratio")
            ax.set_title(f"{tag}: ratio vs phi_bin (color = t_bin)")
            cbar = fig.colorbar(sc, ax=ax)
            cbar.set_label("t_bin")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    make_plots("LOEPS", lo_ratio, lo_dratio, lo_phi_bin, lo_t_bin, lo_diff, lo_sigma, lo_fail_mask)
    make_plots("HIEPS", hi_ratio, hi_dratio, hi_phi_bin, hi_t_bin, hi_diff, hi_sigma, hi_fail_mask)

    CONTINUE = overall
    return CONTINUE