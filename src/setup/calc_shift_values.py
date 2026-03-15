#! /usr/bin/python

#
# Description:
# ================================================================
# Compute MM and t shift values without modifying ROOT files.
# ================================================================
#

import json
import math
import re
import subprocess
import sys

from shift_MM import compute_mm_shift_summary


def to_float(value):
    return float(str(value).replace("p", "."))


def run_t_shift_program(executable, q2, w, theta_cm_deg, mm_shift_apply, beam_energy_gev):
    # tshift_kaonff.f expects the opposite MM-shift convention:
    # MM_expt = MM_correct + MMshift.
    mm_shift_fortran_mev = -1000.0 * mm_shift_apply
    beam_energy_mev = 1000.0 * beam_energy_gev

    input_lines = [
        f"{q2:.6f} {w:.6f} {theta_cm_deg:.6f} {mm_shift_fortran_mev:.6f} {beam_energy_mev:.6f}",
    ]

    process = subprocess.run(
        [executable],
        input="\n".join(input_lines) + "\n",
        capture_output=True,
        text=True,
        check=False,
    )

    if process.returncode != 0:
        raise RuntimeError(
            "tshift_kaonff failed with exit code {}.\nstdout:\n{}\nstderr:\n{}".format(
                process.returncode,
                process.stdout,
                process.stderr,
            )
        )

    output = process.stdout
    match = re.search(
        r"TSHIFT_GEV2\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)",
        output,
    )
    if not match:
        raise RuntimeError(f"Unable to parse t-shift output:\n{output}")

    return {
        "fortran_mm_shift_mev": mm_shift_fortran_mev,
        # This is already the shift in the code's -t convention.
        "t_shift": float(match.group(1)),
        "phi_deg": 0.0,
        "raw_output": output,
    }


def main():
    if len(sys.argv) != 11:
        raise SystemExit(
            "Usage: calc_shift_values.py <particle_type> <simc_root> <data_root> "
            "<q2> <w> <theta_cm_deg> <beam_energy_gev> <mm_min> <mm_max> <tshift_executable_or_none>"
        )

    particle_type = sys.argv[1]
    simc_root = sys.argv[2]
    data_root = sys.argv[3]
    q2 = to_float(sys.argv[4])
    w = to_float(sys.argv[5])
    theta_cm_deg = float(sys.argv[6])
    beam_energy_gev = float(sys.argv[7])
    mm_min = float(sys.argv[8])
    mm_max = float(sys.argv[9])
    tshift_executable = sys.argv[10]

    mm_summary = compute_mm_shift_summary(
        particle_type,
        simc_root,
        data_root,
        plot_xmin=mm_min,
        plot_xmax=mm_max,
    )

    output = dict(mm_summary)
    output.update(
        {
            "particle_type": particle_type,
            "q2": q2,
            "w": w,
            "theta_cm_deg": theta_cm_deg,
            "beam_energy_gev": beam_energy_gev,
            "t_shift": None,
            "phi_deg": None,
        }
    )

    if particle_type == "kaon" and tshift_executable.lower() != "none":
        t_summary = run_t_shift_program(
            tshift_executable,
            q2,
            w,
            theta_cm_deg,
            mm_summary["shift"],
            beam_energy_gev,
        )
        output["t_shift"] = t_summary["t_shift"]
        output["phi_deg"] = t_summary["phi_deg"]
        output["fortran_mm_shift_mev"] = t_summary["fortran_mm_shift_mev"]

    print(json.dumps(output))


if __name__ == "__main__":
    main()
