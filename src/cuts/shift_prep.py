#! /usr/bin/python

#
# Description:
# ================================================================
# Build clean pre-rand-sub MM samples to determine MM/t shifts.
# ================================================================
#

import json
import os
import subprocess
import sys

import ROOT
from ROOT import TFile

from ltsep import Root

lt = Root(os.path.realpath(__file__), "Plot_LTSep")

USER = lt.USER
HOST = lt.HOST
REPLAYPATH = lt.REPLAYPATH
UTILPATH = lt.UTILPATH
LTANAPATH = lt.LTANAPATH
ANATYPE = lt.ANATYPE
OUTPATH = lt.OUTPATH

SETUP_PATH = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "setup"))
SIMC_PATH = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "simc_ana"))
if SETUP_PATH not in sys.path:
    sys.path.append(SETUP_PATH)
if SIMC_PATH not in sys.path:
    sys.path.append(SIMC_PATH)

from calc_shift_values import run_t_shift_program
from shift_MM import (
    build_shifted_hist_from_hist,
    compute_mm_shift_details_from_hists,
    plot_mm_shift_from_hist_details,
    write_t_shift_plots_from_hists,
)


def _to_float(value):
    return float(str(value).replace("p", "."))


def _get_output_base(phi_setting, particle_type, inpDict):
    return "{}/{}_{}_ShiftPrep_Q{}W{}_{}e".format(
        OUTPATH,
        phi_setting,
        particle_type,
        inpDict["Q2"],
        inpDict["W"],
        inpDict["EPSSET"],
    )


def _clone_for_output(hist, clone_name):
    cloned_hist = hist.Clone(clone_name)
    if hasattr(cloned_hist, "SetDirectory"):
        cloned_hist.SetDirectory(0)
    return cloned_hist


def _write_shift_prep_root(output_root, hist_map):
    root_file = TFile.Open(output_root, "RECREATE")
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Unable to create shift-prep ROOT file: {}".format(output_root))

    try:
        for hist_name, hist in hist_map.items():
            if hist is None:
                continue
            cloned_hist = _clone_for_output(hist, hist_name)
            cloned_hist.Write()
    finally:
        root_file.Close()


def _ensure_tshift_executable(particle_type):
    if particle_type != "kaon":
        return None

    src_dir = os.path.join(LTANAPATH, "src")
    exe_name = "tshift_kaonff.exe" if os.name == "nt" else "tshift_kaonff"
    exe_path = os.path.join(src_dir, exe_name)
    src_path = os.path.join(src_dir, "tshift_kaonff.f")
    inc_path = os.path.join(src_dir, "linux_suppl.inc")

    needs_rebuild = (
        (not os.path.exists(exe_path))
        or os.path.getmtime(src_path) > os.path.getmtime(exe_path)
        or os.path.getmtime(inc_path) > os.path.getmtime(exe_path)
    )
    if not needs_rebuild:
        return exe_path

    print("Compiling tshift_kaonff.f...")
    result = subprocess.run(
        ["gfortran", "-o", exe_name, "tshift_kaonff.f"],
        cwd=src_dir,
        capture_output=True,
        text=True,
        check=False,
    )
    if result.stdout:
        print(result.stdout.rstrip())
    if result.stderr:
        print(result.stderr.rstrip())
    if result.returncode != 0:
        raise RuntimeError(
            "Unable to compile tshift_kaonff.f (exit code {}).".format(result.returncode)
        )

    return exe_path


def get_shift_theta_and_beam(phi_setting, inpDict):
    out_filename = inpDict["OutFilename"]
    particle_type = inpDict["ParticleType"]

    config = {
        "Right": {
            "run_nums": [run for run in inpDict["runNumRight"].split(" ") if run],
            "ptheta_vals": inpDict["pThetaValRight"],
            "beam_vals": inpDict["EbeamValRight"],
        },
        "Left": {
            "run_nums": [run for run in inpDict["runNumLeft"].split(" ") if run],
            "ptheta_vals": inpDict["pThetaValLeft"],
            "beam_vals": inpDict["EbeamValLeft"],
        },
        "Center": {
            "run_nums": [run for run in inpDict["runNumCenter"].split(" ") if run],
            "ptheta_vals": inpDict["pThetaValCenter"],
            "beam_vals": inpDict["EbeamValCenter"],
        },
    }

    setting = config[phi_setting]
    run_nums = setting["run_nums"]
    ptheta_vals = setting["ptheta_vals"]
    beam_vals = setting["beam_vals"]
    max_len = min(len(run_nums), len(ptheta_vals), len(beam_vals))
    if max_len == 0:
        raise RuntimeError("No valid theta/beam entries found for {} shift calculation.".format(phi_setting))

    chosen_index = 0
    log_suffix = out_filename.replace("FullAnalysis_", "")
    for idx, run in enumerate(run_nums[:max_len]):
        pid_log = "{}/log/{}_{}_{}_{}.log".format(
            LTANAPATH,
            phi_setting,
            particle_type,
            run,
            log_suffix,
        )
        if os.path.exists(pid_log):
            chosen_index = idx
            break

    return float(ptheta_vals[chosen_index]), float(beam_vals[chosen_index])


def shift_prep(phi_setting, hist_seed, inpDict):
    particle_type = inpDict["ParticleType"]
    q2 = _to_float(inpDict["Q2"])
    w = _to_float(inpDict["W"])
    mm_min = float(inpDict["mm_min"])
    mm_max = float(inpDict["mm_max"])
    tmin = float(inpDict["tmin"])
    tmax = float(inpDict["tmax"])

    sys.path.append(os.path.dirname(__file__))
    from rand_sub import rand_sub
    from compare_simc import compare_simc

    print("\nRunning shift_prep for {} {}...".format(phi_setting, particle_type))

    raw_hist = rand_sub(phi_setting, inpDict, shift_mode="raw", emit_plots=False)
    if len(raw_hist.keys()) <= 1:
        return {
            "phi_setting": phi_setting,
            "mm_shift": None,
            "t_shift": None,
            "artifacts": {},
        }

    simc_request = {
        "phi_setting": phi_setting,
        "normfac_simc": hist_seed["normfac_simc"],
        "InSIMCFilename": hist_seed["InSIMCFilename"],
    }
    simc_hist = compare_simc(simc_request, inpDict, emit_plots=False)

    mm_details = compute_mm_shift_details_from_hists(
        particle_type,
        simc_hist["H_MM_SIMC"],
        raw_hist["H_MM_DATA"],
        plot_xmin=mm_min,
        plot_xmax=mm_max,
    )
    mm_shift = float(mm_details["shift"])

    theta_cm_deg, beam_energy_gev = get_shift_theta_and_beam(phi_setting, inpDict)
    output_base = _get_output_base(phi_setting, particle_type, inpDict)
    mm_plot_filename = output_base + "_MM_Shift.pdf"
    t_plot_filename = output_base + "_T_Shift.pdf"
    output_root = output_base + ".root"
    output_json = output_base + ".json"

    plot_mm_shift_from_hist_details(
        particle_type,
        mm_details,
        mm_plot_filename,
        plot_xmin=mm_min,
        plot_xmax=mm_max,
    )

    mm_shift_entry = {
        "simc_peak": float(mm_details["simc_fit"]["mean"]),
        "simc_peak_err": float(mm_details["simc_fit"]["mean_err"]),
        "data_peak": float(mm_details["data_fit"]["mean"]),
        "data_peak_err": float(mm_details["data_fit"]["mean_err"]),
        "simc_sigma": float(mm_details["simc_fit"]["sigma"]),
        "data_sigma": float(mm_details["data_fit"]["sigma"]),
        "simc_fit_min": float(mm_details["simc_fit"]["fit_min"]),
        "simc_fit_max": float(mm_details["simc_fit"]["fit_max"]),
        "data_fit_min": float(mm_details["data_fit"]["fit_min"]),
        "data_fit_max": float(mm_details["data_fit"]["fit_max"]),
        "shift": mm_shift,
        "theta_cm_deg": theta_cm_deg,
        "beam_energy_gev": beam_energy_gev,
        "plot_filename": mm_plot_filename,
    }

    t_shift_entry = {}
    shifted_t_hist = None
    if particle_type == "kaon":
        tshift_executable = _ensure_tshift_executable(particle_type)
        t_summary = run_t_shift_program(
            tshift_executable,
            q2,
            w,
            theta_cm_deg,
            mm_shift,
            beam_energy_gev,
        )
        t_shift = float(t_summary["t_shift"])
        write_t_shift_plots_from_hists(
            particle_type,
            simc_hist["H_t_SIMC"],
            raw_hist["H_t_DATA"],
            t_shift,
            t_plot_filename,
            plot_xmin=tmin,
            plot_xmax=tmax,
        )
        shifted_t_hist = build_shifted_hist_from_hist(
            raw_hist["H_t_DATA"],
            t_shift,
            hist_name="H_t_DATA_shifted_shiftprep",
            hist_title="Shift-prep shifted -t",
        )
        t_shift_entry = {
            "shift": t_shift,
            "theta_cm_deg": theta_cm_deg,
            "beam_energy_gev": beam_energy_gev,
            "phi_deg": float(t_summary.get("phi_deg", 0.0)),
            "plot_filename": t_plot_filename,
        }

    shifted_mm_hist = build_shifted_hist_from_hist(
        raw_hist["H_MM_DATA"],
        mm_shift,
        hist_name="H_MM_DATA_shifted_shiftprep",
        hist_title="Shift-prep shifted MM",
    )
    _write_shift_prep_root(
        output_root,
        {
            "H_MM_DATA_shiftprep": raw_hist["H_MM_DATA"],
            "H_MM_SIMC_shiftprep": simc_hist["H_MM_SIMC"],
            "H_MM_DATA_shifted_shiftprep": shifted_mm_hist,
            "H_t_DATA_shiftprep": raw_hist["H_t_DATA"],
            "H_t_SIMC_shiftprep": simc_hist["H_t_SIMC"],
            "H_t_DATA_shifted_shiftprep": shifted_t_hist,
        },
    )

    summary_payload = {
        "phi_setting": phi_setting,
        "particle_type": particle_type,
        "mm_shift": mm_shift_entry,
        "t_shift": t_shift_entry,
        "artifacts": {
            "json_filename": output_json,
            "root_filename": output_root,
            "mm_plot_filename": mm_plot_filename,
            "t_plot_filename": t_shift_entry.get("plot_filename"),
        },
    }
    with open(output_json, "w") as output_file:
        json.dump(summary_payload, output_file, indent=2, sort_keys=True)

    return {
        "phi_setting": phi_setting,
        "mm_shift": mm_shift_entry,
        "t_shift": t_shift_entry,
        "artifacts": summary_payload["artifacts"],
        "raw_hist": raw_hist,
        "simc_hist": simc_hist,
    }
