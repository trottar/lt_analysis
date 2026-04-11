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
from ROOT import TFile, TH1F

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
    add_derived_branches_to_file,
    add_shift_branch_to_file,
    build_shifted_hist_from_hist,
    compute_mm_shift_details_from_hists,
    get_hist_nbins,
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


def _make_hist(hist_name, hist_title, hist_xmin, hist_xmax):
    hist_nbins = get_hist_nbins(hist_xmin, hist_xmax)
    hist = TH1F(hist_name, hist_title, hist_nbins, hist_xmin, hist_xmax)
    hist.SetDirectory(0)
    return hist


def _open_tree(filename, tree_name):
    root_file = TFile.Open(filename, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError("Unable to open ROOT file: {}".format(filename))

    tree = root_file.Get(tree_name)
    if not tree:
        root_file.Close()
        raise RuntimeError("Tree '{}' not found in {}".format(tree_name, filename))

    return root_file, tree


def _build_cut_data_hists(phi_setting, particle_type, data_filename, inpDict, mm_min, mm_max, tmin, tmax):
    cuts_path = os.path.dirname(__file__)
    if cuts_path not in sys.path:
        sys.path.append(cuts_path)

    from apply_cuts import apply_data_sub_cuts, set_shift_context, set_val

    set_val(inpDict)
    set_shift_context(phi_setting=phi_setting, shift_mode="raw", mm_shift_summary={}, t_shift_summary={})

    hole_contains = None
    if particle_type == "kaon":
        from hgcer_hole import apply_HGCer_hole_cut
        hgcer_cutg = apply_HGCer_hole_cut(inpDict["Q2"], inpDict["W"], inpDict["EPSSET"])
        hole_contains = hgcer_cutg.IsInside

    tree_name = "Cut_{}_Events_prompt_noRF".format(particle_type.capitalize())
    root_file, tree = _open_tree(data_filename, tree_name)

    mm_hist = _make_hist("H_MM_DATA_shiftprep", "Shift-prep data MM", mm_min, mm_max)
    t_hist = _make_hist("H_t_DATA_shiftprep", "Shift-prep data -t", tmin, tmax)

    try:
        for evt in tree:
            if not apply_data_sub_cuts(evt, mm_min, mm_max):
                continue
            if hole_contains and hole_contains(evt.P_hgcer_xAtCer, evt.P_hgcer_yAtCer):
                continue

            mm_hist.Fill(evt.MM)
            t_hist.Fill(-evt.MandelT)
    finally:
        root_file.Close()

    return mm_hist, t_hist


def _build_cut_simc_hists(phi_setting, simc_filename, inpDict, mm_min, mm_max, tmin, tmax):
    cuts_path = os.path.dirname(__file__)
    if cuts_path not in sys.path:
        sys.path.append(cuts_path)

    from apply_cuts import apply_simc_sub_cuts, set_shift_context, set_val

    set_val(inpDict)
    set_shift_context(phi_setting=phi_setting, shift_mode="raw", mm_shift_summary={}, t_shift_summary={})

    root_file, tree = _open_tree(simc_filename, "h10")

    mm_hist = _make_hist("H_MM_SIMC_shiftprep", "Shift-prep SIMC MM", mm_min, mm_max)
    t_hist = _make_hist("H_t_SIMC_shiftprep", "Shift-prep SIMC -t", tmin, tmax)

    try:
        for evt in tree:
            if not apply_simc_sub_cuts(evt, mm_min, mm_max):
                continue

            mm_hist.Fill(evt.missmass)
            t_hist.Fill(-evt.t)
    finally:
        root_file.Close()

    return mm_hist, t_hist


def _apply_shift_branches(particle_type, data_filename, dummy_filename, mm_shift, t_shift=None):
    tree_names = ["Cut_{}_Events_prompt_noRF".format(particle_type.capitalize())]

    if t_shift is not None:
        shift_branch_specs = [
            ("MM_shift", lambda evt, shift=mm_shift: getattr(evt, "MM") + shift, "Applying"),
            ("t_shift", lambda evt, shift=t_shift: -getattr(evt, "MandelT") + shift, "Applying"),
        ]
        add_derived_branches_to_file(data_filename, tree_names, shift_branch_specs)
        if dummy_filename and os.path.exists(dummy_filename):
            add_derived_branches_to_file(dummy_filename, tree_names, shift_branch_specs)
        return

    add_shift_branch_to_file(data_filename, tree_names, "MM", mm_shift)
    if dummy_filename and os.path.exists(dummy_filename):
        add_shift_branch_to_file(dummy_filename, tree_names, "MM", mm_shift)


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


def shift_prep(phi_setting, inpDict):
    particle_type = inpDict["ParticleType"]
    q2 = _to_float(inpDict["Q2"])
    w = _to_float(inpDict["W"])
    mm_min = float(inpDict["mm_min"])
    mm_max = float(inpDict["mm_max"])
    tmin = float(inpDict["tmin"])
    tmax = float(inpDict["tmax"])

    print("\nRunning shift_prep for {} {}...".format(phi_setting, particle_type))

    data_filename = "{}/{}_{}_{}.root".format(
        OUTPATH,
        phi_setting,
        particle_type,
        inpDict["InDATAFilename"],
    )
    if not os.path.exists(data_filename):
        return {
            "phi_setting": phi_setting,
            "mm_shift": None,
            "t_shift": None,
            "artifacts": {},
        }

    dummy_filename = "{}/{}_{}_{}.root".format(
        OUTPATH,
        phi_setting,
        particle_type,
        inpDict["InDUMMYFilename"],
    )
    simc_filename = "{}/Prod_Coin_Q{}W{}{}_{}e.root".format(
        OUTPATH,
        inpDict["Q2"],
        inpDict["W"],
        phi_setting.lower(),
        inpDict["EPSSET"],
    )
    if not os.path.exists(simc_filename):
        return {
            "phi_setting": phi_setting,
            "mm_shift": None,
            "t_shift": None,
            "artifacts": {},
        }

    data_mm_hist, data_t_hist = _build_cut_data_hists(
        phi_setting,
        particle_type,
        data_filename,
        inpDict,
        mm_min,
        mm_max,
        tmin,
        tmax,
    )
    simc_mm_hist, simc_t_hist = _build_cut_simc_hists(
        phi_setting,
        simc_filename,
        inpDict,
        mm_min,
        mm_max,
        tmin,
        tmax,
    )

    mm_details = compute_mm_shift_details_from_hists(
        particle_type,
        simc_mm_hist,
        data_mm_hist,
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
            simc_t_hist,
            data_t_hist,
            t_shift,
            t_plot_filename,
            plot_xmin=tmin,
            plot_xmax=tmax,
        )
        shifted_t_hist = build_shifted_hist_from_hist(
            data_t_hist,
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
        data_mm_hist,
        mm_shift,
        hist_name="H_MM_DATA_shifted_shiftprep",
        hist_title="Shift-prep shifted MM",
    )

    _apply_shift_branches(
        particle_type,
        data_filename,
        dummy_filename,
        mm_shift,
        t_shift_entry.get("shift"),
    )

    _write_shift_prep_root(
        output_root,
        {
            "H_MM_DATA_shiftprep": data_mm_hist,
            "H_MM_SIMC_shiftprep": simc_mm_hist,
            "H_MM_DATA_shifted_shiftprep": shifted_mm_hist,
            "H_t_DATA_shiftprep": data_t_hist,
            "H_t_SIMC_shiftprep": simc_t_hist,
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
    }
