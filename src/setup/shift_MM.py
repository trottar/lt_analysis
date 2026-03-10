#! /usr/bin/python

#
# Description:
# ================================================================
# Shift experimental missing-mass branches to the fitted SIMC peak.
# ================================================================
#

import os
import sys
from array import array

import ROOT
from ROOT import TFile, TH1F, TTree, kRed


def get_peak_window(particle_type):
    if particle_type == "kaon":
        peak_name = "lambda"
        peak_center = 1.1156
    else:
        peak_name = "proton"
        peak_center = 0.9380

    return peak_name, peak_center - 0.05, peak_center + 0.05


def get_data_tree_names(particle_type):
    particle_name = particle_type.capitalize()
    return [
        f"Uncut_{particle_name}_Events",
        f"Cut_{particle_name}_Events_all_noRF",
        f"Cut_{particle_name}_Events_prompt_noRF",
        f"Cut_{particle_name}_Events_rand_noRF",
        f"Cut_{particle_name}_Events_all_RF",
        f"Cut_{particle_name}_Events_prompt_RF",
        f"Cut_{particle_name}_Events_rand_RF",
    ]


def fit_gaussian(hist, x_min, x_max):
    bin_min = hist.GetXaxis().FindBin(x_min)
    bin_max = hist.GetXaxis().FindBin(x_max)

    max_bin = bin_min
    max_value = hist.GetBinContent(max_bin)
    for bin_idx in range(bin_min, bin_max + 1):
        bin_value = hist.GetBinContent(bin_idx)
        if bin_value > max_value:
            max_bin = bin_idx
            max_value = bin_value

    half_max = max_value * 0.75

    left_bin = max_bin
    right_bin = max_bin
    while left_bin > 1 and hist.GetBinContent(left_bin) > half_max:
        left_bin -= 1
    while right_bin < hist.GetNbinsX() and hist.GetBinContent(right_bin) > half_max:
        right_bin += 1

    fit_min = hist.GetBinCenter(left_bin)
    fit_max = hist.GetBinCenter(right_bin)

    hist.Fit("gaus", "Q", "", fit_min, fit_max)
    fit_func = hist.GetFunction("gaus")
    if not fit_func:
        raise RuntimeError("Gaussian fit failed.")

    fit_func.SetLineColor(kRed)
    return fit_func.GetParameter(1), fit_func.GetParError(1)


def fit_tree_peak(filename, tree_name, branch_name, particle_type):
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    root_file = TFile.Open(filename, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Unable to open ROOT file: {filename}")

    tree = root_file.Get(tree_name)
    if not tree:
        root_file.Close()
        raise RuntimeError(f"Tree '{tree_name}' not found in {filename}")

    _, mm_min, mm_max = get_peak_window(particle_type)
    hist_name = f"hist_{abs(hash((filename, tree_name, branch_name))) & 0xFFFFFFFF}"
    hist = TH1F(hist_name, branch_name, 200, 0.7, 1.5)
    tree.Draw(f"{branch_name}>>{hist_name}", "", "goff")
    hist.SetDirectory(0)

    mean, mean_err = fit_gaussian(hist, mm_min, mm_max)
    root_file.Close()
    return mean, mean_err


def build_shifted_tree(tree, source_branch_name, shift_branch_name, shift):
    shifted_value = array("f", [0.0])
    new_tree = tree.CloneTree(0)
    new_tree.Branch(shift_branch_name, shifted_value, f"{shift_branch_name}/F")

    for evt in tree:
        shifted_value[0] = getattr(evt, source_branch_name) + shift
        new_tree.Fill()

    return new_tree


def add_shift_branch_to_file(filename, tree_names, source_branch_name, shift):
    shift_branch_name = f"{source_branch_name}_shift"
    root_file = TFile.Open(filename, "UPDATE")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Unable to open ROOT file for update: {filename}")

    updated_trees = []
    for tree_name in tree_names:
        tree = root_file.Get(tree_name)
        if not tree:
            print(f"Tree '{tree_name}' not found in {filename}. Skipping.")
            continue

        if tree.GetBranch(shift_branch_name):
            print(f"{filename}:{tree_name} already has '{shift_branch_name}'. Leaving it unchanged.")
            continue

        print(f"Applying shift {shift:+.6f} to {filename}:{tree_name}")
        updated_trees.append(
            (tree_name, build_shifted_tree(tree, source_branch_name, shift_branch_name, shift))
        )

    root_file.cd()
    for tree_name, new_tree in updated_trees:
        new_tree.Write(tree_name, TTree.kOverwrite)

    root_file.Close()
    return len(updated_trees)


def shift_experimental_files_to_simc_peak(particle_type, simc_filename, data_filename, dummy_filename=None):
    ROOT.gROOT.SetBatch(True)

    reference_tree = f"Cut_{particle_type.capitalize()}_Events_prompt_noRF"
    peak_name, _, _ = get_peak_window(particle_type)

    simc_peak, simc_peak_err = fit_tree_peak(simc_filename, "h10", "missmass", particle_type)
    data_peak, data_peak_err = fit_tree_peak(data_filename, reference_tree, "MM", particle_type)
    shift = simc_peak - data_peak

    print(f"\nAligning experimental {peak_name} peak to SIMC...")
    print(f"SIMC peak ({simc_filename}) = {simc_peak:.6f} +/- {simc_peak_err:.6f}")
    print(f"Data peak ({data_filename}) = {data_peak:.6f} +/- {data_peak_err:.6f}")
    print(f"Applying MM shift = {shift:+.6f}")

    tree_names = get_data_tree_names(particle_type)
    add_shift_branch_to_file(data_filename, tree_names, "MM", shift)

    if dummy_filename:
        if os.path.exists(dummy_filename):
            add_shift_branch_to_file(dummy_filename, tree_names, "MM", shift)
        else:
            print(f"Dummy file '{dummy_filename}' not found. Skipping dummy MM shift.")

    return {
        "simc_peak": simc_peak,
        "simc_peak_err": simc_peak_err,
        "data_peak": data_peak,
        "data_peak_err": data_peak_err,
        "shift": shift,
    }


def main():
    if len(sys.argv) not in (4, 5):
        print(
            "Usage: shift_MM.py <particle_type> <simc_root> <data_root> [dummy_root]"
        )
        sys.exit(1)

    particle_type = sys.argv[1]
    simc_filename = sys.argv[2]
    data_filename = sys.argv[3]
    dummy_filename = sys.argv[4] if len(sys.argv) == 5 else None

    shift_experimental_files_to_simc_peak(
        particle_type,
        simc_filename,
        data_filename,
        dummy_filename,
    )


if __name__ == "__main__":
    main()
