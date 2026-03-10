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
from ROOT import TCanvas, TFile, TF1, TH1F, TLegend, TLine, TTree, gStyle, kBlack, kBlue, kGreen, kRed


HIST_NBINS = 200
FIT_HIST_XMIN = 0.7
FIT_HIST_XMAX = 1.5


def get_peak_window(particle_type):
    if particle_type == "kaon":
        peak_name = "lambda"
        peak_center = 1.1156
    else:
        peak_name = "proton"
        peak_center = 0.9380

    return peak_name, peak_center - 0.05, peak_center + 0.05


def get_peak_center(particle_type):
    if particle_type == "kaon":
        return 1.1156
    return 0.9380


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


def fit_gaussian(hist, peak_center, fit_min, fit_max):
    fit_name = f"gaus_{hist.GetName()}"
    fit_func = TF1(fit_name, "gaus", fit_min, fit_max)
    peak_bin = hist.GetXaxis().FindBin(peak_center)
    peak_height = hist.GetBinContent(peak_bin)
    if peak_height <= 0.0:
        peak_height = hist.GetMaximum()
    if peak_height <= 0.0:
        raise RuntimeError("Histogram is empty in the requested fit window.")

    fit_func.SetParameters(peak_height, peak_center, max((fit_max - fit_min) / 6.0, 0.002))
    fit_func.SetParLimits(0, 0.0, peak_height * 10.0)
    fit_func.SetParLimits(1, fit_min, fit_max)
    fit_func.SetParLimits(2, 0.001, max(fit_max - fit_min, 0.01))

    hist.Fit(fit_func, "QR")
    if not fit_func:
        raise RuntimeError("Gaussian fit failed.")

    fit_func.SetLineColor(kRed)
    return fit_func.GetParameter(1), fit_func.GetParError(1)


def build_histogram(
    filename,
    tree_name,
    branch_name,
    hist_name,
    hist_title,
    shift=0.0,
    hist_xmin=FIT_HIST_XMIN,
    hist_xmax=FIT_HIST_XMAX,
):
    if not os.path.exists(filename):
        raise FileNotFoundError(filename)

    root_file = TFile.Open(filename, "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Unable to open ROOT file: {filename}")

    tree = root_file.Get(tree_name)
    if not tree:
        root_file.Close()
        raise RuntimeError(f"Tree '{tree_name}' not found in {filename}")

    expression = branch_name
    if abs(shift) > 0.0:
        expression = f"({branch_name}+({shift:.12g}))"

    hist = TH1F(hist_name, hist_title, HIST_NBINS, hist_xmin, hist_xmax)
    tree.Draw(f"{expression}>>{hist_name}", "", "goff")
    hist.SetDirectory(0)
    root_file.Close()
    return hist


def fit_tree_peak(
    filename,
    tree_name,
    branch_name,
    particle_type,
    hist_xmin=FIT_HIST_XMIN,
    hist_xmax=FIT_HIST_XMAX,
):
    peak_center = get_peak_center(particle_type)
    fit_min = hist_xmin
    fit_max = hist_xmax
    if not (fit_min <= peak_center <= fit_max):
        raise RuntimeError(
            f"Expected peak {peak_center:.4f} is outside histogram range [{fit_min:.4f}, {fit_max:.4f}]."
        )

    hist_name = f"hist_{abs(hash((filename, tree_name, branch_name))) & 0xFFFFFFFF}"
    hist = build_histogram(
        filename,
        tree_name,
        branch_name,
        hist_name,
        branch_name,
        hist_xmin=hist_xmin,
        hist_xmax=hist_xmax,
    )
    mean, mean_err = fit_gaussian(hist, peak_center, fit_min, fit_max)
    return {
        "hist": hist,
        "mean": mean,
        "mean_err": mean_err,
    }


def build_shifted_tree(tree, source_branch_name, shift_branch_name, shift):
    shifted_value = array("f", [0.0])

    # Clone the tree without the old shift branch so the replacement branch
    # is freshly written with the current alignment.
    has_existing_shift_branch = bool(tree.GetBranch(shift_branch_name))
    if has_existing_shift_branch:
        tree.SetBranchStatus(shift_branch_name, 0)

    new_tree = tree.CloneTree(0)

    if has_existing_shift_branch:
        tree.SetBranchStatus(shift_branch_name, 1)

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

        action = "Replacing" if tree.GetBranch(shift_branch_name) else "Applying"
        print(f"{action} shift {shift:+.6f} in {filename}:{tree_name}")
        updated_trees.append(
            (tree_name, build_shifted_tree(tree, source_branch_name, shift_branch_name, shift))
        )

    root_file.cd()
    for tree_name, new_tree in updated_trees:
        new_tree.Write(tree_name, TTree.kOverwrite)

    root_file.Close()
    return len(updated_trees)


def make_text_box(lines):
    text_box = ROOT.TPaveText(0.58, 0.70, 0.89, 0.89, "NDC")
    text_box.SetFillColor(0)
    text_box.SetBorderSize(1)
    text_box.SetTextSize(0.03)
    text_box.SetTextAlign(12)
    for line in lines:
        text_box.AddText(line)
    return text_box


def make_peak_line(peak_value, ymax, color):
    line = TLine(peak_value, 0.0, peak_value, ymax)
    line.SetLineColor(color)
    line.SetLineStyle(2)
    line.SetLineWidth(3)
    return line


def normalize_hist(hist):
    integral = hist.Integral()
    if integral > 0.0:
        hist.Scale(1.0 / integral)


def make_plot_filename(data_filename):
    base_name = os.path.splitext(os.path.basename(data_filename))[0]
    return os.path.join(os.path.dirname(data_filename), f"{base_name}_MM_Shift.pdf")


def write_shift_plots(
    particle_type,
    simc_fit,
    data_fit,
    shifted_hist,
    shift,
    output_pdf,
    plot_xmin=FIT_HIST_XMIN,
    plot_xmax=FIT_HIST_XMAX,
):
    old_opt_stat = gStyle.GetOptStat()
    gStyle.SetOptStat(0)

    canvas_name = f"canvas_mm_shift_{abs(hash(output_pdf)) & 0xFFFFFFFF}"
    canvas = TCanvas(canvas_name, "MM Shift", 900, 700)
    canvas.Print(f"{output_pdf}[")

    simc_hist = simc_fit["hist"].Clone(f"{simc_fit['hist'].GetName()}_plot")
    simc_hist.SetTitle(f"SIMC M_{particle_type[0].upper()} Fit")
    simc_hist.SetLineColor(kRed)
    simc_hist.SetLineWidth(2)
    simc_hist.GetXaxis().SetRangeUser(plot_xmin, plot_xmax)
    simc_hist.Draw("hist")
    simc_fit_line = simc_hist.GetFunction("gaus")
    if simc_fit_line:
        simc_fit_line.Draw("same")
    simc_peak_line = make_peak_line(simc_fit["mean"], simc_hist.GetMaximum(), kRed)
    simc_peak_line.Draw("same")
    simc_text = make_text_box(
        [
            f"SIMC peak = {simc_fit['mean']:.6f}",
            f"Fit err = {simc_fit['mean_err']:.6f}",
        ]
    )
    simc_text.Draw("same")
    canvas.Print(output_pdf)

    data_hist = data_fit["hist"].Clone(f"{data_fit['hist'].GetName()}_plot")
    data_hist.SetTitle(f"Data M_{particle_type[0].upper()} Fit")
    data_hist.SetLineColor(kBlack)
    data_hist.SetLineWidth(2)
    data_hist.GetXaxis().SetRangeUser(plot_xmin, plot_xmax)
    data_hist.Draw("hist")
    data_fit_line = data_hist.GetFunction("gaus")
    if data_fit_line:
        data_fit_line.Draw("same")
    data_peak_line = make_peak_line(data_fit["mean"], data_hist.GetMaximum(), kBlue)
    data_peak_line.Draw("same")
    simc_reference_line = make_peak_line(simc_fit["mean"], data_hist.GetMaximum(), kGreen + 2)
    simc_reference_line.Draw("same")
    data_text = make_text_box(
        [
            f"Data peak = {data_fit['mean']:.6f}",
            f"SIMC peak = {simc_fit['mean']:.6f}",
            f"MM shift = {shift:+.6f}",
        ]
    )
    data_text.Draw("same")
    canvas.Print(output_pdf)

    data_overlay = data_hist.Clone(f"{data_hist.GetName()}_overlay")
    shifted_overlay = shifted_hist.Clone(f"{shifted_hist.GetName()}_overlay")
    simc_overlay = simc_hist.Clone(f"{simc_hist.GetName()}_overlay")

    normalize_hist(data_overlay)
    normalize_hist(shifted_overlay)
    normalize_hist(simc_overlay)

    data_overlay.SetTitle(f"M_{particle_type[0].upper()} Raw / Shifted / SIMC")
    data_overlay.SetLineColor(kBlack)
    shifted_overlay.SetLineColor(kBlue)
    simc_overlay.SetLineColor(kRed)
    data_overlay.SetLineWidth(2)
    shifted_overlay.SetLineWidth(2)
    simc_overlay.SetLineWidth(2)

    ymax = max(
        data_overlay.GetMaximum(),
        shifted_overlay.GetMaximum(),
        simc_overlay.GetMaximum(),
    ) * 1.15
    data_overlay.SetMaximum(ymax)
    data_overlay.GetYaxis().SetTitle("Normalized counts")
    data_overlay.Draw("hist")
    shifted_overlay.Draw("hist same")
    simc_overlay.Draw("hist same")

    overlay_peak_line = make_peak_line(simc_fit["mean"], ymax, kGreen + 2)
    overlay_peak_line.Draw("same")

    legend = TLegend(0.58, 0.55, 0.89, 0.68)
    legend.SetBorderSize(1)
    legend.SetFillStyle(0)
    legend.AddEntry(data_overlay, "Data MM", "l")
    legend.AddEntry(shifted_overlay, "Data MM_shift", "l")
    legend.AddEntry(simc_overlay, "SIMC", "l")
    legend.Draw("same")

    overlay_text = make_text_box(
        [
            f"Data peak = {data_fit['mean']:.6f}",
            f"SIMC peak = {simc_fit['mean']:.6f}",
            f"Applied shift = {shift:+.6f}",
        ]
    )
    overlay_text.Draw("same")
    canvas.Print(output_pdf)

    canvas.Print(f"{output_pdf}]")
    canvas.Close()
    gStyle.SetOptStat(old_opt_stat)


def shift_experimental_files_to_simc_peak(
    particle_type,
    simc_filename,
    data_filename,
    dummy_filename=None,
    plot_xmin=None,
    plot_xmax=None,
):
    ROOT.gROOT.SetBatch(True)

    reference_tree = f"Cut_{particle_type.capitalize()}_Events_prompt_noRF"
    peak_name, _, _ = get_peak_window(particle_type)

    if plot_xmin is None:
        plot_xmin = FIT_HIST_XMIN
    if plot_xmax is None:
        plot_xmax = FIT_HIST_XMAX

    simc_fit = fit_tree_peak(
        simc_filename,
        "h10",
        "missmass",
        particle_type,
        hist_xmin=plot_xmin,
        hist_xmax=plot_xmax,
    )
    data_fit = fit_tree_peak(
        data_filename,
        reference_tree,
        "MM",
        particle_type,
        hist_xmin=plot_xmin,
        hist_xmax=plot_xmax,
    )
    shift = simc_fit["mean"] - data_fit["mean"]

    print(f"\nAligning experimental {peak_name} peak to SIMC...")
    print(f"SIMC peak ({simc_filename}) = {simc_fit['mean']:.6f} +/- {simc_fit['mean_err']:.6f}")
    print(f"Data peak ({data_filename}) = {data_fit['mean']:.6f} +/- {data_fit['mean_err']:.6f}")
    print(f"Applying MM shift = {shift:+.6f}")

    tree_names = get_data_tree_names(particle_type)
    add_shift_branch_to_file(data_filename, tree_names, "MM", shift)

    if dummy_filename:
        if os.path.exists(dummy_filename):
            add_shift_branch_to_file(dummy_filename, tree_names, "MM", shift)
        else:
            print(f"Dummy file '{dummy_filename}' not found. Skipping dummy MM shift.")

    shifted_hist = build_histogram(
        data_filename,
        reference_tree,
        "MM_shift",
        f"hist_shifted_{abs(hash((data_filename, reference_tree))) & 0xFFFFFFFF}",
        "MM_shift",
        hist_xmin=plot_xmin,
        hist_xmax=plot_xmax,
    )
    plot_filename = make_plot_filename(data_filename)
    write_shift_plots(
        particle_type,
        simc_fit,
        data_fit,
        shifted_hist,
        shift,
        plot_filename,
        plot_xmin=plot_xmin,
        plot_xmax=plot_xmax,
    )

    return {
        "simc_peak": simc_fit["mean"],
        "simc_peak_err": simc_fit["mean_err"],
        "data_peak": data_fit["mean"],
        "data_peak_err": data_fit["mean_err"],
        "shift": shift,
        "plot_filename": plot_filename,
    }


def main():
    if len(sys.argv) not in (4, 5, 6, 7):
        print(
            "Usage: shift_MM.py <particle_type> <simc_root> <data_root> [dummy_root] [mm_min mm_max]"
        )
        sys.exit(1)

    particle_type = sys.argv[1]
    simc_filename = sys.argv[2]
    data_filename = sys.argv[3]
    dummy_filename = None
    plot_xmin = None
    plot_xmax = None

    range_arg_start = 4
    if len(sys.argv) in (5, 7):
        dummy_filename = sys.argv[4]
        range_arg_start = 5

    if len(sys.argv) in (6, 7):
        plot_xmin = float(sys.argv[range_arg_start])
        plot_xmax = float(sys.argv[range_arg_start + 1])

    shift_experimental_files_to_simc_peak(
        particle_type,
        simc_filename,
        data_filename,
        dummy_filename,
        plot_xmin=plot_xmin,
        plot_xmax=plot_xmax,
    )


if __name__ == "__main__":
    main()
