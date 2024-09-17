#! /usr/bin/python
###########################################################################################################################
# Created - 26/July/22, Author - Jacob Murphy
# Based on script created - 20/July/21, Author - Muhammad Junaid (mjo147@uregina.ca), University of Regina, Canada (Copyright (c) junaid) #
###########################################################################################################################
# Python version of the pion plotting script. Now utilises uproot to select event of each type and writes them to a root file.
# Python should allow for easier reading of databases storing diferent variables.
# This version of script is for shift workers at JLab
# To run this script, execute: python3 scriptname runnumber

###################################################################################################################################################

# Import relevant packages
import ROOT
import numpy as np
import sys, math, os, subprocess
import array
import re # Regexp package - for string manipulation
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine, TCutG
from ROOT import TExec
from ROOT import kBlack, kBlue, kRed
from typing import Dict, Union, List, Tuple
from array import array
import pandas as pd
import glob

###############################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root

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
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
#################################################################################################################################################

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

def find_input_file(particle_type: str, filename_override: str, phi_setting: str, epsilon: str) -> Union[str, bool]:
    """Find the input file for a given epsilon setting."""
    file_pattern = f"{OUTPATH}/*{particle_type}*{filename_override}*{phi_setting}*.root"
    matching_files = glob.glob(file_pattern)
    
    if not matching_files:
        return False
    
    return min(matching_files, key=len)


def create_histograms(q2_min: float, q2_max: float, w_min: float, w_max: float, t_min: float, t_max: float) -> Dict[str, Union[TH2D, TH1D]]:
    """Create and return necessary histograms."""
    histograms = {
        "Q2vsW_cut": TH2D("Q2vsW_cut", "Q2 vs W Distribution", 400, q2_min, q2_max, 400, w_min, w_max),
        "Q2vsW_high_cut": TH2D("Q2vsW_high_cut", "High Epsilon Q2 vs W Distribution", 400, q2_min, q2_max, 400, w_min, w_max),
        "Q2vsW_mid_cut": TH2D("Q2vsW_mid_cut", "Mid Epsilon Q2 vs W Distribution", 400, q2_min, q2_max, 400, w_min, w_max),
        "Q2vsW_low_cut": TH2D("Q2vsW_low_cut", "Low Epsilon Q2 vs W Distribution", 400, q2_min, q2_max, 400, w_min, w_max),
        "Q2vsW_lolo_cut": TH2D("Q2vsW_lolo_cut", "Low Epsilon Q2 vs W Distribution (Diamond Cut)", 400, q2_min, q2_max, 400, w_min, w_max),
        "Q2vsW_milo_cut": TH2D("Q2vsW_milo_cut", "Mid Epsilon Q2 vs W Distribution (Diamond Cut)", 400, q2_min, q2_max, 400, w_min, w_max),
        "Q2vsW_hilo_cut": TH2D("Q2vsW_hilo_cut", "High Epsilon Q2 vs W Distribution (Diamond Cut)", 400, q2_min, q2_max, 400, w_min, w_max),
        "W_cut": TH1D("W_cut", "W Distribution", 400, w_min, w_max),
        "Q2_cut": TH1D("Q2_cut", "Q2 Distribution", 400, q2_min, q2_max),
        "t_cut": TH1D("t_cut", "-t Distribution", 400, t_min, t_max),
        "t_mi_cut": TH1D("t_mi_cut", "Mid Epsilon -t Distribution", 400, t_min, t_max),
    }
    return histograms

def fill_histograms(input_file: str, histograms: Dict[str, Union[TH2D, TH1D]], epsilon: str, particle_type: str) -> None:
    """Fill histograms with data from the input file."""
    try:
        with TFile.Open(input_file, "READ") as infile:
            tree = infile.Get(f"Cut_{particle_type.capitalize()}_Events_prompt_noRF")
            for event in tree:
                histograms["Q2vsW_cut"].Fill(event.Q2, event.W)
                histograms[f"Q2vsW_{epsilon}_cut"].Fill(event.Q2, event.W)
                if epsilon == "low":
                    histograms["W_cut"].Fill(event.W)
                    histograms["Q2_cut"].Fill(event.Q2)
    except Exception as e:
        print(f"Error processing file {input_file}: {str(e)}")

def fit_diamond(histogram: TH2D, q2_val: float, q2_min: float, q2_max: float, w_min: float, w_max: float) -> Tuple[List[float], List[float], List[float], List[float], List[float], List[float]]:
    """Fit the diamond shape to the histogram."""
    min_q = histogram.FindFirstBinAbove(0)
    max_q = histogram.FindLastBinAbove(0)
    fit_range = max(1, int((max_q - min_q) / 8))  # Ensure fit_range is at least 1
    fit_l = max(1, histogram.FindBin(q2_val) - fit_range * 2)
    fit_r = min(histogram.GetNbinsX(), histogram.FindBin(q2_val) + fit_range)

    lol, hil, lor, hir = [], [], [], []
    xvl, xvr = [], []

    for b in range(fit_range):
        y_l = histogram.ProjectionY("y", b + fit_l, b + fit_l + 1)
        y_r = histogram.ProjectionY("y", b + fit_r, b + fit_r + 1)

        if y_l.GetEntries() > 0 and y_r.GetEntries() > 0:
            min_yl = y_l.FindFirstBinAbove(0) / 400 * (w_max - w_min) + w_min
            max_yl = y_l.FindLastBinAbove(0) / 400 * (w_max - w_min) + w_min
            min_yr = y_r.FindFirstBinAbove(0) / 400 * (w_max - w_min) + w_min
            max_yr = y_r.FindLastBinAbove(0) / 400 * (w_max - w_min) + w_min

            lol.append(min_yl)
            hil.append(max_yl)
            lor.append(min_yr)
            hir.append(max_yr)

            xl = (b + fit_l) / 400 * (q2_max - q2_min) + q2_min
            xr = (b + fit_r) / 400 * (q2_max - q2_min) + q2_min
            xvl.append(xl)
            xvr.append(xr)

    print(f"Debug - fit_diamond: lol: {len(lol)}, hil: {len(hil)}, lor: {len(lor)}, hir: {len(hir)}, xvl: {len(xvl)}, xvr: {len(xvr)}")

    return lol, hil, lor, hir, xvl, xvr

def calculate_fit_parameters(xvl: List[float], xvr: List[float], lol: List[float], hil: List[float], lor: List[float], hir: List[float]) -> Dict[str, float]:
    """Calculate fit parameters for the diamond shape."""
    if not all([xvl, xvr, lol, hil, lor, hir]):
        print("Error: One or more input lists are empty")
        return {"a1": 0, "b1": 0, "a2": 0, "b2": 0, "a3": 0, "b3": 0, "a4": 0, "b4": 0}

    try:
        a1, b1 = np.polyfit(xvl, lol, 1)
        a2, b2 = np.polyfit(xvl, hil, 1)
        a3, b3 = np.polyfit(xvr, lor, 1)
        a4, b4 = np.polyfit(xvr, hir, 1)
    except Exception as e:
        print(f"Error in calculate_fit_parameters: {str(e)}")
        print(f"Debug - xvl: {xvl}")
        print(f"Debug - lol: {lol}")
        print(f"Debug - hil: {hil}")
        print(f"Debug - xvr: {xvr}")
        print(f"Debug - lor: {lor}")
        print(f"Debug - hir: {hir}")
        return {"a1": 0, "b1": 0, "a2": 0, "b2": 0, "a3": 0, "b3": 0, "a4": 0, "b4": 0}

    return {"a1": a1, "b1": b1, "a2": a2, "b2": b2, "a3": a3, "b3": b3, "a4": a4, "b4": b4}

def DiamondPlot(particle_type: str, q2_val: float, q2_min: float, q2_max: float, 
                 w_val: float, w_min: float, w_max: float, phi_setting: str, 
                 t_min: float, t_max: float, inp_dict: Dict) -> Dict[str, float]:
    """Main function for creating diamond plots."""
    
    filename_override = f"Q{str(q2_val).replace('.','p')}W{str(w_val).replace('.','p')}"
    output_file = f"{OUTPATH}/{phi_setting}_{particle_type}_diamond_{filename_override}.pdf"
    
    epsilon_settings = ["high", "mid", "low"]
    input_files = {eps: find_input_file(particle_type, filename_override, phi_setting, eps) 
                   for eps in epsilon_settings}
    
    histograms = create_histograms(q2_min, q2_max, w_min, w_max, t_min, t_max)
    
    for epsilon, input_file in input_files.items():
        if input_file:
            fill_histograms(input_file, histograms, epsilon, particle_type)
    
    if input_files["low"]:
        lol, hil, lor, hir, xvl, xvr = fit_diamond(histograms["Q2vsW_low_cut"], q2_val, q2_min, q2_max, w_min, w_max)
        fit_params = calculate_fit_parameters(xvl, xvr, lol, hil, lor, hir)
    else:
        fit_params = inp_dict

    for epsilon, input_file in input_files.items():
        if input_file:
            apply_diamond_cut(input_file, histograms, particle_type, fit_params, epsilon, t_max)

    draw_histograms(histograms, q2_min, q2_max, w_min, w_max, fit_params, output_file)
    
    return fit_params
