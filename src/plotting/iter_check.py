#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-04 13:56:23 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import uproot as up
import numpy as np
import root_numpy as rnp
import ROOT
import scipy
import scipy.integrate as integrate
from scipy.integrate import quad
import matplotlib.pyplot as plt
from collections import defaultdict
import sys, math, os, subprocess
from array import array
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine, TLatex
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_Prod")

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
################################################################################################################################################

def plot_iteration(histlist, phisetlist, inpDict):

    kinematics = inpDict["kinematics"] 
    W = inpDict["W"] 
    Q2 = inpDict["Q2"] 
    EPSVAL = inpDict["EPSVAL"] 
    OutFilename = inpDict["OutFilename"] 
    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"] 
    NumtBins = inpDict["NumtBins"] 
    NumPhiBins = inpDict["NumPhiBins"] 
    ParticleType = inpDict["ParticleType"]
    closest_date = inpDict["closest_date"]
    formatted_date = inpDict["formatted_date"]

    # Create an empty list to store copied histograms
    histlist_copy = []

    # Check if the value of the dictionary is a TObject that can use .Clone()
    for hist in histlist:
        hist_copy = {}
        for key, val in hist.items():
            if hasattr(val, 'Clone') and callable(getattr(val, 'Clone')):
                # Clone the TObject if it has the 'Clone' method
                hist_copy[key] = val.Clone()
            else:
                # Otherwise, just copy the value
                hist_copy[key] = val
        # Append the copied histogram dictionary to the new list
        histlist_copy.append(hist_copy)

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################

    for hist in histlist_copy:

        print("\n\nReweighing {} histograms for comparison...".format(hist["phi_setting"].lower()))
        for key, val in hist.items():
            if hasattr(val, 'Clone') and callable(getattr(val, 'Clone')):
                hist["{}_OLD".format(key)] = val
                HistNumEvts = hist["{}_OLD".format(key)].GetNbinsX()
                for i, binIndex in enumerate(range(1, HistNumEvts + 1)):
                    # Progress bar
                    Misc.progressBar(i, HistNumEvts + 1,bar_length=25)
                    weight = hist["H_Weight_SIMC"].GetBinContent(binIndex)
                    hist["{}_OLD".format(key)].SetBinContent(binIndex, weight)

    CWeight = TCanvas()
    l_Weight = TLegend(0.1, 0.75, 0.35, 0.95)
    
    for i,hist in enumerate(histlist_copy):
        hist["H_Weight_SIMC"].SetLineColor(i+1)
        hist["H_Weight_SIMC"].Draw("same, E1")
        hist["H_iWeight_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_iWeight_SIMC"].Draw("same, E1")
        hist["H_iWeight_SIMC"].Draw("same, HIST")
        l_iWeight.AddEntry(hist["H_Weight_SIMC"],hist["phi_setting"]+"Data")
        l_iWeight.AddEntry(hist["H_iWeight_SIMC"],hist["phi_setting"]+"Simc")

    l_Weight.Draw()
    CWeight.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date))+'(')
                
    CQ2 = TCanvas()
    
    for i,hist in enumerate(histlist_copy):
        hist["H_Q2_SIMC_OLD"].SetLineColor(i+1)
        hist["H_Q2_SIMC_OLD"].Draw("same, E1")
        hist["H_Q2_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_Q2_SIMC"].Draw("same, E1")
        hist["H_Q2_SIMC"].Draw("same, HIST")

    CQ2.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))

    CW = TCanvas()
    
    for i,hist in enumerate(histlist_copy):
        hist["H_W_SIMC_OLD"].SetLineColor(i+1)
        hist["H_W_SIMC_OLD"].Draw("same, E1")
        hist["H_W_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_W_SIMC"].Draw("same, E1")
        hist["H_W_SIMC"].Draw("same, HIST")

    CW.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))

    Ct = TCanvas()
    
    for i,hist in enumerate(histlist_copy):
        hist["H_t_SIMC_OLD"].SetLineColor(i+1)
        hist["H_t_SIMC_OLD"].Draw("same, E1")
        hist["H_t_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_t_SIMC"].Draw("same, E1")
        hist["H_t_SIMC"].Draw("same, HIST")

    Ct.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))

    Cepsilon = TCanvas()
    
    for i,hist in enumerate(histlist_copy):
        hist["H_epsilon_SIMC_OLD"].SetLineColor(i+1)
        hist["H_epsilon_SIMC_OLD"].Draw("same, E1")
        hist["H_epsilon_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_epsilon_SIMC"].Draw("same, E1")
        hist["H_epsilon_SIMC"].Draw("same, HIST")

    Cepsilon.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))

    CMM = TCanvas()
    
    for i,hist in enumerate(histlist_copy):
        hist["H_MM_SIMC_OLD"].SetLineColor(i+1)
        hist["H_MM_SIMC_OLD"].Draw("same, E1")
        hist["H_MM_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_MM_SIMC"].Draw("same, E1")
        hist["H_MM_SIMC"].Draw("same, HIST")

    CMM.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date)))
    
    Cph_q = TCanvas()
    
    for i,hist in enumerate(histlist_copy):
        hist["H_ph_q_SIMC_OLD"].SetLineColor(i+1)
        hist["H_ph_q_SIMC_OLD"].Draw("same, E1")
        hist["H_ph_q_SIMC"].SetLineColor(i+(len(phisetlist)+1))
        hist["H_ph_q_SIMC"].Draw("same, E1")
        hist["H_ph_q_SIMC"].Draw("same, HIST")

    Cph_q.Print(outputpdf.replace("{}_FullAnalysis_".format(ParticleType),"{}_{}_simc_".format(ParticleType,formatted_date))+')')
