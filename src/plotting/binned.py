#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-21 10:51:14 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
# Import relevant packages
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
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine
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

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import match_to_bin

################################################################################################################################################

def plot_binned(t_bins, phi_bins, histlist, phisetlist, inpDict, yieldDict, ratioDict, aveDict):

    kinematics = inpDict["kinematics"] 
    W = inpDict["W"] 
    Q2 = inpDict["Q2"] 
    EPSVAL = inpDict["EPSVAL"] 
    InDATAFilename = inpDict["InDATAFilename"] 
    InDUMMYFilename = inpDict["InDUMMYFilename"] 
    OutFilename = inpDict["OutFilename"] 
    tmin = inpDict["tmin"] 
    tmax = inpDict["tmax"] 
    NumtBins = inpDict["NumtBins"] 
    NumPhiBins = inpDict["NumPhiBins"] 
    runNumRight = inpDict["runNumRight"] 
    runNumLeft = inpDict["runNumLeft"] 
    runNumCenter = inpDict["runNumCenter"]
    data_charge_right = inpDict["data_charge_right"] 
    data_charge_left = inpDict["data_charge_left"] 
    data_charge_center = inpDict["data_charge_center"] 
    dummy_charge_right = inpDict["dummy_charge_right"] 
    dummy_charge_left = inpDict["dummy_charge_left"] 
    dummy_charge_center = inpDict["dummy_charge_center"] 
    InData_efficiency_right = inpDict["InData_efficiency_right"] 
    InData_efficiency_left = inpDict["InData_efficiency_left"] 
    InData_efficiency_center = inpDict["InData_efficiency_center"] 
    efficiency_table = inpDict["efficiency_table"] 
    ParticleType = inpDict["ParticleType"]

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_" + OutFilename + ".pdf"

    ################################################################################################################################################


    # Find bin centers
    t_bin_centers = (t_bins[:-1] + t_bins[1:]) / 2
    phi_bin_centers = (phi_bins[:-1] + phi_bins[1:]) / 2    

    histbinDict = {}

    # t/phi binned histograms
    histbinDict["H_phibins_DATA"] = TH1D("H_phibins_DATA", "Phi Bins", NumtBins*NumPhiBins, 0, 360.0)
    histbinDict["H_tbins_DATA"] = TH1D("H_tbins_DATA", "t Bins", NumtBins*NumPhiBins, tmin, tmax)

    # Yield histograms
    for phiset in phisetlist:
        histbinDict["H_yield_DATA_{}".format(phiset)] = TH1D("H_yield_DATA_{}".format(phiset), "{} Data Yield".format(phiset), NumtBins*NumPhiBins, 0, 100.0)
        histbinDict["H_yield_SIMC_{}".format(phiset)] = TH1D("H_yield_SIMC_{}".format(phiset), "{} Simc Yield".format(phiset), NumtBins*NumPhiBins, 0, 100.0)

    # Ratio histogram
    for phiset in phisetlist:
        histbinDict["H_ratio_{}".format(phiset)] = TH1D("H_ratio_{}".format(phiset), "{} Ratio".format(phiset), NumtBins*NumPhiBins, -100.0, 100.0)

    # Loop over each tuple key in the dictionary
    for phiset in phisetlist:
        for k, data_key_tuple in enumerate(aveDict["binned_DATA"][phiset]['t']):
            # Access the nested dictionary using the tuple key
            data_nested_dict = aveDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = TH1D("H_Q2_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Data Q2 (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["Q2min"], inpDict["Q2max"])
            histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = TH1D("H_W_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Data W (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["Wmin"], inpDict["Wmax"])
            histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = TH1D("H_t_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Data t (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["tmin"], inpDict["tmax"])

    # Loop over each tuple key in the dictionary
    for phiset in phisetlist:
        for k, simc_key_tuple in enumerate(aveDict["binned_SIMC"][phiset]['t']):
            # Access the nested dictionary using the tuple key
            simc_nested_dict = aveDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin
            histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = TH1D("H_Q2_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Simc Q2 (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["Q2min"], inpDict["Q2max"])
            histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = TH1D("H_W_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Simc W (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["Wmin"], inpDict["Wmax"])
            histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = TH1D("H_t_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Simc t (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, inpDict["tmin"], inpDict["tmax"])

    # Loop over each tuple key in the dictionary
    for phiset in phisetlist:
        for k, data_key_tuple in enumerate(yieldDict["binned_DATA"][phiset]['yield']):
            # Access the nested dictionary using the tuple key
            data_nested_dict = yieldDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = TH1D("H_totevts_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Data Total Events (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, 0.0, 100.0)

    # Loop over each tuple key in the dictionary
    for phiset in phisetlist:
        for k, simc_key_tuple in enumerate(yieldDict["binned_SIMC"][phiset]['yield']):
            # Access the nested dictionary using the tuple key
            simc_nested_dict = yieldDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin
            histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))] = TH1D("H_totevts_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1)), "{} Simc Total Events (t bin {}, phi bin {})".format(phiset, str(i+1), str(j+1)), 100, 0.0, 100.0)

    C_Q2_tbin_DATA = TCanvas()
    C_Q2_tbin_DATA.Divide(NumtBins,NumPhiBins)
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
        for k, data_key_tuple in enumerate(data_key_tuples):
            # Access the nested dictionary using the tuple key
            data_nested_dict = aveDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(data_nested_dict["Q2"][data_key_tuple]["Q2_arr"]), data_nested_dict["Q2"][data_key_tuple]["Q2_ave"]))
            # Fill histogram
            for (itt,jtt), val in np.ndenumerate(data_nested_dict["Q2"][data_key_tuple]["Q2_arr"]):
                histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
            C_Q2_tbin_DATA.cd(k+1)
            histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
            histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
            histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
    C_Q2_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType))+'(')

    C_Q2_tbin_SIMC = TCanvas()
    C_Q2_tbin_SIMC.Divide(NumtBins,NumPhiBins)
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        simc_key_tuples = list(aveDict["binned_SIMC"][phiset]['t'])
        for k, simc_key_tuple in enumerate(simc_key_tuples):
            # Access the nested dictionary using the tuple key
            simc_nested_dict = aveDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin    
            # Fill histogram
            for (itt,jtt), val in np.ndenumerate(simc_nested_dict["Q2"][simc_key_tuple]["Q2_arr"]):
                histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
            C_Q2_tbin_SIMC.cd(k+1)
            histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
            histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
            histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
    C_Q2_tbin_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_W_tbin_DATA = TCanvas()
    C_W_tbin_DATA.Divide(NumtBins,NumPhiBins)
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
        for k, data_key_tuple in enumerate(data_key_tuples):
            # Access the nested dictionary using the tuple key
            data_nested_dict = aveDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin    
            # Fill histogram
            for (itt,jtt), val in np.ndenumerate(data_nested_dict["W"][data_key_tuple]["W_arr"]):
                histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
            C_W_tbin_DATA.cd(k+1)
            histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
            histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
            histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
    C_W_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_W_tbin_SIMC = TCanvas()
    C_W_tbin_SIMC.Divide(NumtBins,NumPhiBins)
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        simc_key_tuples = list(aveDict["binned_SIMC"][phiset]['t'])
        for k, simc_key_tuple in enumerate(simc_key_tuples):
            # Access the nested dictionary using the tuple key
            simc_nested_dict = aveDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin    
            # Fill histogram
            for (itt,jtt), val in np.ndenumerate(simc_nested_dict["W"][simc_key_tuple]["W_arr"]):
                histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
            C_W_tbin_SIMC.cd(k+1)
            histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
            histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
            histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
    C_W_tbin_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_t_tbin_DATA = TCanvas()
    C_t_tbin_DATA.Divide(NumtBins,NumPhiBins)
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
        for k, data_key_tuple in enumerate(data_key_tuples):
            # Access the nested dictionary using the tuple key
            data_nested_dict = aveDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            # Fill histogram
            for (itt,jtt), val in np.ndenumerate(data_nested_dict["t"][data_key_tuple]["t_arr"]):
                histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
            C_t_tbin_DATA.cd(k+1)
            histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
            histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
            histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
    C_t_tbin_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_t_tbin_SIMC = TCanvas()
    C_t_tbin_SIMC.Divide(NumtBins,NumPhiBins)
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        simc_key_tuples = list(aveDict["binned_SIMC"][phiset]['t'])
        for k, simc_key_tuple in enumerate(simc_key_tuples):
            # Access the nested dictionary using the tuple key
            simc_nested_dict = aveDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin    
            # Fill histogram
            for (itt,jtt), val in np.ndenumerate(simc_nested_dict["t"][simc_key_tuple]["t_arr"]):
                histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
            C_t_tbin_SIMC.cd(k+1)
            histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
            histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
            histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
    C_t_tbin_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_t_bins_DATA = TCanvas()
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
        for k, data_key_tuple in enumerate(data_key_tuples):
            # Access the nested dictionary using the tuple key
            data_nested_dict = aveDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            # Fill histogram
            for val in aveDict["binned_DATA"]["t_bins"]:
                histbinDict["H_tbins_DATA"].Fill(float(val))
        histbinDict["H_tbins_DATA"].SetLineColor(it+1)
        histbinDict["H_tbins_DATA"].Draw("same, hist")
    C_t_bins_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_phi_bins_DATA = TCanvas()
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        data_key_tuples = list(aveDict["binned_DATA"][phiset]['t'])
        for k, data_key_tuple in enumerate(data_key_tuples):
            # Access the nested dictionary using the tuple key
            data_nested_dict = aveDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            # Fill histogram
            for val in aveDict["binned_DATA"]["phi_bins"]:
                histbinDict["H_phibins_DATA"].Fill(float(val))
        histbinDict["H_phibins_DATA"].SetLineColor(it+1)            
        histbinDict["H_phibins_DATA"].Draw("same, hist")
    C_phi_bins_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_totevts_DATA = TCanvas()
    C_totevts_DATA.Divide(NumtBins,NumPhiBins)
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        data_key_tuples = list(yieldDict["binned_DATA"][phiset]['yield'])
        for k, data_key_tuple in enumerate(data_key_tuples):
            # Access the nested dictionary using the tuple key
            data_nested_dict = yieldDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(data_nested_dict["yield"][data_key_tuple]["yield_arr"]), data_nested_dict["yield"][data_key_tuple]["yield"]))
            # Fill histogram
            for (itt,jtt), val in np.ndenumerate(data_nested_dict["yield"][data_key_tuple]["yield_arr"]):
                histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
            C_totevts_DATA.cd(k+1)
            histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
            histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
            histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
    C_totevts_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_totevts_SIMC = TCanvas()
    C_totevts_SIMC.Divide(NumtBins,NumPhiBins)
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        simc_key_tuples = list(yieldDict["binned_SIMC"][phiset]['yield'])
        for k, simc_key_tuple in enumerate(simc_key_tuples):
            # Access the nested dictionary using the tuple key
            simc_nested_dict = yieldDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin
            #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(simc_nested_dict["yield"][simc_key_tuple]["yield_arr"]), simc_nested_dict["yield"][simc_key_tuple]["yield"]))
            # Fill histogram
            for (itt,jtt), val in np.ndenumerate(simc_nested_dict["yield"][simc_key_tuple]["yield_arr"]):
                histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Fill(val)
            C_totevts_SIMC.cd(k+1)
            histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].SetLineColor(it+1)
            histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, E1")
            histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset,str(i+1),str(j+1))].Draw("same, hist")
    C_totevts_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    for it,phiset in enumerate(phisetlist):

        C_yieldvsphi_data_plt = TCanvas()
        C_yieldvsphi_data_plt.SetGrid()
        C_yieldvsphi_data_plt.Divide(1,NumtBins)
        l_yieldvsphi_data_plt = TLegend(0.115,0.35,0.33,0.5)

        yield_data = []
        yield_simc = []
        phibins_data = []
        phibins_simc = []    
        data_key_tuples = list(yieldDict["binned_DATA"][phiset]['yield'])
        simc_key_tuples = list(yieldDict["binned_SIMC"][phiset]['yield'])
        for data_key_tuple,simc_key_tuple in zip(data_key_tuples,simc_key_tuples):
            tmp_yield_data = [[],[]]
            tmp_yield_simc = [[],[]]
            tmp_phibins_data = [[],[]]
            tmp_phibins_simc = [[],[]]
            # Access the nested dictionary using the tuple key
            data_nested_dict = yieldDict["binned_DATA"][phiset]
            simc_nested_dict = yieldDict["binned_SIMC"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            tmp_yield_data[0].append(yieldDict["binned_DATA"]["t_bins"][i])
            tmp_yield_data[1].append(data_nested_dict["yield"][data_key_tuple]["yield"])
            tmp_yield_simc[0].append(yieldDict["binned_SIMC"]["t_bins"][i])
            tmp_yield_simc[1].append(simc_nested_dict["yield"][simc_key_tuple]["yield"])
            tmp_phibins_data[0].append(yieldDict["binned_DATA"]["t_bins"][i])
            tmp_phibins_data[1].append(yieldDict["binned_DATA"]["phi_bins"][j])
            tmp_phibins_simc[0].append(yieldDict["binned_SIMC"]["t_bins"][i])
            tmp_phibins_simc[1].append(yieldDict["binned_SIMC"]["phi_bins"][j])
            yield_data.append(tmp_yield_data)
            yield_simc.append(tmp_yield_simc)
            phibins_data.append(tmp_phibins_data)
            phibins_simc.append(tmp_phibins_simc)

        # Match t-bins with list of yields
        yield_data = match_to_bin(yield_data)
        yield_simc = match_to_bin(yield_simc)
        phibins_data = match_to_bin(phibins_data)
        phibins_simc = match_to_bin(phibins_simc)

        multiDict = {}
        for i, val in enumerate(t_bin_centers):

            multiDict["G_yieldvsphi_plt_{}".format(i)] = TMultiGraph()

            G_yieldvsphi_data = TGraphErrors(len(yield_data[i][1]),phibins_data[i][1],yield_data[i][1],np.array([0]*len(phibins_data[i][1])),np.array([0]*len(yield_data[i][1])))
            G_yieldvsphi_simc = TGraphErrors(len(yield_simc[i][1]),phibins_simc[i][1],yield_simc[i][1],np.array([0]*len(phibins_simc[i][1])),np.array([0]*len(yield_simc[i][1])))

            G_yieldvsphi_data.SetMarkerStyle(21)
            G_yieldvsphi_data.SetMarkerSize(1)
            G_yieldvsphi_data.SetMarkerColor(1)
            multiDict["G_yieldvsphi_plt_{}".format(i)].Add(G_yieldvsphi_data)

            G_yieldvsphi_simc.SetMarkerStyle(21)
            G_yieldvsphi_simc.SetMarkerSize(1)
            G_yieldvsphi_simc.SetMarkerColor(2)
            multiDict["G_yieldvsphi_plt_{}".format(i)].Add(G_yieldvsphi_simc)

            C_yieldvsphi_data_plt.cd(i+1)

            multiDict["G_yieldvsphi_plt_{}".format(i)].Draw("AP")
            multiDict["G_yieldvsphi_plt_{}".format(i)].SetTitle("{} t = {};#phi; Yield".format(phiset,val))

            multiDict["G_yieldvsphi_plt_{}".format(i)].GetYaxis().SetTitleOffset(1.5)
            multiDict["G_yieldvsphi_plt_{}".format(i)].GetXaxis().SetTitleOffset(1.5)
            multiDict["G_yieldvsphi_plt_{}".format(i)].GetXaxis().SetLabelSize(0.04)
    
        l_yieldvsphi_data_plt.AddEntry(G_yieldvsphi_data,"Data")
        l_yieldvsphi_data_plt.AddEntry(G_yieldvsphi_simc,"Simc")
        l_yieldvsphi_data_plt.Draw()

        C_yieldvsphi_data_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    for it,phiset in enumerate(phisetlist):

        C_ratiovsphi_plt = TCanvas()
        C_ratiovsphi_plt.SetGrid()
        C_ratiovsphi_plt.Divide(1,NumtBins)

        ratio = []
        phibins = []
        key_tuples = list(ratioDict["binned"][phiset]['ratio'])
        for key_tuple in key_tuples:
            tmp_ratio = [[],[]]
            tmp_phibins = [[],[]]
            # Access the nested dictionary using the tuple key
            nested_dict = ratioDict["binned"][phiset]
            i = key_tuple[0] # t bin
            j = key_tuple[1] # phi bin
            tmp_ratio[0].append(ratioDict["binned"]["t_bins"][i])
            tmp_ratio[1].append(nested_dict["ratio"][key_tuple]["ratio"])
            tmp_phibins[0].append(ratioDict["binned"]["t_bins"][i])
            tmp_phibins[1].append(ratioDict["binned"]["phi_bins"][j])
            ratio.append(tmp_ratio)
            phibins.append(tmp_phibins)

        # Match t-bins with list of ratios
        ratio = match_to_bin(ratio)
        phibins = match_to_bin(phibins)

        multiDict = {}
        for i, val in enumerate(t_bin_centers):

            multiDict["G_ratiovsphi_plt_{}".format(i)] = TMultiGraph()

            G_ratiovsphi = TGraphErrors(len(ratio[i][1]),phibins[i][1],ratio[i][1],np.array([0]*len(phibins[i][1])),np.array([0]*len(ratio[i][1])))

            G_ratiovsphi.SetMarkerStyle(21)
            G_ratiovsphi.SetMarkerSize(1)
            G_ratiovsphi.SetMarkerColor(1)
            multiDict["G_ratiovsphi_plt_{}".format(i)].Add(G_ratiovsphi)

            C_ratiovsphi_plt.cd(i+1)

            multiDict["G_ratiovsphi_plt_{}".format(i)].Draw("AP")
            multiDict["G_ratiovsphi_plt_{}".format(i)].SetTitle("{} t = {};#phi; Ratio".format(phiset,val))

            multiDict["G_ratiovsphi_plt_{}".format(i)].GetYaxis().SetTitleOffset(1.5)
            multiDict["G_ratiovsphi_plt_{}".format(i)].GetXaxis().SetTitleOffset(1.5)
            multiDict["G_ratiovsphi_plt_{}".format(i)].GetXaxis().SetLabelSize(0.04)

        C_ratiovsphi_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_yield_DATA = TCanvas()
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        data_key_tuples = list(yieldDict["binned_DATA"][phiset]['yield'])
        for k, data_key_tuple in enumerate(data_key_tuples):
            # Access the nested dictionary using the tuple key
            data_nested_dict = yieldDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(data_nested_dict["yield"][data_key_tuple]["yield"]), data_nested_dict["yield"][data_key_tuple]["yield"]))
            # Fill histogram
            val = data_nested_dict["yield"][data_key_tuple]["yield"]
            histbinDict["H_yield_DATA_{}".format(phiset)].Fill(val)
        histbinDict["H_yield_DATA_{}".format(phiset)].SetLineColor(it+1)
        histbinDict["H_yield_DATA_{}".format(phiset)].Draw("same, E1")
        histbinDict["H_yield_DATA_{}".format(phiset)].Draw("same, hist")
    C_yield_DATA.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_yield_SIMC = TCanvas()
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        simc_key_tuples = list(yieldDict["binned_SIMC"][phiset]['yield'])
        for k, simc_key_tuple in enumerate(simc_key_tuples):
            # Access the nested dictionary using the tuple key
            simc_nested_dict = yieldDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin
            #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(simc_nested_dict["yield"][simc_key_tuple]["yield"]), simc_nested_dict["yield"][simc_key_tuple]["yield"]))
            # Fill histogram
            val = simc_nested_dict["yield"][simc_key_tuple]["yield"]
            histbinDict["H_yield_SIMC_{}".format(phiset)].Fill(val)
        histbinDict["H_yield_SIMC_{}".format(phiset)].SetLineColor(it+1)
        histbinDict["H_yield_SIMC_{}".format(phiset)].Draw("same, E1")
        histbinDict["H_yield_SIMC_{}".format(phiset)].Draw("same, hist")
    C_yield_SIMC.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_ratio = TCanvas()
    # Loop over each tuple key in the dictionary
    for it,phiset in enumerate(phisetlist):
        key_tuples = list(ratioDict["binned"][phiset]['ratio'])
        for k, key_tuple in enumerate(key_tuples):
            # Access the nested dictionary using the tuple key
            nested_dict = ratioDict["binned"][phiset]
            i = key_tuple[0] # t bin
            j = key_tuple[1] # phi bin
            #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(nested_dict["ratio"][key_tuple]["ratio"]), nested_dict["ratio"][key_tuple]["ratio"]))
            # Fill histogram
            val = nested_dict["ratio"][key_tuple]["ratio"]
            histbinDict["H_ratio_{}".format(phiset)].Fill(val)
        histbinDict["H_ratio_{}".format(phiset)].SetLineColor(it+1)
        histbinDict["H_ratio_{}".format(phiset)].Draw("same, E1")
        histbinDict["H_ratio_{}".format(phiset)].Draw("same, hist")
    C_ratio.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_yield_data_plt = TCanvas()
    G_yield_data_plt = TMultiGraph()
    l_yield_data_plt = TLegend(0.115,0.35,0.33,0.5)

    C_yield_data_plt.SetGrid()

    yield_data = np.array([])
    yield_simc = np.array([])
    setting = np.array([])
    for it,phiset in enumerate(phisetlist):
        data_key_tuples = list(yieldDict["binned_DATA"][phiset]['yield'])
        simc_key_tuples = list(yieldDict["binned_SIMC"][phiset]['yield'])
        for data_key_tuple,simc_key_tuple in zip(data_key_tuples,simc_key_tuples):
            # Access the nested dictionary using the tuple key
            data_nested_dict = yieldDict["binned_DATA"][phiset]        
            simc_nested_dict = yieldDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin
            #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(simc_nested_dict["yield"][simc_key_tuple]["yield"]), simc_nested_dict["yield"][simc_key_tuple]["yield"]))
            # Fill histogram
            yield_data = np.append(yield_data, [data_nested_dict["yield"][data_key_tuple]["yield"]])        
            yield_simc = np.append(yield_simc, [simc_nested_dict["yield"][simc_key_tuple]["yield"]])
            if phiset == "Center": setting = np.append(setting,0)
            elif phiset == "Left": setting = np.append(setting,1)
            else: setting = np.append(setting,2)

    G_yield_data = TGraphErrors(len(yield_data),setting,yield_data,np.array([0]*len(setting)),np.array([0]*len(yield_data)))
    G_yield_simc = TGraphErrors(len(yield_simc),setting,yield_simc,np.array([0]*len(setting)),np.array([0]*len(yield_simc)))

    G_yield_data.SetMarkerStyle(21)
    G_yield_data.SetMarkerSize(1)
    G_yield_data.SetMarkerColor(1)
    G_yield_data_plt.Add(G_yield_data)

    G_yield_simc.SetMarkerStyle(21)
    G_yield_simc.SetMarkerSize(1)
    G_yield_simc.SetMarkerColor(2)
    G_yield_data_plt.Add(G_yield_simc)

    G_yield_data_plt.Draw("AP")
    G_yield_data_plt.SetTitle(" ;Setting; Yield")

    for i,hist in enumerate(histlist):
        bin_ix = G_yield_data_plt.GetXaxis().FindBin(i)
        G_yield_data_plt.GetXaxis().SetBinLabel(bin_ix,hist["phi_setting"])

    G_yield_data_plt.GetYaxis().SetTitleOffset(1.5)
    G_yield_data_plt.GetXaxis().SetTitleOffset(1.5)
    G_yield_data_plt.GetXaxis().SetLabelSize(0.04)

    l_yield_data_plt.AddEntry(G_yield_data,"Data")
    l_yield_data_plt.AddEntry(G_yield_simc,"Simc")
    l_yield_data_plt.Draw()

    C_yield_data_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType)))

    C_ratio_plt = TCanvas()
    G_ratio_plt = TMultiGraph()

    C_ratio_plt.SetGrid()

    ratio_data = np.array([])
    setting = np.array([])
    for it,phiset in enumerate(phisetlist):
        key_tuples = list(ratioDict["binned"][phiset]['ratio'])
        for k, key_tuple in enumerate(key_tuples):
            # Access the nested dictionary using the tuple key
            nested_dict = ratioDict["binned"][phiset]
            i = key_tuple[0] # t bin
            j = key_tuple[1] # phi bin
            #print("~~~~~~~~~~~~~~~~~~~~~~",(k, i, j, len(nested_dict["ratio"][key_tuple]["ratio"]), nested_dict["ratio"][key_tuple]["ratio"]))
            ratio_data = np.append(ratio_data, [nested_dict["ratio"][key_tuple]["ratio"]])
            if phiset == "Center": setting = np.append(setting,0)
            elif phiset == "Left": setting = np.append(setting,1)
            else: setting = np.append(setting,2)

    G_ratio = TGraphErrors(len(ratio_data),setting,ratio_data,np.array([0]*len(setting)),np.array([0]*len(ratio_data)))

    G_ratio.SetMarkerStyle(21)
    G_ratio.SetMarkerSize(1)
    G_ratio.SetMarkerColor(1)
    G_ratio_plt.Add(G_ratio)

    G_ratio_plt.Draw("AP")

    G_ratio_plt.SetTitle(" ;Setting; Ratio")

    for i,hist in enumerate(histlist):
        bin_ix = G_ratio_plt.GetXaxis().FindBin(i)
        G_ratio_plt.GetXaxis().SetBinLabel(bin_ix,hist["phi_setting"])

    G_ratio_plt.GetYaxis().SetTitleOffset(1.5)
    G_ratio_plt.GetXaxis().SetTitleOffset(1.5)
    G_ratio_plt.GetXaxis().SetLabelSize(0.04)

    C_ratio_plt.Print(outputpdf.replace("{}_".format(ParticleType),"{}_binned_".format(ParticleType))+')')

    # Save histograms to list of dictionaries
    for hist in histlist:

        phiset = hist["phi_setting"]
        
        # t/phi binned histograms
        hist["H_phibins_DATA"] = histbinDict["H_phibins_DATA"]
        hist["H_tbins_DATA"] = histbinDict["H_tbins_DATA"]

        hist["H_yield_DATA"]     = histbinDict["H_yield_DATA_{}".format(phiset)]
        hist["H_yield_SIMC"]     = histbinDict["H_yield_SIMC_{}".format(phiset)]

        hist["H_ratio"]     = histbinDict["H_ratio_{}".format(phiset)]

        for k, data_key_tuple in enumerate(aveDict["binned_DATA"][phiset]['t']):
            # Access the nested dictionary using the tuple key
            data_nested_dict = aveDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            hist["H_Q2_tbin_DATA_{}_{}".format(str(i+1), str(j+1))] = histbinDict["H_Q2_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))]
            hist["H_W_tbin_DATA_{}_{}".format(str(i+1), str(j+1))] = histbinDict["H_W_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))]
            hist["H_t_tbin_DATA_{}_{}".format(str(i+1), str(j+1))] = histbinDict["H_t_tbin_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))]

        for k, simc_key_tuple in enumerate(aveDict["binned_SIMC"][phiset]['t']):
            # Access the nested dictionary using the tuple key
            simc_nested_dict = aveDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin
            hist["H_Q2_tbin_SIMC_{}_{}".format(str(i+1), str(j+1))] = histbinDict["H_Q2_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))]
            hist["H_W_tbin_SIMC_{}_{}".format(str(i+1), str(j+1))] = histbinDict["H_W_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))]
            hist["H_t_tbin_SIMC_{}_{}".format(str(i+1), str(j+1))] = histbinDict["H_t_tbin_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))]

        for k, data_key_tuple in enumerate(yieldDict["binned_DATA"][phiset]['yield']):
            # Access the nested dictionary using the tuple key
            data_nested_dict = yieldDict["binned_DATA"][phiset]
            i = data_key_tuple[0] # t bin
            j = data_key_tuple[1] # phi bin
            hist["H_totevts_DATA_{}_{}".format(str(i+1), str(j+1))] = histbinDict["H_totevts_DATA_{}_{}_{}".format(phiset, str(i+1), str(j+1))]

        for k, simc_key_tuple in enumerate(yieldDict["binned_SIMC"][phiset]['yield']):
            # Access the nested dictionary using the tuple key
            simc_nested_dict = yieldDict["binned_SIMC"][phiset]
            i = simc_key_tuple[0] # t bin
            j = simc_key_tuple[1] # phi bin
            hist["H_totevts_SIMC_{}_{}".format(str(i+1), str(j+1))] = histbinDict["H_totevts_SIMC_{}_{}_{}".format(phiset, str(i+1), str(j+1))]
