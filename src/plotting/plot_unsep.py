#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-09-20 16:03:42 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd
import root_pandas as rpd
import numpy as np
import ROOT
from ROOT import TCanvas, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TGaxis, TLine, TMath, TPaveText, TArc, TGraphPolar, TLatex, TH2Poly
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta
from functools import reduce
import re
import sys, math, os, subprocess
    
###############################################################################################################################################

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

def file_to_df(f_name, columns):
    '''
    Read in file and convert to dataframe with custom column names
    '''

    lineskip=False
    
    # Open the file for reading
    with open(f_name, 'r') as file:
        lines = file.readlines()
        if '1.000000\n' in lines:
            lineskip=True

    if lineskip:
        df = pd.read_csv(f_name, header=None, sep=' ', skiprows=1, skipfooter=1, engine='python')
    else:
        df = pd.read_csv(f_name, header=None, sep=' ')    
    df.columns = columns
    return df

################################################################################################################################################
'''
def fix_spacing(f_name):
    ''
    Fortran created files are bad with spacing. This fixes it.
    ''

    # Open the file for reading
    with open(f_name, 'r') as file:

        # Read the lines of the file and split based on whitespace
        lines = file.readlines()
        lines = [re.split(r'\s+', line.strip()) for line in lines]

        # Join the split lines with a single space
        lines = [' '.join(line) for line in lines]

        # Write the lines back to the file
        with open(f_name, 'w') as output:
            output.write('\n'.join(lines))

# Fix file spacing to work in pandas
fix_spacing(LTANAPATH+"/src/averages/avek.{}.dat".format(Q2.replace("p","")))
fix_spacing(LTANAPATH+"/src/xsects/x_unsep.{}_{}_{:.0f}".format(PID, Q2.replace("p",""), float(LOEPS)*100))
fix_spacing(LTANAPATH+"/src/xsects/x_unsep.{}_{}_{:.0f}".format(PID, Q2.replace("p",""), float(HIEPS)*100))

'''
################################################################################################################################################
# Read in files and convert to dataframes

def unsep_xsect(inpDict):

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
    POL = inpDict["POL"]
    
    if int(POL) == 1:
        pol_str = "pl"
    elif int(POL) == -1:
        pol_str = "mn"
    else:
        print("ERROR: Invalid polarity...must be +1 or -1")
        sys.exit(2)

    ################################################################################################################################################

    foutname = OUTPATH + "/" + ParticleType + "_unsep_" + OutFilename + ".root"
    fouttxt  = OUTPATH + "/" + ParticleType + "_unsep_" + OutFilename + ".txt"
    outputpdf  = OUTPATH + "/" + ParticleType + "_unsep_" + OutFilename + ".pdf"

    ################################################################################################################################################
    
    file_df_dict = {}

    setting_file = LTANAPATH+"/src/{}/list.settings".format(ParticleType)
    file_df_dict['setting_df'] = file_to_df(setting_file, ['POL', 'Q2', 'EPSVAL', 'thpq', 'TMIN', 'TMAX', 'NumtBins'])

    for i,row in file_df_dict['setting_df'].iterrows():
        if row['Q2'] == float(Q2.replace("p",".")):
            file_df_dict['beam_file'] = file_to_df(LTANAPATH+"/src/{}/beam/Eb_KLT.dat", ['ebeam', 'Q2', 'EPSVAL'])
            file_df_dict['avek_file'] = file_to_df(LTANAPATH+"/src/{}/averages/avek.{}.dat".format(Q2.replace("p","")) \
                                                   , ['W', 'dW', 'Q2', 'dQ2', 't', 'dt', 'th_pos', "tbin"])

            if row['EPSVAL'] == float(EPSVAL):
                if row['thpq'] < 0.0:
                    file_df_dict['aver_eps_{}'.format('right')] = file_to_df( \
                                                                                LTANAPATH+"/src/{}/averages/aver.{}_{}_{:.0f}.dat" \
                                                                                .format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100) \
                                                                                , ['ratio', 'dratio', 'phibin', 'tbin'])
                    file_df_dict['kindata_eps_{}'.format('right')] = file_to_df( \
                                                                                   LTANAPATH+"/src/{}/kindata/kindata.{}_{}_{:.0f}_{}.dat" \
                                                                                   .format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100, int(row['thpq']*1000)) \
                                                                                   , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
                if row['thpq'] > 0.0:
                    file_df_dict['aver_eps_{}'.format('left')] = file_to_df( \
                                                                               LTANAPATH+"/src/{}/averages/aver.{}_{}_{:.0f}.dat" \
                                                                               .format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100) \
                                                                               , ['ratio', 'dratio', 'phibin', 'tbin'])
                    file_df_dict['kindata_eps_{}'.format('left')] = file_to_df( \
                                                                                  LTANAPATH+"/src/{}/kindata/kindata.{}_{}_{:.0f}_+{}.dat" \
                                                                                  .format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100, int(row['thpq']*1000)) \
                                                                                  , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
                if row['thpq'] == 0.0:
                    file_df_dict['aver_eps_{}'.format('center')] = file_to_df( \
                                                                                 LTANAPATH+"/src/{}/averages/aver.{}_{}_{:.0f}.dat" \
                                                                                 .format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100) \
                                                                                 , ['ratio', 'dratio', 'phibin', 'tbin'])
                    file_df_dict['kindata_eps_{}'.format('center')] = file_to_df( \
                                                                                    LTANAPATH+"/src/{}/kindata/kindata.{}_{}_{:.0f}_+0000.dat" \
                                                                                    .format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100) \
                                                                                    , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
                file_df_dict['xsects_file_eps'] = file_to_df( \
                                                                LTANAPATH+"/src/{}/xsects/x_unsep.{}_{}_{:.0f}" \
                                                                .format(ParticleType, pol_str, Q2.replace("p",""), float(EPSVAL)*100) \
                                                                , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 'tm', 'W', 'Q2'])

    ################################################################################################################################################
    ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
    ################################################################################################################################################

    C_ratio_phi = TCanvas()
    C_ratio_phi.SetGrid()

    G_ratio_phi = ROOT.TGraphErrors()

    G_ratio_phi.SetTitle("eps = %s ; #phi_{bin}; Ratio" % LOEPS)

    phi_setting = ['left', 'center']

    for ps in phi_setting:
        for i in range(len(file_df_dict['aver_loeps_{}'.format(ps)]['ratio'].tolist())):
            G_ratio_phi.SetPoint(i, np.array(file_df_dict['aver_loeps_{}'.format(ps)]['phibin'].tolist())[i], np.array(file_df_dict['aver_loeps_{}'.format(ps)]['ratio'].tolist())[i])
            G_ratio_phi.SetPointError(i, 0, np.array(file_df_dict['aver_loeps_{}'.format(ps)]['dratio'].tolist())[i])
            G_ratio_phi.SetMarkerColor(i+1)

    G_ratio_phi.SetMarkerStyle(21)
    G_ratio_phi.SetMarkerSize(1)

    G_ratio_phi.Draw('AP')

    C_ratio_phi.Print(outputpdf + '(')

    C_ratio_t = TCanvas()
    C_ratio_t.SetGrid()

    G_ratio_t = ROOT.TGraphErrors()

    G_ratio_t.SetTitle("eps = %s ; t_{bin}; Ratio" % LOEPS)

    phi_setting = ['left', 'center']

    for ps in phi_setting:
        for i in range(len(file_df_dict['aver_loeps_{}'.format(ps)]['ratio'].tolist())):
            G_ratio_t.SetPoint(i, np.array(file_df_dict['aver_loeps_{}'.format(ps)]['tbin'].tolist())[i], np.array(file_df_dict['aver_loeps_{}'.format(ps)]['ratio'].tolist())[i])
            G_ratio_t.SetPointError(i, 0, np.array(file_df_dict['aver_loeps_{}'.format(ps)]['dratio'].tolist())[i])
            G_ratio_t.SetMarkerColor(i+1)

    G_ratio_t.SetMarkerStyle(21)
    G_ratio_t.SetMarkerSize(1)

    G_ratio_t.Draw('AP')

    C_ratio_t.Print(outputpdf)

    C_Q2_tbin = TCanvas()
    C_Q2_tbin.SetGrid()

    G_Q2_tbin = ROOT.TGraph()

    l_Q2_tbin = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_Q2_tbin.SetTitle("eps = %s ; #theta_{cm}; Q^{2} Average" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['Q2'].tolist())):
        G_Q2_tbin.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['Q2'].tolist())[i])
        l_Q2_tbin.AddEntry(G_Q2_tbin, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_Q2_tbin.SetMarkerColor(i+1)

    G_Q2_tbin.SetMarkerStyle(21)
    G_Q2_tbin.SetMarkerSize(1)

    G_Q2_tbin.Draw('AP')
    l_Q2_tbin.Draw()

    C_Q2_tbin.Print(outputpdf)

    C_W_tbin = TCanvas()
    C_W_tbin.SetGrid()

    G_W_tbin = ROOT.TGraph()

    l_W_tbin = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_W_tbin.SetTitle("eps = %s ; #theta_{cm}; W Average" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['W'].tolist())):
        G_W_tbin.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['W'].tolist())[i])
        l_W_tbin.AddEntry(G_W_tbin, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_W_tbin.SetMarkerColor(i+1)

    G_W_tbin.SetMarkerStyle(21)
    G_W_tbin.SetMarkerSize(1)

    G_W_tbin.Draw('AP')
    l_W_tbin.Draw()

    C_W_tbin.Print(outputpdf)

    C_t_tbin = TCanvas()
    C_t_tbin.SetGrid()

    G_t_tbin = ROOT.TGraph()

    l_t_tbin = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_t_tbin.SetTitle("eps = %s ; #theta_{cm}; t Average" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['tm'].tolist())):
        G_t_tbin.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i])
        l_t_tbin.AddEntry(G_t_tbin, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_t_tbin.SetMarkerColor(i+1)

    G_t_tbin.SetMarkerStyle(21)
    G_t_tbin.SetMarkerSize(1)

    G_t_tbin.Draw('AP')
    l_t_tbin.Draw()

    C_t_tbin.Print(outputpdf)

    C_Q2_phitbin = TCanvas()
    C_Q2_phitbin.SetGrid()

    G_Q2_phitbin = ROOT.TGraph()

    l_Q2_phitbin = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_Q2_phitbin.SetTitle("eps = %s ; #phi; Q^{2} Average" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['Q2'].tolist())):
        G_Q2_phitbin.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['Q2'].tolist())[i])
        l_Q2_phitbin.AddEntry(G_Q2_phitbin, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_Q2_phitbin.SetMarkerColor(i+1)

    G_Q2_phitbin.SetMarkerStyle(21)
    G_Q2_phitbin.SetMarkerSize(1)

    G_Q2_phitbin.Draw('AP')
    l_Q2_phitbin.Draw()

    C_Q2_phitbin.Print(outputpdf)

    C_W_phitbin = TCanvas()
    C_W_phitbin.SetGrid()

    G_W_phitbin = ROOT.TGraph()

    l_W_phitbin = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_W_phitbin.SetTitle("eps = %s ; #phi; W Average" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['W'].tolist())):
        G_W_phitbin.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['W'].tolist())[i])
        l_W_phitbin.AddEntry(G_W_phitbin, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_W_phitbin.SetMarkerColor(i+1)

    G_W_phitbin.SetMarkerStyle(21)
    G_W_phitbin.SetMarkerSize(1)

    G_W_phitbin.Draw('AP')
    l_W_phitbin.Draw()

    C_W_phitbin.Print(outputpdf)

    C_t_phitbin = TCanvas()
    C_t_phitbin.SetGrid()

    G_t_phitbin = ROOT.TGraph()

    l_t_phitbin = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_t_phitbin.SetTitle("eps = %s ; #phi; t Average" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['tm'].tolist())):
        G_t_phitbin.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i])
        l_t_phitbin.AddEntry(G_t_phitbin, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_t_phitbin.SetMarkerColor(i+1)

    G_t_phitbin.SetMarkerStyle(21)
    G_t_phitbin.SetMarkerSize(1)

    G_t_phitbin.Draw('AP')
    l_t_phitbin.Draw()

    C_t_phitbin.Print(outputpdf)

    C_xreal_th = TCanvas()
    C_xreal_th.SetGrid()

    G_xreal_th = ROOT.TGraphErrors()

    l_xreal_th = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_xreal_th.SetTitle("eps = %s ; #theta_{cm}; x_real" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['x_real'].tolist())):
        G_xreal_th.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['x_real'].tolist())[i])
        G_xreal_th.SetPointError(i, 0, np.array(file_df_dict['xsects_file_loeps']['dx_real'].tolist())[i])
        l_xreal_th.AddEntry(G_xreal_th, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_xreal_th.SetMarkerColor(i+1)

    G_xreal_th.SetMarkerStyle(21)
    G_xreal_th.SetMarkerSize(1)

    G_xreal_th.Draw('AP')
    l_xreal_th.Draw()

    C_xreal_th.Print(outputpdf)

    C_xmod_th = TCanvas()
    C_xmod_th.SetGrid()

    G_xmod_th = ROOT.TGraph()

    l_xmod_th = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_xmod_th.SetTitle("eps = %s ; #theta_{cm}; x_mod" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['x_mod'].tolist())):
        G_xmod_th.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['th_cm'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['x_mod'].tolist())[i])
        l_xmod_th.AddEntry(G_xmod_th, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_xmod_th.SetMarkerColor(i+1)

    G_xmod_th.SetMarkerStyle(21)
    G_xmod_th.SetMarkerSize(1)

    G_xmod_th.Draw('AP')
    l_xmod_th.Draw()

    C_xmod_th.Print(outputpdf)

    C_xreal_phi = TCanvas()
    C_xreal_phi.SetGrid()

    G_xreal_phi = ROOT.TGraphErrors()

    l_xreal_phi = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_xreal_phi.SetTitle("eps = %s ; #phi; x_real" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['x_real'].tolist())):
        G_xreal_phi.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['x_real'].tolist())[i])
        G_xreal_phi.SetPointError(i, 0, np.array(file_df_dict['xsects_file_loeps']['dx_real'].tolist())[i])
        l_xreal_phi.AddEntry(G_xreal_phi, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_xreal_phi.SetMarkerColor(i+1)

    G_xreal_phi.SetMarkerStyle(21)
    G_xreal_phi.SetMarkerSize(1)

    G_xreal_phi.Draw('AP')
    l_xreal_phi.Draw()

    C_xreal_phi.Print(outputpdf)

    C_xmod_phi = TCanvas()
    C_xmod_phi.SetGrid()

    G_xmod_phi = ROOT.TGraph()

    l_xmod_phi = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_xmod_phi.SetTitle("eps = %s ; #phi; x_mod" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['x_mod'].tolist())):
        G_xmod_phi.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['phi'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['x_mod'].tolist())[i])
        l_xmod_phi.AddEntry(G_xmod_phi, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_xmod_phi.SetMarkerColor(i+1)

    G_xmod_phi.SetMarkerStyle(21)
    G_xmod_phi.SetMarkerSize(1)

    G_xmod_phi.Draw('AP')
    l_xmod_phi.Draw()

    C_xmod_phi.Print(outputpdf)

    C_xreal_t = TCanvas()
    C_xreal_t.SetGrid()

    G_xreal_t = ROOT.TGraphErrors()

    l_xreal_t = ROOT.TLegend(0.8,0.8,0.95,0.95)

    G_xreal_t.SetTitle("eps = %s ; t; x_real" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['x_real'].tolist())):
        G_xreal_t.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['x_real'].tolist())[i])
        G_xreal_t.SetPointError(i, 0, np.array(file_df_dict['xsects_file_loeps']['dx_real'].tolist())[i])
        l_xreal_t.AddEntry(G_xreal_t, "t = {:.4f}".format(np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i]))
        G_xreal_t.SetMarkerColor(i+1)

    G_xreal_t.SetMarkerStyle(21)
    G_xreal_t.SetMarkerSize(1)

    G_xreal_t.Draw('AP')
    l_xreal_t.Draw()

    C_xreal_t.Print(outputpdf)

    C_xmod_t = TCanvas()
    C_xmod_t.SetGrid()

    G_xmod_t = ROOT.TGraph()

    G_xmod_t.SetTitle("eps = %s ; t; x_mod" % LOEPS)

    for i in range(len(file_df_dict['xsects_file_loeps']['x_mod'].tolist())):
        G_xmod_t.SetPoint(i, np.array(file_df_dict['xsects_file_loeps']['tm'].tolist())[i], np.array(file_df_dict['xsects_file_loeps']['x_mod'].tolist())[i])
        G_xmod_t.SetMarkerColor(i+1)

    G_xmod_t.SetMarkerStyle(21)
    G_xmod_t.SetMarkerSize(1)

    G_xmod_t.Draw('AP')

    C_xmod_t.Print(outputpdf + ')')
