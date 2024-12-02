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
from array import array
import pandas as pd
from scipy.optimize import curve_fit
import glob

###############################################################################################################################################
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

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import open_root_file, apply_bin_threshold

################################################################################################################################################
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
#################################################################################################################################################

print("Running as %s on %s, hallc_replay_lt path assumed as %s" % (USER, HOST, REPLAYPATH))

def diamond_fit(Q2vsW_hist, Q2Val, fitrange=10):
    def line(x, a, b):
        return a * x + b

    if fitrange <= 1:
        print("Fit range too small, setting to default of 10.")
        fitrange = 10
    
    lol, hil, lor, hir = [], [], [], []
    xvl, xvr = [], []
    
    fitl = Q2vsW_hist.GetXaxis().FindBin(Q2Val) - int(fitrange*5)
    fitr = Q2vsW_hist.GetXaxis().FindBin(Q2Val) + int(fitrange)

    print("\nFinding diamond fits...")
    for b in range(fitrange*2):
        # Progress bar
        Misc.progressBar(b, fitrange*2-1, bar_length=25)
        
        # Left side
        proj_l = Q2vsW_hist.ProjectionY("y", b+fitl, b+fitl)
        minYl = proj_l.GetBinCenter(proj_l.FindFirstBinAbove(0))
        maxYl = proj_l.GetBinCenter(proj_l.FindLastBinAbove(0))
        lol.append(minYl)
        hil.append(maxYl)
        
        # Right side
        proj_r = Q2vsW_hist.ProjectionY("y", b+fitr, b+fitr+1)
        minYr = proj_r.GetBinCenter(proj_r.FindFirstBinAbove(0))
        maxYr = proj_r.GetBinCenter(proj_r.FindLastBinAbove(0))
        lor.append(minYr)
        hir.append(maxYr)
        
        xl = Q2vsW_hist.GetXaxis().GetBinCenter(b+fitl)
        xr = Q2vsW_hist.GetXaxis().GetBinCenter(b+fitr)
        xvl.append(xl)
        xvr.append(xr)

    # Fit all four sides
    popt_l_low, _ = curve_fit(line, xvl, lol)
    popt_l_high, _ = curve_fit(line, xvl, hil)
    popt_r_low, _ = curve_fit(line, xvr, lor)
    popt_r_high, _ = curve_fit(line, xvr, hir)

    return popt_l_low, popt_l_high, popt_r_low, popt_r_high

def DiamondPlot(ParticleType, Q2Val, Q2min, Q2max, WVal, Wmin, Wmax, phi_setting, tmin, tmax, inpDict):

    ##############
    # HARD CODED #
    ##############
    if Q2Val == 0.4:
        Q2Val = 0.385
        WVal = 2.20
        Qs = "0p4"
        Ws = "2p20"
    elif Q2Val == "2p1":
        Q2Val = 2.115
        WVal = 2.95
        Qs = "2p1"
        Ws = "2p95"
    else:
        Qs = str(Q2Val).replace('.','p')
        Ws = str(WVal).replace('.','p')
    ##############
    ##############
    ##############
    
    FilenameOverride = 'Q'+Qs+'W'+Ws
    
    Analysis_Distributions = OUTPATH+"/{}_{}_diamond_{}.pdf".format(phi_setting, ParticleType, FilenameOverride)

    nbins = 200
    
    lowe_input = False
    mide_input = False
    highe_input = False

    # Arbitrary lengths far longer than any file name
    lenh = 10000
    lenm = 10000
    lenl = 10000
    if(phi_setting == '0'): phi_setting = ""
    print("\n\nKinematics: ",FilenameOverride,"\nPhi Setting: ",phi_setting)
    for file in glob.glob(OUTPATH+'/*'+phi_setting+'*'+ParticleType+'*'+FilenameOverride+'*.root'):
	# Searches through OUTPUT recursively for files matching the wild card format, taking the shortest one
        # Shortest file assumed to be full analyisis as it will not have "part" or "week" or "dummy" labels
        #print(file)
        if "high" in file:
            if (len(file) < lenh):
                highe_input = file
                lenh = len(file)
        
        if "mid" in file:
            if (len(file) < lenm):
                mide_input = file
                lenm = len(file)

        if "low" in file:
            if (len(file) < lenl):
                lowe_input = file
                lenl = len(file)

    if (highe_input == False and mide_input == False and lowe_input == False):
        print("!!!!! ERROR !!!!!\n No valid file found! \n!!!!! ERROR !!!!!")
        sys.exit(1)

    ##############################################################################################################################################
    labelh = ""
    labelm = ""
    labell = ""
    if (highe_input !=False):
        #print("test high")
        labelh = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{High (Blue) Epsilon}; Q2; W"
        if (mide_input !=False):
            labelh = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{High (Blue) and Mid (Red) Epsilon}; Q2; W"
            #print("test high mid")
            if (lowe_input !=False):
                labelh = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{High (Blue), Mid (Red), and Low (Green) Epsilon}; Q2; W"
                #print("test high mid low")
        elif (lowe_input !=False):
            #print("test high low")
            labelh = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{High (Blue) and Low (Red) Epsilon}; Q2; W"
    elif (mide_input !=False):
        #print("test mid")
        labelm = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{Mid (Blue) Epsilon}; Q2; W"
        if (lowe_input !=False):
            #print("test mid low")
            labelm = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{Mid (Blue) and Low (Red) Epsilon}; Q2; W"
    elif (lowe_input !=False):
        #print("test low")
        labell = "#splitline{Q2 vs W Dist for Prompt Events (Prompt Cut)}{Low (Blue) Epsilon}; Q2; W"


    Title = ""
    Q2vsW_cut = TH2D("Q2vsW_cut", labelh, nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_mide_cut = TH2D("Q2vsW_mide_cut",labelm, nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_lowe_cut = TH2D("Q2vsW_lowe_cut", labell, nbins, Q2min, Q2max, nbins, Wmin, Wmax)    

    Q2vsW_hi_cut = TH2D("Q2vsW_high_cut", "High Epsilon Q2 vs W Dist for Prompt Events (Prompt Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_mi_cut = TH2D("Q2vsW_middle_cut","Mid Epsilon Q2 vs W Dist for Prompt Events (Prompt Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_lo_cut = TH2D("Q2vsW_low_cut", "Low Epsilon Q2 vs W Dist for Prompt Events (Prompt Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)

    W_cut = TH1D("W_cut", "High Epsilon W Dist for Prompt Events (Prompt Cut); W", nbins, Wmin, Wmax)
    Q2_cut = TH1D("Q2_cut", "High Epsilon Q2 Dist for Prompt Events (Prompt  Cut); Q2", nbins, Q2min, Q2max)
    t_cut = TH1D("t_cut", "High Epsilon -t Dist for Prompt Events (t-Range  Cut); -t", nbins, tmin, tmax)
    t_mi_cut = TH1D("t_mi_cut", "Mid Epsilon -t Dist for Prompt Events (t-Range  Cut); -t", nbins, tmin, tmax)

    Q2vsW_lolo_cut = TH2D("Q2vsW_low_lowcut", "Low Epsilon Q2 vs W Dist for Prompt Events (Diamond Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_hilo_cut = TH2D("Q2vsW_high_lowcut", "High Epsilon Q2 vs W Dist for Prompt Events (Diamond and t Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_milo_cut = TH2D("Q2vsW_mid_lowcut","Mid Epsilon Q2 vs W Dist for Prompt Events (Diamond and t Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)
    Q2vsW_himi_cut = TH2D("Q2vsW_high_midcut", "High Epsilon Q2 vs W Dist for Prompt Events (Mid-Diamond and t Cut); Q2; W", nbins, Q2min, Q2max, nbins, Wmin, Wmax)

    a1 = 0
    b1 = 0
    a2 = 0
    b2 = 0
    a3 = 0
    b3 = 0
    a4 = 0
    b4 = 0
    minQ = 0
    maxQ = 0
    fitl = 0
    fitr = 0
    fitrange = 0
    lol = []
    lor = []
    hil = []
    hir = []
    xvl = []
    xvr = []
    for k in range(0, 3):

	# Construct the name of the rootfile based upon the info we provided
        if (k==2):
            if (highe_input == False): 
                continue
            elif (highe_input != False): # Special condition, with 5th arg, use 5th arg as file name
                rootName = highe_input
                print("\n***** High Epsilon File Found! *****\n")
        if (k==1):
            if (mide_input == False): 
                continue
            elif (mide_input != False): # Special condition, with 5th arg, use 5th arg as file name
                rootName = mide_input
                print("\n*****Mid Epsilon File Found! *****\n")
        if (k==0):
            if (lowe_input == False): 
                continue
            elif (lowe_input != False): # Special condition, with 5th arg, use 5th arg as file name
                rootName = lowe_input
                print("\n***** Low Epsilon File Found! *****\n")
        print ("Attempting to process %s" %(rootName))

	###############################################################################################################################################
        ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen

	###############################################################################################################################################

	# Read stuff from the main event tree
        infile = open_root_file(rootName, "READ")

	# Assumes 2021 trees do not have Prompt MM cut, as some do not right now. *** NEED TO BE REPLAYED AGAIN WITH THIS BRANCH ***
        Cut_Events_all_noRF_tree = infile.Get("Cut_{}_Events_prompt_noRF".format(ParticleType.capitalize()))

	##############################################################################################################################################
        countB = 0
        countA = 0
        badfit = True
        if (k==2): # High
            # for event in Cut_Events_Prompt_tree:
            for event in Cut_Events_all_noRF_tree:
                Q2vsW_cut.Fill(event.Q2, event.W)
                Q2vsW_hi_cut.Fill(event.Q2, event.W)
        elif (k==1): # Mid
            # for event in Cut_Events_Prompt_tree:
            for event in Cut_Events_all_noRF_tree:
                Q2vsW_mide_cut.Fill(event.Q2, event.W)
                Q2vsW_mi_cut.Fill(event.Q2, event.W)
        elif (k==0): # Low
            # for event in Cut_Events_Prompt_tree:
            for event in Cut_Events_all_noRF_tree:
                Q2vsW_lowe_cut.Fill(event.Q2, event.W)
                Q2vsW_lo_cut.Fill(event.Q2, event.W)
                W_cut.Fill(event.W)
                Q2_cut.Fill(event.Q2)
                countB +=1

        # Apply the threshold to the histograms
        ##############
        # HARD CODED #
        ##############
        event_threshold = 5 # Q2=2.115,W=2.95 | Q2=3.0,W=2.32 | Q2=4.4,W=2.74
        #event_threshold = 15 # Q2=3.0,W=3.14
        ##############
        ##############
        ##############
        if (k==2): # High
            apply_bin_threshold(Q2vsW_cut, event_threshold)
            apply_bin_threshold(Q2vsW_hi_cut, event_threshold)
        elif (k==1): # Mid
            apply_bin_threshold(Q2vsW_mide_cut, event_threshold)
            apply_bin_threshold(Q2vsW_mi_cut, event_threshold)
        elif (k==0): # Low
            apply_bin_threshold(Q2vsW_lowe_cut, event_threshold)
            apply_bin_threshold(Q2vsW_lo_cut, event_threshold)
                
        #Does assume nbins bins for Q2 and W, centered at kinematic values
        minQ = Q2_cut.FindFirstBinAbove(0)
        maxQ = Q2_cut.FindLastBinAbove(0)
        fitrange = int((maxQ-minQ)/100)
        print("fitrange: ",fitrange)                
        if (k == 0):  # Low epsilon
            # Replace the existing diamond fitting code with:
            fit_results = diamond_fit(Q2vsW_lowe_cut, Q2Val, fitrange)
            a1, b1 = fit_results[0]
            a2, b2 = fit_results[1]
            a3, b3 = fit_results[2]
            a4, b4 = fit_results[3]

            for event in Cut_Events_all_noRF_tree:
                if (event.W > a1*event.Q2 + b1 and
                    event.W < a2*event.Q2 + b2 and
                    event.W > a3*event.Q2 + b3 and
                    event.W < a4*event.Q2 + b4):
                    Q2vsW_lolo_cut.Fill(event.Q2, event.W)
                    countA += 1
                    
        if (lowe_input != False and k>0):
            print("\n\n")
            if (k==2):
                for event in Cut_Events_all_noRF_tree:
                    if (event.W/event.Q2>a1+b1/event.Q2 and event.W/event.Q2<a2+b2/event.Q2 and event.W/event.Q2>a3+b3/event.Q2 and event.W/event.Q2<a4+b4/event.Q2):
                        if (tmax != False):
                            Q2vsW_hilo_cut.Fill(event.Q2, event.W)
                            t_cut.Fill(-event.MandelT)
                        else:
                            print("!!!!! Error! tmax not found! Skipping t-range cut !!!!!")
                            Q2vsW_hilo_cut.Fill(event.Q2, event.W)
            elif (k==1):
                for event in Cut_Events_all_noRF_tree:
                    if (event.W/event.Q2>a1+b1/event.Q2 and event.W/event.Q2<a2+b2/event.Q2 and event.W/event.Q2>a3+b3/event.Q2 and event.W/event.Q2<a4+b4/event.Q2):
                        if (tmax != False):
                            Q2vsW_milo_cut.Fill(event.Q2, event.W)
                            t_mi_cut.Fill(-event.MandelT)
                        else:
                            print("!!!!! Error! tmax not found! Skipping t-range cut !!!!!")
                            Q2vsW_milo_cut.Fill(event.Q2, event.W)        

        print("Histograms filled")

        infile.Close()

    if phi_setting == "Center":        
        paramDict = {

            "a1" : a1,
            "b1" : b1,
            "a2" : a2,
            "b2" : b2,
            "a3" : a3,
            "b3" : b3,
            "a4" : a4,
            "b4" : b4
        }

    else:

        paramDict = {}

    ##############################################################################################################################################
    c1_kin = TCanvas("c1_kin", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
    gStyle.SetTitleFontSize(0.03)
    gStyle.SetPalette(86)
    ex1 = TExec("ex1","gStyle->SetPalette(86)")
    ex2 = TExec("ex2","gStyle->SetPalette(75)")
    ex3 = TExec("ex3","gStyle->SetPalette(68)")
    gStyle.SetOptStat(0)
    pages = 2
    if (highe_input !=False):
        #print("test high")
        Q2vsW_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
        Q2vsW_cut.Draw("col")
        ex1.Draw()
        Q2vsW_cut.Draw("col same")
        if (mide_input !=False):
            #print("test high mid")
            Q2vsW_mide_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_mide_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
            ex2.Draw()
            Q2vsW_mide_cut.Draw("col same")
            pages = 3
            if (lowe_input !=False):
                #print("test high mid low")
                Q2vsW_lowe_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
                Q2vsW_lowe_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)                
                ex3.Draw()
                Q2vsW_lowe_cut.Draw("col same")
                pages = 6
        elif (lowe_input !=False):
            #print("test high low")
            Q2vsW_lowe_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_lowe_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
            ex2.Draw()
            Q2vsW_lowe_cut.Draw("col same")
            pages = 4
    elif (mide_input !=False):
	#print("test mid")
        Q2vsW_mide_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_mide_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_mide_cut.Draw("colz")
        ex1.Draw()
        Q2vsW_mide_cut.Draw("col same")
        if (lowe_input !=False):
            #print("test mid low")
            Q2vsW_lowe_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_lowe_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)            
            ex2.Draw()
            Q2vsW_lowe_cut.Draw("col same")
            pages = 4
    elif (lowe_input !=False):
        #print("test low")
        Q2vsW_lowe_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_lowe_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_lowe_cut.Draw("colz")
        ex1.Draw()
        Q2vsW_lowe_cut.Draw("col same")

    c1_kin.Print(Analysis_Distributions + "(")

    #############################################################################################################################

    end = ""
    endm = ""
    endc = ""
    endf = ""
    if (pages==2): end = ")"
    if (pages==3): endm = ")"
    if (pages==4): endc = ")"
    if (pages==6): endf = ")"

    gStyle.SetOptStat(1)
    gStyle.SetPalette(55)

    if (tmax != False):
        c1_kint = TCanvas("c1_kint", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        t_cut.Draw("colz")
        c1_kint.Print(Analysis_Distributions)
    if (highe_input != False):
        c1_kinh = TCanvas("c1_kinh", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        Q2vsW_hi_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_hi_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
        Q2vsW_hi_cut.Draw("colz")
        c1_kinh.Print(Analysis_Distributions+end)
    if (mide_input != False):
        c1_kinm = TCanvas("c1_kinm", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        Q2vsW_mi_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_mi_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_mi_cut.Draw("colz")
        c1_kinm.Print(Analysis_Distributions+end+endm)
    if (lowe_input != False):
        c1_kinl = TCanvas("c1_kinl", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        Q2vsW_lo_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_lo_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_lo_cut.Draw("colz")
        c1_kinl.Print(Analysis_Distributions)
        c1_kinll = TCanvas("c1_kinll", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
        Q2vsW_lolo_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
        Q2vsW_lolo_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)        
        Q2vsW_lolo_cut.Draw("colz")
        #    lol.clear()
        #   lor.clear()
        #  hil.clear()
        # hir.clear()
        # xvl.clear()
        # xvr.clear()
        #print (xvl[0], xvl[-1])
        line1 = TLine(Q2min, a1*Q2min+b1, Q2max, a1*Q2max+b1)   
        line2 = TLine(Q2min, a2*Q2min+b2, Q2max, a2*Q2max+b2) 
        line3 = TLine(Q2min, a3*Q2min+b3, Q2max, a3*Q2max+b3) 
        line4 = TLine(Q2min, a4*Q2min+b4, Q2max, a4*Q2max+b4)  
        line1.SetLineColor(1)
        line2.SetLineColor(2)
        line3.SetLineColor(3)
        line4.SetLineColor(4)  
        line1.SetLineWidth(5)
        line2.SetLineWidth(5)
        line3.SetLineWidth(5)
        line4.SetLineWidth(5)
        line1.Draw()
        line2.Draw()
        line3.Draw()
        line4.Draw()
        x1 = (nbins*0.25)*(Q2max-Q2min)+Q2min
        x2 = (nbins*0.75)*(Q2max-Q2min)+Q2min
        line1f = TLine(x1,a1*x1+b1,x2,a1*x2+b1)   
        line2f = TLine(x1,a2*x1+b2,x2,a2*x2+b2) 
        line3f = TLine(x1,a3*x1+b3,x2,a3*x2+b3) 
        line4f = TLine(x1,a4*x1+b4,x2,a4*x2+b4)  
        line1f.SetLineColor(1)
        line2f.SetLineColor(2)
        line3f.SetLineColor(3)
        line4f.SetLineColor(4)  
        line1f.SetLineWidth(2)
        line2f.SetLineWidth(2)
        line3f.SetLineWidth(2)
        line4f.SetLineWidth(2)
        line1f.Draw()
        line2f.Draw()
        line3f.Draw()
        line4f.Draw()
        c1_kinll.Print(Analysis_Distributions+end)

        if (mide_input != False):
            c1_kinml = TCanvas("c1_kinml", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
            Q2vsW_milo_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_milo_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)            
            Q2vsW_milo_cut.Draw("colz")
            c1_kinml.Print(Analysis_Distributions+endc)

        if (highe_input != False):
            c1_kinhl = TCanvas("c1_kinhl", "%s Kinematic Distributions" % ParticleType, 100, 0, 1000, 900)
            Q2vsW_hilo_cut.GetXaxis().SetRangeUser(Q2min-Q2min*0.1, Q2max+Q2max*0.1)
            Q2vsW_hilo_cut.GetYaxis().SetRangeUser(Wmin-Wmin*0.1, Wmax+Wmax*0.1)
            Q2vsW_hilo_cut.Draw("colz")
            c1_kinhl.Print(Analysis_Distributions+endc+endf)
	
            
    return paramDict
