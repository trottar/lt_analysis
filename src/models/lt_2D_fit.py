#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-05-07 02:50:07 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import ROOT
from ROOT import TGraphErrors, TF1, TF2, TGraph2DErrors, TCanvas
from ROOT import TString, TNtuple, TMinuit
from array import array
import math
import ctypes
import os, sys

ParticleType = sys.argv[1]
POL = float(sys.argv[2])

if POL > 0:
    polID = 'pl'
else:
    polID = 'mn'

Q2 = sys.argv[3]
W = sys.argv[4]

LOEPS = float(sys.argv[5])
HIEPS = float(sys.argv[6])

################################################################################################################################################
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

################################################################################################################################################

outputpdf  = "{}/{}_lt_fit_Q{}W{}.pdf".format(OUTPATH, ParticleType, Q2, W)

################################################################################################################################################
# Suppressing the terminal splash of Print()
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kError;")
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Constants
#pt_to_pt_systematic_error = 2.9 # Percent, just matching Bill's for now
pt_to_pt_systematic_error = 3.6 # In percent, matches PAC propsal projections (https://redmine.jlab.org/attachments/download/635/k12_proposal.pdf)
PI = math.pi

###############################################################################################################################################
# Import separated xsects models

from lt_active import LT_sep_x_lo_fun, LT_sep_x_lo_fun_unsep, LT_sep_x_hi_fun, LT_sep_x_hi_fun_unsep, set_val

###############################################################################################################################################

def single_setting(q2_set, w_set, fn_lo, fn_hi):

    # Set HIEPS for lt_active script
    set_val(LOEPS, HIEPS)

    eps_diff = HIEPS-LOEPS
    
    sig_L_g  = TGraphErrors()
    sig_T_g  = TGraphErrors()
    sig_LT_g = TGraphErrors()
    sig_TT_g = TGraphErrors()

    sig_lo = TGraphErrors()
    sig_hi = TGraphErrors()
    sig_diff_g = TGraphErrors()
    
    nlo = TNtuple("nlo", "nlo", "x/F:dx:x_mod:eps:theta:phi:t:w:Q2")
    nlo.ReadFile(fn_lo)
    
    nhi = TNtuple("nhi", "nhi", "x/F:dx:x_mod:eps:theta:phi:t:w:Q2")
    nhi.ReadFile(fn_hi)

    '''
    # Print values for nlo
    print("Values for nlo:")
    for evt in nlo:
        print("x: {}, dx: {}, x_mod: {}, eps: {}, theta: {}, phi: {}, t: {}, w: {}, Q2: {}".format(
            evt.x, evt.dx, evt.x_mod, evt.eps, evt.theta, evt.phi, evt.t, evt.w, evt.Q2
        ))

    # Print values for nhi
    print("\nValues for nhi:")
    for evt in nhi:
        print("x: {}, dx: {}, x_mod: {}, eps: {}, theta: {}, phi: {}, t: {}, w: {}, Q2: {}".format(
            evt.x, evt.dx, evt.x_mod, evt.eps, evt.theta, evt.phi, evt.t, evt.w, evt.Q2
        ))
    '''
    
    q2_list = []
    w_list = []
    theta_list = []
    t_list = []
    lo_eps_list = []
    hi_eps_list = []

    for evt in nlo:
        if evt.t not in t_list:
            q2_list.append(evt.Q2)
            w_list.append(evt.w)
            theta_list.append(evt.theta)
            t_list.append(evt.t)
            lo_eps_list.append(evt.eps)

    tmp_t_list = []
    for evt in nhi:
        if evt.t not in tmp_t_list:
            tmp_t_list.append(evt.t)
            hi_eps_list.append(evt.eps)
        
    t_bin_num = len(t_list)

    lo_cross_sec = np.zeros(t_bin_num, dtype=float)
    hi_cross_sec = np.zeros(t_bin_num, dtype=float)
    lo_cross_sec_err = np.zeros(t_bin_num, dtype=float)
    hi_cross_sec_err = np.zeros(t_bin_num, dtype=float)

    # Thesis, center only (removed last bin because it sucks)
    #for i in range(0, t_bin_num-1):
    for i in range(0, t_bin_num):    

        print("\n/*--------------------------------------------------*/")
        print(" Starting t-bin {0} (t={1:.4f})...".format(i+1, float(t_list[i])))
        print("\n/*--------------------------------------------------*/\n\n")

        #tcut = "t=={0} && x!=0.0 && dx>0.0".format(float(t_list[i]))
        tcut = "t=={0} && x!=0.0".format(float(t_list[i]))
        #tcut = "t=={0} && x!=0.0 && (phi>25)".format(float(t_list[i]))
        print(tcut)
        
        lo_eps = lo_eps_list[i]
        hi_eps = hi_eps_list[i]

        nlo.Draw("x:phi:dx", tcut, "goff")

        #glo_tmp = TGraphErrors(nlo.GetSelectedRows(), nlo.GetV2(), nlo.GetV1(), [0]*nlo.GetSelectedRows(), nlo.GetV3())
        glo_tmp = TGraphErrors()
        # Use SetPoint and SetPointError to fill the graph
        for j in range(nlo.GetSelectedRows()):
            
            glo_tmp.SetPoint(j, nlo.GetV2()[j], nlo.GetV1()[j])
            glo_tmp.SetPointError(j, 0, nlo.GetV3()[j])

        flo = TF1("lo_eps_fit", LT_sep_x_lo_fun, 0, 360, 4)  # 4 is the number of input params
        flo_unsep = TF1("lo_eps_unsep", LT_sep_x_lo_fun_unsep, 0, 2*PI, 4) # 4 is the number of input params

        glo = glo_tmp.Clone("glo")
            
        ave_sig_lo = glo.GetMean(2)
        err_sig_lo = glo.GetRMS(2)

        sig_lo.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_lo.GetYaxis().SetTitle("#bar{#it{#sigma}_{Low}}")
        sig_lo.SetTitle("t = {:.3f}".format(t_list[i]))
        
        sig_lo.SetPoint(sig_lo.GetN(), float(t_list[i]), ave_sig_lo)
        sig_lo.SetPointError(sig_lo.GetN()-1, 0, err_sig_lo)
        
        nhi.Draw("x:phi:dx", tcut, "goff")

        #ghi_tmp = TGraphErrors(nhi.GetSelectedRows(), nhi.GetV2(), nhi.GetV1(), 0, nhi.GetV3())
        ghi_tmp = TGraphErrors()
        # Use SetPoint and SetPointError to fill the graph
        for j in range(nhi.GetSelectedRows()):
                        
            ghi_tmp.SetPoint(j, nhi.GetV2()[j], nhi.GetV1()[j])
            ghi_tmp.SetPointError(j, 0, nhi.GetV3()[j])

        fhi = TF1("hi_eps_fit", LT_sep_x_hi_fun, 0, 360, 4)  # 4 is the number of input params
        fhi_unsep = TF1("hi_eps_unsep", LT_sep_x_hi_fun_unsep, 0, 2*PI, 4)  # 4 is the number of input params
            
        ghi = ghi_tmp.Clone("ghi")

        ave_sig_hi = ghi.GetMean(2)
        err_sig_hi = ghi.GetRMS(2)

        sig_hi.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_hi.GetYaxis().SetTitle("#bar{#it{#sigma}_{High}}")
        sig_hi.SetTitle("t = {:.3f}".format(t_list[i]))
        
        sig_hi.SetPoint(sig_hi.GetN(), float(t_list[i]), ave_sig_hi)
        sig_hi.SetPointError(sig_hi.GetN()-1, 0, err_sig_hi)
        
        g_plot_err = TGraph2DErrors()
        g_xx, g_yy, g_yy_err = ctypes.c_double(0),ctypes.c_double(0),ctypes.c_double(0)

        for ii in range(glo.GetN()):
            
            glo.GetPoint(ii, g_xx, g_yy)
            g_yy_err = math.sqrt((glo.GetErrorY(ii) / g_yy.value)**2 + (pt_to_pt_systematic_error/100)**2) * g_yy.value

            lo_cross_sec_err[i] += 1 / (g_yy_err**2)

            g_plot_err.SetPoint(g_plot_err.GetN(), g_xx, lo_eps, g_yy)
            g_plot_err.SetPointError(g_plot_err.GetN()-1, 0.0, 0.0, math.sqrt((glo.GetErrorY(ii))**2 + (pt_to_pt_systematic_error/100)**2))
            #g_plot_err.SetPointError(g_plot_err.GetN()-1, 0.0, 0.0, g_yy_err)

        for ii in range(ghi.GetN()):
            
            ghi.GetPoint(ii, g_xx, g_yy)
            g_yy_err = math.sqrt((ghi.GetErrorY(ii) / g_yy.value)**2 + (pt_to_pt_systematic_error/100)**2) * g_yy.value

            hi_cross_sec_err[i] += 1 / (g_yy_err**2)

            g_plot_err.SetPoint(g_plot_err.GetN(), g_xx, hi_eps, g_yy)
            g_plot_err.SetPointError(g_plot_err.GetN()-1, 0.0, 0.0, math.sqrt((ghi.GetErrorY(ii))**2 + (pt_to_pt_systematic_error/100)**2))
            #g_plot_err.SetPointError(g_plot_err.GetN()-1, 0.0, 0.0, g_yy_err)

        try:
            lo_cross_sec_err[i] = 1/math.sqrt(lo_cross_sec_err[i])            
            hi_cross_sec_err[i] = 1/math.sqrt(hi_cross_sec_err[i])
        except ZeroDivisionError:
            lo_cross_sec_err[i] = -1000
            hi_cross_sec_err[i] = -1000
            
        g_plot_err.SetFillColor(29)
        g_plot_err.SetMarkerSize(0.8)
        g_plot_err.SetMarkerStyle(20)
        g_plot_err.SetMarkerColor(ROOT.kRed)
        g_plot_err.SetLineColor(ROOT.kBlue-3)
        g_plot_err.SetLineWidth(2)
        
        # x->phi, y->eps, PI/180 = 0.017453
        fff2 = TF2("fff2", "[0] + y*[1] + sqrt(2*y*(1+y))*cos(x*0.017453)*[2] + y*cos(2*x*0.017453)*[3]", 0, 360, 0.0, 1.0)
        #fff2 = TF2("fff2", "[0] + y*[1] + sqrt(2*y*(1+y))*cos(x*0.017453)*[2] + y*cos(2*x*0.017453)*[3]", 0, 360, LOEPS-0.2*LOEPS, HIEPS+0.2*HIEPS)

        sigL_change = TGraphErrors()
        sigT_change = TGraphErrors()
        sigLT_change = TGraphErrors()
        sigTT_change = TGraphErrors()

        #########
        # Fit 1 #
        #########

        #print("\n/*--------------------------------------------------*/")
        #print(" Fitting Step 1")
        #print(" Fit L and T, while Fix LT and TT")

        # Set parameter 0 and 1
        fff2.SetParameter(0, 1)
        fff2.SetParLimits(0, 0.0, 1000)

        fff2.SetParameter(1, 1)
        fff2.SetParLimits(1, 0.0, 1000)

        # Fix parameter 2 and 3
        fff2.FixParameter(2, 0.0)
        fff2.FixParameter(3, 0.0)

        # Options: M-> Improve fit info splash, R-> Use range specified, Q-> Quiet splash
        g_plot_err.Fit(fff2, "MRQ")

        sigL_change.SetTitle("t = {:.3f}".format(t_list[i]))
        sigL_change.GetXaxis().SetTitle("Fit Step")
        sigL_change.GetYaxis().SetTitle("#it{#sigma}_{L}")
        
        # sigL_change
        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN() + 1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN() - 1, 0, fff2.GetParError(1))

        sigT_change.SetTitle("t = {:.3f}".format(t_list[i]))
        sigT_change.GetXaxis().SetTitle("Fit Step")
        sigT_change.GetYaxis().SetTitle("#it{#sigma}_{T}")
        
        # sigT_change
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN() + 1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN() - 1, 0, fff2.GetParError(0))

        #########            
        # Fit 2 #
        #########

        #print("\n/*--------------------------------------------------*/")
        #print(" Fitting Step 2")
        #print(" Fit LT, while Fix L, T, and TT")

        # Fix parameter 0, 1, and 3
        fff2.FixParameter(0, fff2.GetParameter(0))
        fff2.FixParameter(1, fff2.GetParameter(1))
        fff2.FixParameter(3, fff2.GetParameter(3))

        # Release parameter 2
        fff2.ReleaseParameter(2)

        # Set parameter 2
        fff2.SetParameter(2, 0.0)
        fff2.SetParLimits(2, -10, 50)

        # Options: M-> Improve fit info splash, R-> Use range specified, Q-> Quiet splash
        g_plot_err.Fit(fff2, "MRQ")

        # Update sigL_change and sigT_change
        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN() + 1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN() - 1, 0, fff2.GetParError(1))

        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN() + 1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN() - 1, 0, fff2.GetParError(0))
        
        #########
        # Fit 3 #
        #########

        #print("\n/*--------------------------------------------------*/")
        #print(" Fitting Step 3")
        #print(" Fit L and T, while Fix LT and TT")

        # Release parameter 0 and 1
        fff2.ReleaseParameter(0)
        fff2.ReleaseParameter(1)

        # Set parameter 0 and 1
        fff2.SetParameter(0, fff2.GetParameter(0))
        fff2.SetParameter(1, fff2.GetParameter(1))

        # Fix parameter 2 and 3
        fff2.FixParameter(2, fff2.GetParameter(2))
        fff2.FixParameter(3, fff2.GetParameter(3))

        # Options: M-> Improve fit info splash, R-> Use range specified, Q-> Quiet splash
        g_plot_err.Fit(fff2, "MRQ")

        # sigL_change
        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN() + 1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN() - 1, 0, fff2.GetParError(1))

        # sigT_change
        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN() + 1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN() - 1, 0, fff2.GetParError(0))

        #########
        # Fit 4 #
        #########

        #print("\n/*--------------------------------------------------*/")
        #print(" Fitting Step 4")
        #print(" Fit TT, while Fix T, L, and LT")

        # Fix parameter 0, 1, and 2
        fff2.FixParameter(0, fff2.GetParameter(0))
        fff2.FixParameter(1, fff2.GetParameter(1))
        fff2.FixParameter(2, fff2.GetParameter(2))

        # Release parameter 3
        fff2.ReleaseParameter(3)

        # Set parameter 3
        fff2.SetParameter(3, 0.0)
        fff2.SetParLimits(3, -10, 50)

        # Options: M-> Improve fit info splash, R-> Use range specified, Q-> Quiet splash
        g_plot_err.Fit(fff2, "MRQ")

        # Update sigL_change and sigT_change
        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN() + 1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN() - 1, 0, fff2.GetParError(1))

        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN() + 1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN() - 1, 0, fff2.GetParError(0))
        
        #########
        # Fit 5 #
        #########

        #print("\n/*--------------------------------------------------*/")
        #print(" Fitting Step 5")
        #print(" Fit T and L, while Fix LT and TT")

        # Release parameter 0 and 1
        fff2.ReleaseParameter(0)
        fff2.ReleaseParameter(1)

        # Set parameter 0 and 1
        fff2.SetParameter(0, fff2.GetParameter(0))
        fff2.SetParameter(1, fff2.GetParameter(1))

        # Fix parameter 2 and 3
        fff2.FixParameter(2, fff2.GetParameter(2))
        fff2.FixParameter(3, fff2.GetParameter(3))

        # Options: M-> Improve fit info splash, R-> Use range specified, Q-> Quiet splash
        g_plot_err.Fit(fff2, "MRQ")

        # Update sigL_change and sigT_change
        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN() + 1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN() - 1, 0, fff2.GetParError(1))

        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN() + 1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN() - 1, 0, fff2.GetParError(0))

        #########
        # Fit 6 #
        #########

        #print("\n/*--------------------------------------------------*/")
        #print(" Fitting Step 4")
        #print(" Fit LT and TT, while Fix L and T")

        # Fix parameter 0 and 1
        fff2.FixParameter(0, fff2.GetParameter(0))
        fff2.FixParameter(1, fff2.GetParameter(1))

        # Release and set parameters
        fff2.ReleaseParameter(2)
        fff2.ReleaseParameter(3)
        fff2.SetParameter(2, fff2.GetParameter(2))
        fff2.SetParameter(3, fff2.GetParameter(3))

        # Fit
        # Options: M-> Improve fit info splash, R-> Use range specified, Q-> Quiet splash
        g_plot_err.Fit("fff2", "MRQ")
        
        # Update sigL_change and sigT_change
        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN() + 1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN() - 1, 0, fff2.GetParError(1))

        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN() + 1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN() - 1, 0, fff2.GetParError(0))
        
        ####################
        # Last Step, Fit 7 #
        ####################

        #print("\n/*--------------------------------------------------*/")
        #print(" Last Step")
        #print(" Fit All")

        # Release all parameters
        fff2.ReleaseParameter(0)
        fff2.ReleaseParameter(1)
        fff2.ReleaseParameter(2)
        fff2.ReleaseParameter(3)

        # Set parameter limits
        fff2.SetParLimits(0, 0.0, 1000)
        fff2.SetParLimits(1, 0.0, 1000)
        fff2.SetParLimits(2, -10, 50)
        fff2.SetParLimits(3, -10, 50)

        # Fit
        # Options: M-> Improve fit info splash, R-> Use range specified, Q-> Quiet splash
        g_plot_err.Fit(fff2, "MRQ")

        # Update sigL_change and sigT_change
        sigL_change.SetPoint(sigL_change.GetN(), sigL_change.GetN() + 1, fff2.GetParameter(1))
        sigL_change.SetPointError(sigL_change.GetN() - 1, 0, fff2.GetParError(1))

        sigT_change.SetPoint(sigT_change.GetN(), sigT_change.GetN() + 1, fff2.GetParameter(0))
        sigT_change.SetPointError(sigT_change.GetN() - 1, 0, fff2.GetParError(0))
        
        c1 =  TCanvas()

        c1.Update()
                
        # Print plots for c1 canvases
        g_plot_err.Draw("perr")
        g_plot_err.SetTitle("t = {:.3f}".format(t_list[i]))
        g_plot_err.GetXaxis().SetTitle("#it{#phi} [degree]")
        g_plot_err.GetXaxis().CenterTitle()
        g_plot_err.GetXaxis().SetTitleOffset(1.5)
        g_plot_err.GetYaxis().SetTitle("Epsilon")
        g_plot_err.GetYaxis().CenterTitle()
        g_plot_err.GetYaxis().SetTitleOffset(1.5)
        g_plot_err.GetZaxis().SetTitle("#it{#sigma}")
        g_plot_err.GetZaxis().CenterTitle()
        g_plot_err.GetZaxis().SetTitleOffset(1.5)
        if i == 0:
            c1.Print(outputpdf+'(')
        else:
            c1.Print(outputpdf)
        c1.Clear()

        c2 =  TCanvas()
                
        # Update c2
        c2.Update()

        # Go to the c2 canvas
        c2.cd()
        
        # Set properties for glo
        glo.SetMarkerStyle(5)
        glo.GetXaxis().SetLimits(0, 360)

        # Set properties for ghi
        ghi.SetMarkerColor(2)
        ghi.SetLineColor(2)
        ghi.SetMarkerStyle(4)

        # Create TMultiGraph and add glo, ghi
        g = ROOT.TMultiGraph()
        g.Add(glo)
        g.Add(ghi)

        # Draw TMultiGraph
        g.Draw("AP")

        # Set properties for the TMultiGraph
        #g.GetHistogram().SetMinimum(glo.GetMinimum() * 0.8)
        #g.GetHistogram().SetMaximum(glo.GetMaximum() * 1.2)
        g.GetYaxis().SetTitle("Unseparated Cross Section [nb/GeV^{2}]")
        g.GetYaxis().CenterTitle()
        g.GetYaxis().SetTitleOffset(1.2)
        g.GetYaxis().SetTitleSize(0.04)

        g.GetXaxis().SetTitle("#it{#phi} [degree]")
        g.GetXaxis().CenterTitle()
        g.GetXaxis().SetLimits(0, 360)
        g.GetXaxis().SetTitleSize(0.04)
        
        # Update canvas c2
        c2.Update()

        # Fix parameters for flo, flo_unsep, fhi, and fhi_unsep
        flo.FixParameter(0, fff2.GetParameter(0)) # sigT
        flo.FixParameter(1, fff2.GetParameter(1)) # sigL
        flo.FixParameter(2, fff2.GetParameter(2)) # sigLT
        flo.FixParameter(3, fff2.GetParameter(3)) # sigTT

        flo_unsep.FixParameter(0, fff2.GetParameter(0)) # sigT
        flo_unsep.FixParameter(1, fff2.GetParameter(1)) # sigL
        flo_unsep.FixParameter(2, fff2.GetParameter(2)) # sigLT
        flo_unsep.FixParameter(3, fff2.GetParameter(3)) # sigTT

        fhi.FixParameter(0, fff2.GetParameter(0)) # sigT
        fhi.FixParameter(1, fff2.GetParameter(1)) # sigL
        fhi.FixParameter(2, fff2.GetParameter(2)) # sigLT
        fhi.FixParameter(3, fff2.GetParameter(3)) # sigTT

        fhi_unsep.FixParameter(0, fff2.GetParameter(0)) # sigT
        fhi_unsep.FixParameter(1, fff2.GetParameter(1)) # sigL
        fhi_unsep.FixParameter(2, fff2.GetParameter(2)) # sigLT
        fhi_unsep.FixParameter(3, fff2.GetParameter(3)) # sigTT

        # Options: M-> Improve fit info splash, R-> Use range specified, Q-> Quiet splash
        glo.Fit(flo, "MRQ")
        ghi.Fit(fhi, "MRQ")
        
        # Set line properties for flo and fhi
        flo.SetLineColor(1)
        fhi.SetLineColor(2)
        flo.SetLineWidth(2)
        fhi.SetLineWidth(2)
        fhi.SetLineStyle(2)

        # Set line color for ghi
        ghi.SetLineColor(2)

        # Draw flo and fhi on the same canvas
        flo.Draw("same")
        fhi.Draw("same")

        # Create and draw TLegend
        leg = ROOT.TLegend(0.85, 0.85, 0.99, 0.99)
        leg.SetFillColor(0)
        leg.SetMargin(0.4)
        leg.SetTextSize(0.04)
        leg.AddEntry(glo, "Low #it{#font[120]{e}}", "pl")
        leg.AddEntry(ghi, "High #it{#font[120]{e}}", "pl")
        leg.Draw()

        # Check the fit status for 'flo'
        flo_status = flo.GetNDF()  # GetNDF() returns the number of degrees of freedom
        flo_status_message = "Not Fitted" if flo_status == 0 else "Fit Successful"

        # Check the fit status for 'fhi'
        fhi_status = fhi.GetNDF()  # GetNDF() returns the number of degrees of freedom
        fhi_status_message = "Not Fitted" if fhi_status == 0 else "Fit Successful"        
        
        # Create TText for fit status
        fit_status = ROOT.TLatex()
        fit_status.SetTextSize(0.04)
        fit_status.DrawTextNDC(0.15, 0.85, "t={:.3f}, Q2={:.1f}, W={:.2f}".format(t_list[i], float(q2_set.replace("p",".")), float(w_set.replace("p","."))))
        #fit_status.DrawTextNDC(0.15, 0.80, "Low Fit Status: " + flo_status_message)
        #fit_status.DrawTextNDC(0.15, 0.75, "High Fit Status: " + fhi_status_message)

        # Adjust the maximum and minimum of glo based on ghi values
        if ghi.GetMaximum() > glo.GetMaximum():
            glo.SetMaximum(ghi.GetMaximum() * 1.1)

        if ghi.GetMinimum() < glo.GetMinimum():
            glo.SetMinimum(ghi.GetMinimum() * 0.9)

        sig_diff_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_diff_g.GetYaxis().SetTitle("#Delta#sigma_{ave}/#Delta#epsilon")

        try:
            sig_diff = (ave_sig_hi-ave_sig_lo)/eps_diff
            sig_diff_g.SetPoint(sig_diff_g.GetN(), float(t_list[i]), sig_diff)

            sig_diff_err = sig_diff*math.sqrt((err_sig_hi/ave_sig_hi)**2+(err_sig_lo/ave_sig_lo)**2)
            sig_diff_g.SetPointError(sig_diff_g.GetN()-1, 0, sig_diff_err)
        except ZeroDivisionError:
            sig_diff = 0.0
            sig_diff_g.SetPoint(sig_diff_g.GetN(), float(t_list[i]), sig_diff)

            sig_diff_err = 0.0
            sig_diff_g.SetPointError(sig_diff_g.GetN()-1, 0, sig_diff_err)            
            
        # Define variables for cross sections and errors
        sig_l, sig_t, sig_lt, sig_tt = fff2.GetParameter(1), fff2.GetParameter(0), fff2.GetParameter(2), fff2.GetParameter(3)
        sig_l_err, sig_t_err, sig_lt_err, sig_tt_err = fff2.GetParError(1), fff2.GetParError(0), fff2.GetParError(2), fff2.GetParError(3)

        # Print values to console
        print("\nBin {}: Outputting...  ".format(i+1), "sig_l: ", sig_l, "sig_t: ", sig_t, \
              "sig_lt: ", sig_lt, "sig_tt: ", sig_tt, \
              "t: ", t_list[i], "W: ", w_list[i], "Q2:", q2_list[i], \
              "eps_lo: ", lo_eps_list[i], "eps_hi: ", hi_eps_list[i])

        fn_sep = "{}/src/{}/xsects/x_sep.{}_Q{}W{}.dat".format(LTANAPATH, ParticleType, polID, Q2.replace("p",""), W.replace("p",""))
        try:
            mode = 'w' if i == 0 else 'a'
            with open(fn_sep, mode) as f:
                # Write values to output file
                f.write("{} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
                    sig_l, sig_l_err, sig_t, sig_t_err, sig_lt, sig_lt_err, sig_tt, sig_tt_err,

                    fff2.GetChisquare(), t_list[i], w_list[i], q2_list[i], theta_list[i]
                ))

        except IOError:
            print("Error writing to file {}.".format(fn_sep))
            
        # Delete g_plot_err
        del g_plot_err
        
        g_sig_l_total.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        g_sig_l_total.GetYaxis().SetTitle("#it{#sigma}_{L} [nb/GeV^{2}]")
        
        g_sig_t_total.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        g_sig_t_total.GetYaxis().SetTitle("#it{#sigma}_{T} [nb/GeV^{2}]")

        g_sig_lt_total.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        g_sig_lt_total.GetYaxis().SetTitle("#it{#sigma}_{LT} [nb/GeV^{2}]")

        g_sig_tt_total.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        g_sig_tt_total.GetYaxis().SetTitle("#it{#sigma}_{TT} [nb/GeV^{2}]")
        
        # Set points and errors for g_sig_l_total, g_sig_t_total, g_sig_lt_total, and g_sig_tt_total
        g_sig_l_total.SetPoint(g_sig_l_total.GetN(), float(t_list[i]), sig_l)
        g_sig_l_total.SetPointError(g_sig_l_total.GetN() - 1, 0, sig_l_err)

        g_sig_t_total.SetPoint(g_sig_t_total.GetN(), float(t_list[i]), sig_t)
        g_sig_t_total.SetPointError(g_sig_t_total.GetN() - 1, 0, sig_t_err)

        g_sig_lt_total.SetPoint(g_sig_lt_total.GetN(), float(t_list[i]), sig_lt)
        g_sig_lt_total.SetPointError(g_sig_lt_total.GetN() - 1, 0, sig_lt_err)

        g_sig_tt_total.SetPoint(g_sig_tt_total.GetN(), float(t_list[i]), sig_tt)
        g_sig_tt_total.SetPointError(g_sig_tt_total.GetN() - 1, 0, sig_tt_err)

        sig_L_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_L_g.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{L} [nb/GeV^{2}]")
        sig_L_g.SetTitle("t = {:.3f}".format(t_list[i]))
        
        sig_T_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_T_g.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{T} [nb/GeV^{2}]")
        sig_T_g.SetTitle("t = {:.3f}".format(t_list[i]))

        sig_LT_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_LT_g.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{LT} [nb/GeV^{2}]")
        sig_LT_g.SetTitle("t = {:.3f}".format(t_list[i]))

        sig_TT_g.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
        sig_TT_g.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right)_{TT} [nb/GeV^{2}]")
        sig_TT_g.SetTitle("t = {:.3f}".format(t_list[i]))        
        
        # Set points and errors for sig_L_g, sig_T_g, sig_LT_g, and sig_TT_g
        sig_L_g.SetPoint(i, float(t_list[i]), sig_l)
        sig_T_g.SetPoint(i, float(t_list[i]), sig_t)
        sig_LT_g.SetPoint(i, float(t_list[i]), sig_lt)
        sig_TT_g.SetPoint(i, float(t_list[i]), sig_tt)

        sig_L_g.SetPointError(i, 0.0, sig_l_err)
        sig_T_g.SetPointError(i, 0.0, sig_t_err)
        sig_LT_g.SetPointError(i, 0.0, sig_lt_err)
        sig_TT_g.SetPointError(i, 0.0, sig_tt_err)

        # Create TCanvas
        c4 = TCanvas()

        # Draw and save plots
        sigL_change.Draw("a*")
        c4.Print(outputpdf)

        sigT_change.Draw("a*")
        c4.Print(outputpdf)    
        c4.Clear()

        # Adjust top and right margins for c2 canvas
        c2.SetTopMargin(0.03)
        c2.SetRightMargin(0.03)

        # Print plots for c2 canvases
        c2.Print(outputpdf)
        c2.Clear()
            
        # Create TCanvas
        c3 = TCanvas()

        # Draw and save plots for sig_L_g, sig_T_g, sig_LT_g, and sig_TT_g
        sig_L_g.Draw("a*")
        c3.Print(outputpdf)

        sig_T_g.Draw("a*")
        c3.Print(outputpdf)

        sig_LT_g.Draw("a*")
        c3.Print(outputpdf)

        sig_TT_g.Draw("a*")
        c3.Print(outputpdf)
        c3.Clear()

        # Create TCanvas
        c5 = TCanvas()

        sig_lo.Draw("a*")
        c5.Print(outputpdf)
        
        sig_hi.Draw("a*")
        c5.Print(outputpdf)

        sig_diff_g.Draw("a*")
        c5.Print(outputpdf)
        c5.Clear()
        
        # Calculate integrated cross sections
        lo_cross_sec[i] = flo_unsep.Integral(0, 2*PI) / (2*PI)
        hi_cross_sec[i] = fhi_unsep.Integral(0, 2*PI) / (2*PI)
        
        g_unsep_lo.SetPoint(g_unsep_lo.GetN(), LOEPS, lo_cross_sec[i])
        g_unsep_lo.SetPointError(g_unsep_lo.GetN()-1, 0, lo_cross_sec_err[i])        
        
        g_unsep_hi.SetPoint(g_unsep_hi.GetN(), HIEPS, hi_cross_sec[i])
        g_unsep_hi.SetPointError(g_unsep_hi.GetN()-1, 0, hi_cross_sec_err[i])                
        
        fn_unsep = "{}/src/{}/xsects/unsep_Q{}W{}.csv".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))
        try:
            mode = 'w' if i == 0 else 'a'
            with open(fn_unsep, mode) as f:
                # Write values to output file
                f.write("{:.3f} {:.3f} {:.3f} {:.3f} {:.4f}\n".format(
                    lo_cross_sec[i], lo_cross_sec_err[i], hi_cross_sec[i], hi_cross_sec_err[i], t_list[i]
                ))

        except IOError:
            print("Error writing to file {}.".format(fn_unsep))

        # Delete canvases
        del c1
        del c2
        del c3
        del c4
        del c5
            
    return t_list

fn_lo =  "{}/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(LTANAPATH, ParticleType, polID, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100)
fn_hi =  "{}/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(LTANAPATH, ParticleType, polID, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100)

g_sig_l_total = TGraphErrors()
g_sig_t_total = TGraphErrors()
g_sig_lt_total = TGraphErrors()
g_sig_tt_total = TGraphErrors()

g_unsep_lo = TGraphErrors()
g_unsep_hi = TGraphErrors()

t_list = single_setting(Q2, W, fn_lo, fn_hi) # Main function that performs fitting

f_lin_l = TF1("f_lin_l", "[0]+[1]*x", 0.0, 1.0)
f_lin_t = TF1("f_lin_t", "[0]+[1]*x", 0.0, 1.0)

# Define the number of events in 'lo'
num_events = g_unsep_lo.GetN()

# Loop over each event in 'lo'
for i in range(num_events):

    # Create a canvas
    c = ROOT.TCanvas("c", "c", 800, 600)
    
    # Create a new TMultiGraph for each event
    g_unsep_mult = ROOT.TMultiGraph()
    
    # Create TGraphErrors for 'lo' event
    x_lo, y_lo = ctypes.c_double(0), ctypes.c_double(0)
    g_unsep_lo.GetPoint(i, x_lo, y_lo)
    x_err_lo = g_unsep_lo.GetErrorX(i)
    y_err_lo = g_unsep_lo.GetErrorY(i)
    g_lo_event = ROOT.TGraphErrors(1, array('d', [x_lo.value]), array('d', [y_lo.value]), array('d', [x_err_lo]), array('d', [y_err_lo]))
    
    # Add 'lo' event to the TMultiGraph
    g_unsep_mult.Add(g_lo_event)
    
    # Create TGraphErrors for 'hi' event
    x_hi, y_hi = ctypes.c_double(0), ctypes.c_double(0)
    g_unsep_hi.GetPoint(i, x_hi, y_hi)
    x_err_hi = g_unsep_hi.GetErrorX(i)
    y_err_hi = g_unsep_hi.GetErrorY(i)
    g_hi_event = ROOT.TGraphErrors(1, array('d', [x_hi.value]), array('d', [y_hi.value]), array('d', [x_err_hi]), array('d', [y_err_hi]))
    
    # Add 'hi' event to the TMultiGraph
    g_unsep_mult.Add(g_hi_event)

    g_lo_event.SetMarkerStyle(5)
    g_hi_event.SetMarkerStyle(5)
    
    # Draw TMultiGraph for the current event
    g_unsep_mult.Draw("AP")
    
    # Set axis titles and offsets
    g_unsep_mult.GetYaxis().SetTitle("Unseparated Cross Section [nb/GeV^{2}]")
    g_unsep_mult.GetYaxis().SetTitleOffset(1.2)
    g_unsep_mult.GetXaxis().SetTitle("#epsilon")
    g_unsep_mult.GetXaxis().SetTitleOffset(1.2)
    
    # Fit functions to 'lo' and 'hi' events
    f_lin = ROOT.TF1("f_lin", "[0]*x + [1]", 0, 1)  # Define fit function for 'lo'
    g_unsep_mult.Fit(f_lin, "MRQ")
        
    # Set line properties for 'lo' and 'hi' fits
    f_lin.SetLineColor(2)
    f_lin.SetLineWidth(2)    
    f_lin.Draw("same")
    
    # Create and draw TLegend
    leg = ROOT.TLegend(0.85, 0.85, 0.99, 0.99)
    leg.SetFillColor(0)
    leg.SetMargin(0.4)
    leg.AddEntry(g_hi_event, "t = {:.3f}".format(t_list[i]), "")
    leg.AddEntry(g_lo_event, "Low #it{#font[120]{e}}", "p")
    leg.AddEntry(g_hi_event, "High #it{#font[120]{e}}", "p")
    leg.Draw()
    
    c.Print(outputpdf)
    c.Clear()

    del c
    
f_exp_l = TF1("f_exp_l", "[0]*exp(-[1]*x)", 0.0, 2.0)
f_exp_t = TF1("f_exp_t", "[0]*exp(-[1]*x)", 0.0, 2.0)

c_total_l_t = TCanvas()

# Create TMultiGraph and add glo, ghi
g_sig_mult = ROOT.TMultiGraph()
g_sig_mult.Add(g_sig_l_total)
g_sig_mult.Add(g_sig_t_total)

#g_sig_l_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_l_total.GetYaxis().SetRangeUser(-5,30)
#g_sig_t_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_t_total.GetYaxis().SetRangeUser(-5,30)

g_sig_mult.Draw("AP")

g_sig_mult.GetYaxis().SetTitle("#left(#frac{#it{d#sigma}}{#it{dt}}#right) [nb/GeV^{2}]")
g_sig_mult.GetYaxis().SetTitleOffset(1.2)

g_sig_mult.GetXaxis().SetTitle("#it{-t} [GeV^{2}]")
g_sig_mult.GetXaxis().SetTitleOffset(1.2)

g_sig_l_total.Fit(f_exp_l, "MRQ")
g_sig_t_total.Fit(f_exp_t, "MRQ")

# Set properties for g_sig_l_total and g_sig_t_total
g_sig_l_total.SetLineColor(1)
g_sig_l_total.SetMarkerStyle(5)
g_sig_t_total.SetLineColor(2)
g_sig_t_total.SetMarkerColor(2)
g_sig_t_total.SetMarkerStyle(4)

# Set line properties for g_sig_l_total and g_sig_t_total
f_exp_l.SetLineColor(1)
f_exp_l.SetLineWidth(2)
f_exp_t.SetLineColor(2)
f_exp_t.SetLineWidth(2)
f_exp_t.SetLineStyle(2)

# Draw f_exp_l and f_exp_t on the same canvas
f_exp_l.Draw("same")
f_exp_t.Draw("same")

# Create and draw TLegend
leg = ROOT.TLegend(0.85, 0.85, 0.99, 0.99)
leg.SetFillColor(0)
leg.SetMargin(0.4)
leg.AddEntry(g_sig_l_total, "#it{#sigma}_{L}", "p")
leg.AddEntry(g_sig_t_total, "#it{#sigma}_{T}", "p")
leg.Draw()

c_total_l_t.Print(outputpdf)
c_total_l_t.Clear()

ROOT.gStyle.SetOptFit(1)

# Set the size of the statbox
ROOT.gStyle.SetStatW(0.15)  # Set width of statbox
ROOT.gStyle.SetStatH(0.15)  # Set height of statbox

# Set the position of the statbox
#ROOT.gStyle.SetStatX(0.8)  # Set x position of statbox (right edge)
#ROOT.gStyle.SetStatY(0.8)  # Set y position of statbox (top edge)

f_exp = TF1("f_exp", "[0]*exp(-[1]*x)", 0.0, 2.0)

c_total = TCanvas()

#g_sig_l_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_l_total.GetYaxis().SetRangeUser(-5,30)
g_sig_l_total.Draw("A*")
g_sig_l_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf)
c_total.Clear()

#g_sig_t_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_t_total.GetYaxis().SetRangeUser(-5,30)
g_sig_t_total.SetMarkerColor(1)
g_sig_t_total.SetLineColor(1)
g_sig_t_total.Draw("A*")
g_sig_t_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf)
c_total.Clear()

#g_sig_lt_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_lt_total.GetYaxis().SetRangeUser(-5,3)
g_sig_lt_total.Draw("A*")
g_sig_lt_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf)
c_total.Clear()

#g_sig_tt_total.GetXaxis().SetRangeUser(0.0,0.4)
#g_sig_tt_total.GetYaxis().SetRangeUser(-3,1)
g_sig_tt_total.Draw("A*")
g_sig_tt_total.Fit(f_exp, "MRQ")
c_total.Print(outputpdf+')')
c_total.Clear()
