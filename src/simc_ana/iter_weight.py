#! /usr/bin/python

#
# Description: Adapted from fortran code wt28_3.f
# ================================================================
# Time-stamp: "2023-09-18 11:51:28 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TCanvas, TH1D, TH2D, gStyle, gPad, TPaveText, TArc, TGraphErrors, TGraphPolar, TFile, TLegend, TMultiGraph, TLine
import array
import sys, math, os, subprocess

##################################################################################################################################################
# Importing utility functions

sys.path.append("utility")
from utility import run_fortran

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

def iter_weight(param_file, fort_param, simc_root, inpDict, phi_setting):
    '''
    # Fortran script converted to python
    
    # Define constants
    pi = 3.14159
    mtar_gev = 0.93827231
    
    # Define parameters
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12 = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    
    if abs(q2_set - 245) < 1:
        p1 =  0.25961E+02 
        p2 = -0.10000E+02 
        p3 = -0.15838E+02 
        p4 =  0.00000E+00 
        p5 =  0.46859E+02 
        p6 = -0.30000E+02 
        p7 = -0.33572E+01 
        p8 =  0.00000E+00 
        p9 =  0.10000E+04 
        p10 = -0.28000E+02 
        p11 =  0.35000E+01 
        p12 = -0.67276E+02 
    else:
        print('wtn: q2 error', q2_set)
        return

    # Parameterization based upon Fpi-1 pi+ IT25, 12.04.18
    # Revised for IT21, 12.11.06
    q2_gev = float(q2_set) / 100.0
    tav = (0.0735 + 0.028 * math.log(q2_gev)) * q2_gev
    ftav = (ti - tav) / tav
    ft = ti / ((ti + 0.139570**2)**2)

    nsigl = (p1 + p2 * math.log(q2_gev)) * math.exp((p3 + p4 * math.log(q2_gev)) * (ti - 0.2))
    nsigt = p5 + p6 * math.log(q2_gev) + (p7 + p8 * math.log(q2_gev)) * ftav

    nsiglt = (p9 * math.exp(p10 * ti) + p11 / ti) * math.sin(thetacmi)
    nsigtt = (p12 * q2_gev * math.exp(-q2_gev)) * ft * math.sin(thetacmi)**2

    nsig219 = (nsigt + epsiloni * nsigl + epsiloni * math.cos(2.0 * phicmi) * nsigtt
              + math.sqrt(2.0 * epsiloni * (1.0 + epsiloni)) * math.cos(phicmi) * nsiglt) / 1.0

    wfactor = 1.0 / (Wcmi**2 - mtar_gev**2)**2
    nsig = nsig219 * wfactor

    nsig = nsig / (2.0 * pi * 1.0E+06)  # dsig/dtdphicm in microbarns/MeV**2/rad

    wtn = Weight * nsig / dsigdt

    wtn_limit = 0.20
    if wtn < wtn_limit and wtn > 0.0:
        continue
    else:
        wtn = 0.0

    return
    '''

    formatted_date  = inpDict["formatted_date"]
    Q2 = inpDict["Q2"].replace("p","")

    # Define diamond cut parameters
    a1 = inpDict["a1"]
    b1 = inpDict["b1"]
    a2 = inpDict["a2"]
    b2 = inpDict["b2"]
    a3 = inpDict["a3"]
    b3 = inpDict["b3"]
    a4 = inpDict["a4"]
    b4 = inpDict["b4"]
    
    param_arr = []
    with open(param_file, 'r') as f:
        for i, line in enumerate(f):
            columns = line.split()
            param_arr.append(str(columns[0]))    

    if not os.path.isfile(simc_root):
        print("\n\nERROR: No simc file found called {}\n\n".format(simc_root))

    InFile_SIMC = TFile.Open(simc_root, "OPEN")
    TBRANCH_SIMC  = InFile_SIMC.Get("h10")

    ################################################################################################################################################
    # Run over simc root branch to determine new weight

    print("\nGrabbing %s simc..." % phi_setting)
    for i,evt in enumerate(TBRANCH_SIMC):

      # Progress bar
      Misc.progressBar(i, TBRANCH_SIMC.GetEntries(),bar_length=25)

      # Define the acceptance cuts  
      SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
      HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)
      if ( a1 == 0.0 and  b1 == 0.0 and  a2 == 0.0 and  b2 == 0.0 and  a3 == 0.0 and  b3 == 0.0 and  a4 == 0.0 and  b4 == 0.0):
          Diamond = True
      else:
          try:
              Diamond = (evt.W/evt.Q2>a1+b1/evt.Q2) & (evt.W/evt.Q2<a2+b2/evt.Q2) & (evt.W/evt.Q2>a3+b3/evt.Q2) & (evt.W/evt.Q2<a4+b4/evt.Q2)
          except ZeroDivisionError:
              Diamond = False

      #........................................

      #Fill SIMC events
      if(HMS_Acceptance & SHMS_Acceptance & Diamond):

          # thetacm and phicm are correct, the next line is just for testingx
          #inp_fort_param = '{} {} {} {} {} {} {} {} {} '.format(Q2, evt.Q2, evt.W, evt.t, evt.epsilon, evt.thetacm, evt.phicm, evt.sigcm, evt.Weight)+' '.join(param_arr)
          inp_fort_param = '{} {} {} {} {} {} {} {} {} '.format(Q2, evt.Q2, evt.W, evt.t, evt.epsilon, evt.thetapq, evt.phipq, evt.sigcm, evt.Weight)+ \
                           ' '.join(param_arr)
          print(inp_fort_param)
              
          # Get the standard output and standard error
          stdout, stderr = run_fortran(fort_param, inp_fort_param)

          # Print the output
          print(stdout)
          #print(stderr)
    
    H_Weight_SIMC  = TH1D("H_Weight_SIMC","{} Weight".format(formatted_date), 500, 0, 1e-8)    
    
        
