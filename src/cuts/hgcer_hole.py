#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-13 23:39:44 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TCutG

def apply_HGCer_hole_cut(Q2, W, EPSSET, simc=False):

    if "high" in EPSSET:
        # Defined HGCer Geometric cuts
        cutg = TCutG("cutg",21)
        cutg.SetVarX("P_hgcer_xAtCer")
        cutg.SetVarY("P_hgcer_yAtCer")
        cutg.SetPoint(0, 2+10, -25)
        cutg.SetPoint(1, 2+10, -2)
        cutg.SetPoint(2, 2+10, -1.5)
        cutg.SetPoint(3, 3+10, 0)
        cutg.SetPoint(4, 3+10, 1)
        cutg.SetPoint(5, 3+10, 2.3)
        cutg.SetPoint(6, 3+10, 3.0)
        cutg.SetPoint(7, 2+10, 4.5)
        cutg.SetPoint(8, 2+10, 5)
        cutg.SetPoint(9, 2+10, 25)
        cutg.SetPoint(10, 0+10, 25.5)
        cutg.SetPoint(11, 0+10, 5.5)
        cutg.SetPoint(12, 1+10, 4)
        cutg.SetPoint(13, -1+10, 3)
        cutg.SetPoint(14, -2+10, 2)
        cutg.SetPoint(15, -2.3+10, 1)
        cutg.SetPoint(16, -1.5+10, 0)
        cutg.SetPoint(17, -1+10, -1)
        cutg.SetPoint(18, 0.5+10, -2)
        cutg.SetPoint(19, 0.5+10, -25)
        cutg.SetPoint(20, 2+10, -25)
        
    else:
        # Defined HGCer Geometric cuts
        cutg = TCutG("cutg",21)
        cutg.SetVarX("P_hgcer_xAtCer")
        cutg.SetVarY("P_hgcer_yAtCer")
        cutg.SetPoint(0, 2, -25)
        cutg.SetPoint(1, 2, -2)
        cutg.SetPoint(2, 2, -1.5)
        cutg.SetPoint(3, 3, 0)
        cutg.SetPoint(4, 3, 1)
        cutg.SetPoint(5, 3, 2.3)
        cutg.SetPoint(6, 3, 3.0)
        cutg.SetPoint(7, 2, 4.5)
        cutg.SetPoint(8, 2, 5)
        cutg.SetPoint(9, 2, 25)
        cutg.SetPoint(10, 0, 25.5)
        cutg.SetPoint(11, 0, 5.5)
        cutg.SetPoint(12, 1, 4)
        cutg.SetPoint(13, -1, 3)
        cutg.SetPoint(14, -2, 2)
        cutg.SetPoint(15, -2.3, 1)
        cutg.SetPoint(16, -1.5, 0)
        cutg.SetPoint(17, -1, -1)
        cutg.SetPoint(18, 0.5, -2)
        cutg.SetPoint(19, 0.5, -25)
        cutg.SetPoint(20, 2, -25)    

    '''
    # Defined HGCer Geometric cuts
    cutg = TCutG("cutg")
    cutg.SetVarX("P_hgcer_yAtCer")
    cutg.SetVarY("P_hgcer_xAtCer")
    '''
    return cutg
