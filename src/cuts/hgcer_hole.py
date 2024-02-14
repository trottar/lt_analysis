#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-13 22:57:28 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TCutG

def apply_HGCer_hole_cut(Q2, W, EPSSET, simc=False):

    #'''
    # Defined HGCer Geometric cuts
    cutg = TCutG("cutg",21)
    cutg.SetVarX("P_hgcer_yAtCer")
    cutg.SetVarY("P_hgcer_xAtCer")
    cutg.SetPoint(0,-25,2)
    cutg.SetPoint(1,-2,2)
    cutg.SetPoint(2,-1,2.5)
    cutg.SetPoint(3,0,3)
    cutg.SetPoint(4,1,3)
    cutg.SetPoint(5,2,3.3)
    cutg.SetPoint(6,3,3.0)
    cutg.SetPoint(7,4,2.5)
    cutg.SetPoint(8,5,2)
    cutg.SetPoint(9,25,2)
    cutg.SetPoint(10,25,0.5)
    cutg.SetPoint(11,5,0.5)
    cutg.SetPoint(12,4,1)
    cutg.SetPoint(13,3,-1)
    cutg.SetPoint(14,2,-2)
    cutg.SetPoint(15,1,-2.3)
    cutg.SetPoint(16,0,-1.5)
    cutg.SetPoint(17,-1,-1)
    cutg.SetPoint(18,-2,0.5)
    cutg.SetPoint(19,-25,0.5)
    cutg.SetPoint(20,-25,2)

    '''
    # Defined HGCer Geometric cuts
    cutg = TCutG("cutg")
    cutg.SetVarX("P_hgcer_yAtCer")
    cutg.SetVarY("P_hgcer_xAtCer")
    '''
    return cutg
