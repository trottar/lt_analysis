#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-13 22:42:24 trottar"
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
    cutg.SetPoint(0,-25+10,2)
    cutg.SetPoint(1,-2+10,2)
    cutg.SetPoint(2,-1+10,2.5)
    cutg.SetPoint(3,0+10,3)
    cutg.SetPoint(4,1+10,3)
    cutg.SetPoint(5,2+10,3.3)
    cutg.SetPoint(6,3+10,3.0)
    cutg.SetPoint(7,4+10,2.5)
    cutg.SetPoint(8,5+10,2)
    cutg.SetPoint(9,25+10,2)
    cutg.SetPoint(10,25+10,0.5)
    cutg.SetPoint(11,5+10,0.5)
    cutg.SetPoint(12,4+10,1)
    cutg.SetPoint(13,3+10,-1)
    cutg.SetPoint(14,2+10,-2)
    cutg.SetPoint(15,1+10,-2.3)
    cutg.SetPoint(16,0+10,-1.5)
    cutg.SetPoint(17,-1+10,-1)
    cutg.SetPoint(18,-2+10,0.5)
    cutg.SetPoint(19,-25+10,0.5)
    cutg.SetPoint(20,-25+10,2)

    '''
    # Defined HGCer Geometric cuts
    cutg = TCutG("cutg")
    cutg.SetVarX("P_hgcer_yAtCer")
    cutg.SetVarY("P_hgcer_xAtCer")
    '''
    return cutg
