#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-02 06:05:13 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import numpy as np
import sys, os, math

BINNING_PATH = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "binning"))
if BINNING_PATH not in sys.path:
    sys.path.append(BINNING_PATH)

from theta_cm import calculate_tmin

###############################################################################################################################################
# Establish global variables so that they're not called each iteration of loop

# First, define empty variables
W = ""
Q2 = ""
EPSSET = ""
POL = ""
tmin = ""
tmax = ""
cut_poly = []
c0_dict = {}
cut_poly_xmin = None
cut_poly_xmax = None
cut_poly_ymin = None
cut_poly_ymax = None
c0_value = 0.0
#TMIN_RESOLUTION_THRESHOLD = -1e-6 # 0.001 MeV^2, adjust as needed based on resolution studies
TMIN_RESOLUTION_THRESHOLD = 0.0


def _sort_ccw_points(points):
    pts = [(float(p[0]), float(p[1])) for p in points]
    if len(pts) < 3:
        return pts
    cx = sum(p[0] for p in pts) / float(len(pts))
    cy = sum(p[1] for p in pts) / float(len(pts))
    pts = sorted(pts, key=lambda p: math.atan2(p[1] - cy, p[0] - cx))
    return pts


def _point_in_convex_poly(x, y, poly, eps=1e-12):
    if poly is None or len(poly) < 3:
        return False
    px, py = float(x), float(y)
    if (
        cut_poly_xmin is not None
        and (px < cut_poly_xmin or px > cut_poly_xmax or py < cut_poly_ymin or py > cut_poly_ymax)
    ):
        return False
    for i in range(len(poly)):
        a = poly[i]
        b = poly[(i + 1) % len(poly)]
        if ((b[0] - a[0]) * (py - a[1]) - (b[1] - a[1]) * (px - a[0])) < -eps:
            return False
    return True

# Then, set global variables which is called with arguments
def set_val(inpDict):
    
    global W, Q2, EPSSET, POL, ParticleType
    global tmin, tmax
    global cut_poly
    global cut_poly_xmin, cut_poly_xmax, cut_poly_ymin, cut_poly_ymax
    global c0_value
    
    W = inpDict["W"] 
    Q2 = inpDict["Q2"] 
    EPSSET = inpDict["EPSSET"]
    POL = inpDict["POL"]
    ParticleType = inpDict["ParticleType"]

    tmin = float(inpDict["tmin"] )
    tmax = float(inpDict["tmax"] )

    if inpDict.get("cut_mode") != "poly":
        print("ERROR: Invalid cut_mode. Expected 'poly'.")
        sys.exit(2)

    points = inpDict.get("poly_points", [])
    if len(points) < 3:
        print("ERROR: Invalid polygon cut. Expected at least 3 vertices.")
        sys.exit(2)
    cut_poly = _sort_ccw_points(points)
    cut_poly_xmin = min(p[0] for p in cut_poly)
    cut_poly_xmax = max(p[0] for p in cut_poly)
    cut_poly_ymin = min(p[1] for p in cut_poly)
    cut_poly_ymax = max(p[1] for p in cut_poly)
    
    ##############
    # HARD CODED #
    ##############

    global c0_dict
    
    # Adjusted HMS delta to fix hsxfp correlation
    # See Dave Gaskell's slides for more info: https://redmine.jlab.org/attachments/2316
    # Note: these momenta are from Dave's slides and may not reflect what is used here
    h_momentum_list = [0.889, 0.968, 2.185, 2.328, 3.266, 4.2, 4.712, 5.292, 6.59]
    c0_list = [-1.0, -2.0, -2.0, -2.0, -3.0, -5.0, -6.0, -6.0, -3.0]

    if ParticleType == "kaon":
        for c0, p in zip(c0_list, h_momentum_list):
            if p == 0.889:
                c0_dict["Q2p1W2p95_lowe"] = c0 # Proper value 0.888
            elif p == 0.968:
                c0_dict["Q0p5W2p40_lowe"] = c0
                c0_dict["Q3p0W3p14_lowe"] = c0 # Proper value 1.821
                c0_dict["Q5p5W3p02_lowe"] = c0 # Proper value 0.962
            elif p == 2.185:
                c0_dict["Q0p5W2p40_highe"] = c0 # Proper value 2.066
                c0_dict["Q3p0W2p32_lowe"] = c0
            elif p == 2.328:
                c0_dict["Q4p4W2p74_lowe"] = c0
            elif p == 3.266:
                c0_dict["Q5p5W3p02_highe"] = c0            
            elif p == 4.2:
                c0_dict["Q3p0W3p14_highe"] = c0 # Proper value 4.204
            elif p == 4.712:
                c0_dict["Q4p4W2p74_highe"] = c0            
            elif p == 5.292:
                c0_dict["Q2p1W2p95_highe"] = c0
            elif p == 6.59:
                c0_dict["Q3p0W2p32_highe"] = c0
    else:
        c0_dict["Q0p4W2p20_lowe"] = 0.0
        c0_dict["Q0p4W2p20_highe"] = 0.0

    try:
        c0_value = c0_dict["Q{}W{}_{}e".format(Q2, W, EPSSET)]
    except KeyError:
        c0_value = 0.0
            
    ##############
    ##############        
    ##############
    
###############################################################################################################################################

def _compute_base_data_cut_state(evt):

    ##############
    # HARD CODED #
    ##############

    adj_hsdelta = evt.hsdelta + c0_value * evt.hsxpfp

    ##############
    ##############        
    ##############

    #CUTs Definations 
    if evt.P_hod_goodstarttime != 1 or evt.P_dc_InsideDipoleExit != 1:
        return False, adj_hsdelta
    if evt.ssdelta < -10.0 or evt.ssdelta > 20.0 or evt.ssxptar < -0.06 or evt.ssxptar > 0.06 or evt.ssyptar < -0.04 or evt.ssyptar > 0.04:
        return False, adj_hsdelta

    if evt.H_hod_goodstarttime != 1 or evt.H_dc_InsideDipoleExit != 1:
        return False, adj_hsdelta
    if adj_hsdelta < -8.0 or adj_hsdelta > 8.0 or evt.hsxptar < -0.08 or evt.hsxptar > 0.08 or evt.hsyptar < -0.045 or evt.hsyptar > 0.045:
        return False, adj_hsdelta

    return _point_in_convex_poly(evt.Q2, evt.W, cut_poly), adj_hsdelta


def _in_window(value, lower, upper):
    return (lower <= value) and (value < upper)


def _passes_tmin_resolution(minus_t, w, q2):
    minus_tmin = calculate_tmin(ParticleType, POL, float(w), float(q2))
    if not (math.isfinite(minus_t) and math.isfinite(minus_tmin)):
        return False
    return (minus_t - minus_tmin) > TMIN_RESOLUTION_THRESHOLD


def get_shifted_mm(evt, mm_offset=0.0):
    try:
        return evt.MM_shift
    except AttributeError:
        return evt.MM + mm_offset


def _compute_data_cut_state(evt, mm_min=0.7, mm_max=1.5, mm_offset=0.0):
    base_cuts, adj_hsdelta = _compute_base_data_cut_state(evt)
    if not base_cuts:
        return False, False, adj_hsdelta

    adj_t = get_shifted_t(evt)
    adj_mm = get_shifted_mm(evt, mm_offset=mm_offset)
    t_in_range = _in_window(adj_t, tmin, tmax)
    tmin_resolved = _passes_tmin_resolution(adj_t, evt.W, evt.Q2)
    mm_in_range = _in_window(adj_mm, mm_min, mm_max)
    return base_cuts and t_in_range and tmin_resolved and mm_in_range, base_cuts and t_in_range and tmin_resolved, adj_hsdelta


def evaluate_data_event(evt, mm_min=0.7, mm_max=1.5, mm_offset=0.0):
    return _compute_data_cut_state(evt, mm_min, mm_max, mm_offset)


def evaluate_data_cut_bools(evt, mm_min=0.7, mm_max=1.5, mm_offset=0.0):
    allcuts, subcuts, _ = _compute_data_cut_state(evt, mm_min, mm_max, mm_offset)
    return allcuts, subcuts


def apply_data_cuts(evt, mm_min=0.7, mm_max=1.5, mm_offset=0.0):
    return _compute_data_cut_state(evt, mm_min, mm_max, mm_offset)[0]

###############################################################################################################################################

# Prefer the stored shifted branch when available; otherwise fall back to the
# legacy convention of deriving -t from MandelT on the fly.
def get_shifted_t(evt):
    try:
        return evt.t_shift
    except AttributeError:
        return -evt.MandelT

# Subtraction cuts
def apply_data_sub_cuts(evt, mm_min=0.7, mm_max=1.5, mm_offset=0.0):
    return _compute_data_cut_state(evt, mm_min, mm_max, mm_offset)[1]

################################################################################################################################################

def apply_simc_cuts(evt, mm_min=0.7, mm_max=1.5):

    ##############
    # HARD CODED #
    ##############

    adj_missmass = evt.missmass

    ##############
    ##############        
    ##############        
    
    # Define the acceptance cuts  
    SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssxptar>=-0.06) & (evt.ssxptar<=0.06) & (evt.ssyptar>=-0.04) & (evt.ssyptar<=0.04)
    HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsxptar>=-0.08) & (evt.hsxptar<=0.08) & (evt.hsyptar>=-0.045) & (evt.hsyptar<=0.045)

    ## testing...
    #SHMS_Acceptance = (evt.ssdelta>=-10.0) & (evt.ssdelta<=20.0) & (evt.ssyptar>=-0.06) & (evt.ssyptar<=0.06) & (evt.ssxptar>=-0.04) & (evt.ssxptar<=0.04)
    #HMS_Acceptance = (evt.hsdelta>=-8.0) & (evt.hsdelta<=8.0) & (evt.hsyptar>=-0.08) & (evt.hsyptar<=0.08) & (evt.hsxptar>=-0.045) & (evt.hsxptar<=0.045)

    Diamond = _point_in_convex_poly(evt.Q2, evt.W, cut_poly)

    minus_t = -evt.t
    t_in_range = _in_window(minus_t, tmin, tmax)
    tmin_resolved = _passes_tmin_resolution(minus_t, evt.W, evt.Q2)
    mm_in_range = _in_window(adj_missmass, mm_min, mm_max)
    ALLCUTS = HMS_Acceptance and SHMS_Acceptance and Diamond and t_in_range and tmin_resolved and mm_in_range
    #ALLCUTS = (-100000<=-evt.t) # No cuts

    return ALLCUTS
