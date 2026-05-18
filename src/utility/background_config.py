#! /usr/bin/python

"""
Shared analysis configuration for stage-2 background tuning and bin finding.

Edit the constants in this file when you want one place to control the
fallback BG scale, binning thresholds, and optimizer search behavior.
"""

from __future__ import annotations


BG_STAT_SCALE1 = 0.5
BG_STAT_SCALE2 = 0.5

# Background-fit shape guards. These are used to reject pathological
# empirical-fit curves before they can be subtracted or chosen by the
# optimizer.
BG_FIT_NEG_TOL = 1e-8
BG_FIT_SHAPE_SAMPLES = 128
BG_FIT_CONCAVITY_REL_TOL = 0.02
BG_FIT_INTERIOR_MIN_REL_TOL = 0.05
BG_FIT_MAX_NONPOSITIVE_FRACTION = 0.05
BG_FIT_MAX_NONPOSITIVE_SPAN = 0.006

# Candidate selection strategy for the Step-4 optimizer.
# "weighted" ranks candidates by a config-controlled weighted score.
# "lexicographic" preserves the older strict fail->mean_dev->rms->kin order.
BG_OPT_SELECTION_MODE = "weighted"

# Relative weights for the weighted selection mode. Lower is better for the
# first four metrics; higher is better for valid_ratio_bins.
BG_OPT_METRIC_WEIGHTS = {
    "kinematic_score": 0.5,
    "ratio_rms": 0.25,
    "ratio_mean_dev": 0.15,
    "ratio_fail_count": 0.05,
    "valid_ratio_bins": 0.05,
}

# Sigma cut used when ranking data/SIMC ratio agreement.
RATIO_SIGMA_THRESHOLD = 3.0

# Shared bin-finding controls migrated from find_bins.py.
PHI_BIN_MIN_EVENTS = 10
T_BIN_MIN_EVENTS = 500
MIN_PHI_BINS = 5
MIN_T_BINS = 2
PHI_BIN_MIN_DEG = -180.0
PHI_BIN_MAX_DEG = 180.0
T_BIN_ADJUST_TOLERANCE = 1e-3
T_BIN_ADJUST_MAX_ITERATIONS = 50000
T_BIN_EDGE_BIAS = 2.0

# Joint optimizer search space.
# The bash-script bin counts are assumed to be close, so only scan local
# neighbors within +/- 1 bin for the shared bin-count optimization.
NUM_T_BIN_SEARCH_OFFSETS = (-1, 0)
NUM_PHI_BIN_SEARCH_OFFSETS = (-1, 0)

BG_STAT_SCALE1_COARSE_MIN = 0.0
BG_STAT_SCALE1_COARSE_MAX = 1.0
BG_STAT_SCALE1_COARSE_STEP = 0.1
BG_STAT_SCALE1_REFINE_WINDOW = 0.1
BG_STAT_SCALE1_REFINE_STEP = 0.025

BG_STAT_SCALE2_COARSE_MIN = 0.0
BG_STAT_SCALE2_COARSE_MAX = 1.0
BG_STAT_SCALE2_COARSE_STEP = 0.1
BG_STAT_SCALE2_REFINE_WINDOW = 0.1
BG_STAT_SCALE2_REFINE_STEP = 0.025
BG_STAT_SCALE2_FINALIST_COUNT = 2

# Full-spectrum MM diagnostic plot controls used by the Step-4 optimizer PDF.
# These are check-plot only and do not change the MM cut window used elsewhere.
BG_OPT_MM_PLOT_MIN = 0.7
BG_OPT_MM_PLOT_MAX = 1.5
BG_OPT_MM_PLOT_NBINS = 100
# window: SIMC scaled to data inside MM cut window
# proper: SIMC kept at proper archived normalization
BG_OPT_MM_SIMC_SCALE_MODE = "proper"

# Variables used for lightweight SIMC-vs-data kinematic scoring.
KINEMATIC_SCORE_VARS = ("Q2", "W", "t", "mm", "xptar", "yptar")

# -------------------------------------------------------------------------
# Catalogue of background-shape models and their default side-bands
# -------------------------------------------------------------------------
BG_MODELS = {
    "Q3p0W2p32" : {

            "Fit 1" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_lowe": {
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },
            },

            "Fit 2" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.22), #Q4p4W2p74
                        #"left":  (1.06, 1.11),
                        #"right": (1.14, 1.20),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.23),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.22), #Q4p4W2p74
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.20),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },
        },

        # --- Sigma peak 2nd-order Chebyshev ---------------------------
        "sigma_peak": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                #"left":  (1.00, 1.06),
                #"right": (1.14, 1.20),
                "left":  (1.165, 1.20),
                "right": (1.20, 1.21),
            }
        },
    },
    "Q4p4W2p74" : {

            "Fit 1" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_lowe": {
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.15, 1.20),
                        "right": (1.30, 1.40),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.15, 1.20),
                        "right": (1.30, 1.40),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.15, 1.20),
                        "right": (1.30, 1.40),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.15, 1.20),
                        "right": (1.30, 1.40),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.15, 1.20),
                        "right": (1.30, 1.40),
                    }
                },
            },

            "Fit 2" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
 #                       "left":  (1.05, 1.10), #Q4p4W2p74
#                        "right": (1.30, 1.35), #Q4p4W2p74
                        "left":  (1.05, 1.10),                       
                        "right": (1.10, 1.20),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
 #                       "left":  (1.05, 1.10), #Q4p4W2p74
#                        "right": (1.30, 1.35), #Q4p4W2p74
                        "left":  (1.05, 1.10),                       
                        "right": (1.10, 1.20),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
 #                       "left":  (1.05, 1.10), #Q4p4W2p74
#                        "right": (1.30, 1.35), #Q4p4W2p74
                        "left":  (1.05, 1.10),                       
                        "right": (1.10, 1.20),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
 #                       "left":  (1.05, 1.10), #Q4p4W2p74
#                        "right": (1.30, 1.35), #Q4p4W2p74
                        "left":  (1.05, 1.10),                       
                        "right": (1.10, 1.20),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
 #                       "left":  (1.05, 1.10), #Q4p4W2p74
#                        "right": (1.30, 1.35), #Q4p4W2p74
                        "left":  (1.05, 1.10),                       
                        "right": (1.10, 1.20),
                    }
                },
        },

        # --- Sigma peak 2nd-order Chebyshev ---------------------------
        "sigma_peak": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                #"left":  (1.00, 1.06),
                #"right": (1.14, 1.20),
                "left":  (1.165, 1.20),
                "right": (1.20, 1.21),
            }
        },
    },
    "Q5p5W3p02" : {

            "Fit 1" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_lowe": {
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.25), #Q5p5W3p02
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.25), #Q5p5W3p02
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.25), #Q5p5W3p02
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.25), #Q5p5W3p02
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.25), #Q5p5W3p02
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },
            },

            "Fit 2" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.25), #Q5p5W3p02
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.22), #Q5p5W3p02
                        #"left":  (1.06, 1.11),
                        #"right": (1.14, 1.20),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.25), #Q5p5W3p02
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.23),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.22), #Q5p5W3p02
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.20),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q5p5W3p02
                        #"right": (1.20, 1.25), #Q5p5W3p02
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },
        },

        # --- Sigma peak 2nd-order Chebyshev ---------------------------
        "sigma_peak": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                #"left":  (1.00, 1.06),
                #"right": (1.14, 1.20),
                "left":  (1.165, 1.20),
                "right": (1.20, 1.21),
            }
        },
    },
    "Q2p1W2p95" : { 

            "Fit 1" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_lowe": {
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },
            },

            "Fit 2" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.22), #Q4p4W2p74
                        #"left":  (1.06, 1.11),
                        #"right": (1.14, 1.20),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.23),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.22), #Q4p4W2p74
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.20),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.05, 1.10),
                        "right": (1.12, 1.17),
                    }
                },
        },

        # --- Sigma peak 2nd-order Chebyshev ---------------------------
        "sigma_peak": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                #"left":  (1.00, 1.06),
                #"right": (1.14, 1.20),
                "left":  (1.165, 1.20),
                "right": (1.20, 1.21),
            }
        },
    },
    "Q3p0W3p14" : {

            "Fit 1" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_lowe": {
                    # quadratic, forced to be 0 at x = 1.12
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    # quadratic, forced to be 0 at x = 1.15
                    "func_expr": "[0]*(x-1.12) + [1]*(x-1.12)*(x-1.12)",
                    "n_par":      2,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.12, 1.15),
                        "right": (1.16, 1.21),
                    }
                },
            },

            "Fit 2" : {
                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.05, 1.10),
                        "right": (1.20, 1.25),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.22), #Q4p4W2p74
                        #"left":  (1.06, 1.11),
                        #"right": (1.14, 1.20),
                        "left":  (1.05, 1.10),
                        "right": (1.20, 1.25),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.23),
                        "left":  (1.05, 1.10),
                        "right": (1.20, 1.25),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.05, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.22), #Q4p4W2p74
                        #"left":  (1.08, 1.10),
                        #"right": (1.18, 1.20),
                        "left":  (1.05, 1.10),
                        "right": (1.20, 1.25),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "cheb2_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.05, 1.10),
                        "right": (1.20, 1.25),
                    }
                },
        },

        # --- Sigma peak 2nd-order Chebyshev ---------------------------
        "sigma_peak": {
            "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
            "n_par":      3,
            "sidebands": {
                #"left":  (1.00, 1.06),
                #"right": (1.14, 1.20),
                "left":  (1.165, 1.20),
                "right": (1.20, 1.21),
            }
        },
    },
}

################################################################################################################################################


def _round_scale(value):
    return round(float(value), 6)


def get_bg_scale_setting_key(epsset, phi_setting):
    return "{}:{}".format(str(epsset).lower(), str(phi_setting))


def resolve_bg_stat_scale1(inpDict, phi_setting=None):
    if phi_setting is None:
        phi_setting = inpDict.get("phi_setting")
    key = get_bg_scale_setting_key(inpDict.get("EPSSET", ""), phi_setting)
    scale_map = inpDict.get("bg_stat_scale1_by_setting", {})
    if key in scale_map:
        return float(scale_map[key])
    if "bg_stat_scale1" in inpDict:
        return float(inpDict["bg_stat_scale1"])
    return float(BG_STAT_SCALE1)


def resolve_bg_stat_scale2(inpDict, phi_setting=None):
    if phi_setting is None:
        phi_setting = inpDict.get("phi_setting")
    key = get_bg_scale_setting_key(inpDict.get("EPSSET", ""), phi_setting)
    scale_map = inpDict.get("bg_stat_scale2_by_setting", {})
    if key in scale_map:
        return float(scale_map[key])
    if "bg_stat_scale2" in inpDict:
        return float(inpDict["bg_stat_scale2"])
    return float(BG_STAT_SCALE2)


def build_bin_count_candidates(requested_t_bins, requested_phi_bins):
    candidates = []
    seen = set()
    for t_offset in NUM_T_BIN_SEARCH_OFFSETS:
        for phi_offset in NUM_PHI_BIN_SEARCH_OFFSETS:
            t_bins = max(MIN_T_BINS, int(requested_t_bins) + int(t_offset))
            phi_bins = max(MIN_PHI_BINS, int(requested_phi_bins) + int(phi_offset))
            key = (t_bins, phi_bins)
            if key not in seen:
                seen.add(key)
                candidates.append(key)
    return candidates


def _build_scale_coarse_candidates(default_value, coarse_min, coarse_max, coarse_step):
    values = []
    candidate = coarse_min
    while candidate <= coarse_max + 1e-12:
        values.append(_round_scale(candidate))
        candidate += coarse_step
    if _round_scale(default_value) not in values:
        values.append(_round_scale(default_value))
    return sorted(set(values))


def _build_scale_refined_candidates(best_value, coarse_min, coarse_max, refine_window, refine_step):
    lo = max(coarse_min, float(best_value) - refine_window)
    hi = min(coarse_max, float(best_value) + refine_window)
    values = []
    candidate = lo
    while candidate <= hi + 1e-12:
        values.append(_round_scale(candidate))
        candidate += refine_step
    values.append(_round_scale(best_value))
    return sorted(set(values))


def get_bg_scale1_coarse_candidates():
    return _build_scale_coarse_candidates(
        BG_STAT_SCALE1,
        BG_STAT_SCALE1_COARSE_MIN,
        BG_STAT_SCALE1_COARSE_MAX,
        BG_STAT_SCALE1_COARSE_STEP,
    )


def get_bg_scale1_refined_candidates(best_value):
    return _build_scale_refined_candidates(
        best_value,
        BG_STAT_SCALE1_COARSE_MIN,
        BG_STAT_SCALE1_COARSE_MAX,
        BG_STAT_SCALE1_REFINE_WINDOW,
        BG_STAT_SCALE1_REFINE_STEP,
    )


def get_bg_scale_coarse_candidates():
    return _build_scale_coarse_candidates(
        BG_STAT_SCALE2,
        BG_STAT_SCALE2_COARSE_MIN,
        BG_STAT_SCALE2_COARSE_MAX,
        BG_STAT_SCALE2_COARSE_STEP,
    )


def get_bg_scale_refined_candidates(best_value):
    return _build_scale_refined_candidates(
        best_value,
        BG_STAT_SCALE2_COARSE_MIN,
        BG_STAT_SCALE2_COARSE_MAX,
        BG_STAT_SCALE2_REFINE_WINDOW,
        BG_STAT_SCALE2_REFINE_STEP,
    )
