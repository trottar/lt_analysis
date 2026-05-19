#! /usr/bin/python

"""
Shared analysis configuration for stage-2 background tuning and bin finding.

Edit the constants in this file when you want one place to control the
fallback BG scale, binning thresholds, and optimizer search behavior.
"""

from __future__ import annotations

import json
import os
from copy import deepcopy


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

BG_OPT_ACTIVE_PROFILE = "nominal_weighted"

BG_OPT_PROFILES = {
    "nominal_weighted": {
        "selection_mode": "weighted",
        "metric_weights": deepcopy(BG_OPT_METRIC_WEIGHTS),
    },
    "ratio_first": {
        "selection_mode": "lexicographic",
    },
    "kinematic_diagnostic_only": {
        "selection_mode": "weighted",
        "metric_weights": {
            "kinematic_score": 0.00,
            "ratio_rms": 0.40,
            "ratio_mean_dev": 0.35,
            "ratio_fail_count": 0.20,
            "valid_ratio_bins": 0.05,
        },
    },
    "no_empirical_residual": {
        "force_bg_stat_scale1": 0.0,
        "force_bg_stat_scale2": 0.0,
    },
    "fit1_only": {
        "force_bg_stat_scale2": 0.0,
    },
    "fit2_only": {
        "force_bg_stat_scale1": 0.0,
    },
    "high_uses_low_epsilon_scales": {
        "selection_mode": "lexicographic",
        "common_epsilon_scale_behavior": "reuse_low_epsilon_scales_for_high",
    },
}

BG_SYSTEMATIC_PROFILES = [
    "nominal_weighted",
    "ratio_first",
    "kinematic_diagnostic_only",
    "no_empirical_residual",
    "fit1_only",
    "fit2_only",
    "high_uses_low_epsilon_scales",
]

# Sigma cut used when ranking data/SIMC ratio agreement.
RATIO_SIGMA_THRESHOLD = 3.0

ALLOW_CONFIG_DRIFT = False
EMP_EPS_DIFF_WARN_THRESHOLD = 0.03
EMP_EPS_DIFF_FAIL_THRESHOLD = None
BG_OVERSUB_WARN_FRACTION = 0.02
BG_OVERSUB_WARN_MAX_RATIO = 1.10
BG_SYSTEMATICS_USE_ISOLATED_WORKSPACE = False
REGRESSION_ABS_TOL = 1e-9
REGRESSION_REL_TOL = 1e-6

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

_PROFILE_OVERRIDE_ENV = "LT_BG_PROFILE_SNAPSHOT_JSON"

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


def _load_runtime_overrides():
    override_path = os.environ.get(_PROFILE_OVERRIDE_ENV, "").strip()
    if not override_path:
        return {}
    if not os.path.exists(override_path):
        raise FileNotFoundError(
            "Background-config override file from {} was not found: {}".format(
                _PROFILE_OVERRIDE_ENV,
                override_path,
            )
        )
    with open(override_path, "r") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(
            "Background-config override file {} must contain a JSON object".format(
                override_path
            )
        )
    return payload


_RUNTIME_BG_CONFIG_OVERRIDES = _load_runtime_overrides()
for _override_key, _override_value in _RUNTIME_BG_CONFIG_OVERRIDES.items():
    if isinstance(_override_key, str) and _override_key.isupper():
        globals()[_override_key] = _override_value


def get_available_bg_profile_names():
    return sorted(BG_OPT_PROFILES.keys())


def get_active_bg_profile_name():
    return str(BG_OPT_ACTIVE_PROFILE).strip()


def _normalize_selection_mode(value):
    mode = str(value).strip().lower()
    if mode not in {"weighted", "lexicographic"}:
        raise ValueError(
            "Invalid background selection mode '{}' in profile '{}'. "
            "Allowed values are 'weighted' or 'lexicographic'.".format(
                value,
                get_active_bg_profile_name(),
            )
        )
    return mode


def _normalize_metric_weights(weights):
    if weights is None:
        return deepcopy(BG_OPT_METRIC_WEIGHTS)
    if not isinstance(weights, dict):
        raise ValueError(
            "metric_weights for background profile '{}' must be a dictionary".format(
                get_active_bg_profile_name()
            )
        )
    normalized = deepcopy(BG_OPT_METRIC_WEIGHTS)
    for key, value in weights.items():
        normalized[str(key)] = float(value)
    return normalized


def _normalize_common_epsilon_scale_behavior(raw_profile, profile_name):
    if profile_name == "common_epsilon_scales":
        raise ValueError(
            "Background profile 'common_epsilon_scales' was renamed to "
            "'high_uses_low_epsilon_scales'. The current implementation reuses "
            "the frozen low-epsilon scale map on the high-epsilon pass rather "
            "than performing a true joint low/high optimization."
        )

    behavior = raw_profile.get("common_epsilon_scale_behavior")
    if behavior is None:
        if raw_profile.get("use_low_epsilon_scales_for_high"):
            behavior = "reuse_low_epsilon_scales_for_high"
        elif raw_profile.get("use_common_epsilon_scales"):
            behavior = "reuse_low_epsilon_scales_for_high"
        else:
            behavior = "independent"

    behavior = str(behavior).strip().lower()
    allowed = {"independent", "reuse_low_epsilon_scales_for_high"}
    if behavior not in allowed:
        raise ValueError(
            "Invalid common_epsilon_scale_behavior '{}' in profile '{}'. "
            "Allowed values are {}.".format(
                behavior,
                profile_name,
                ", ".join(sorted(allowed)),
            )
        )
    return behavior


def get_active_bg_profile():
    profile_name = get_active_bg_profile_name()
    if profile_name == "common_epsilon_scales":
        raise ValueError(
            "Background profile 'common_epsilon_scales' was renamed to "
            "'high_uses_low_epsilon_scales'. The current implementation reuses "
            "the frozen low-epsilon scale map on the high-epsilon pass rather "
            "than performing a true joint low/high optimization."
        )
    if profile_name not in BG_OPT_PROFILES:
        raise ValueError(
            "Invalid BG_OPT_ACTIVE_PROFILE '{}'. Allowed profiles: {}".format(
                profile_name,
                ", ".join(get_available_bg_profile_names()),
            )
        )
    raw_profile = deepcopy(BG_OPT_PROFILES[profile_name])
    resolved = {
        "name": profile_name,
        "selection_mode": _normalize_selection_mode(
            raw_profile.get("selection_mode", BG_OPT_SELECTION_MODE)
        ),
        "metric_weights": _normalize_metric_weights(
            raw_profile.get("metric_weights", BG_OPT_METRIC_WEIGHTS)
        ),
        "force_bg_stat_scale1": raw_profile.get("force_bg_stat_scale1"),
        "force_bg_stat_scale2": raw_profile.get("force_bg_stat_scale2"),
        "common_epsilon_scale_behavior": _normalize_common_epsilon_scale_behavior(
            raw_profile,
            profile_name,
        ),
        "constraints": deepcopy(raw_profile.get("constraints", {})),
        "raw_profile": raw_profile,
    }
    resolved["use_common_epsilon_scales"] = bool(
        resolved["common_epsilon_scale_behavior"] != "independent"
    )
    return resolved


def get_active_selection_mode():
    return get_active_bg_profile()["selection_mode"]


def get_active_metric_weights():
    return deepcopy(get_active_bg_profile()["metric_weights"])


def get_forced_bg_scale1():
    value = get_active_bg_profile().get("force_bg_stat_scale1")
    return None if value is None else float(value)


def get_forced_bg_scale2():
    value = get_active_bg_profile().get("force_bg_stat_scale2")
    return None if value is None else float(value)


def use_common_epsilon_scales():
    return bool(get_active_bg_profile().get("use_common_epsilon_scales", False))


def get_common_epsilon_scale_behavior():
    return str(get_active_bg_profile().get("common_epsilon_scale_behavior", "independent"))


def get_resolved_bg_profile_settings():
    return deepcopy(get_active_bg_profile())


def get_runtime_bg_config_overrides():
    return deepcopy(_RUNTIME_BG_CONFIG_OVERRIDES)


_ACTIVE_BG_PROFILE_VALIDATION = get_active_bg_profile()


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
