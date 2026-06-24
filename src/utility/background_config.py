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


#BG_STAT_SCALE1 = 0.5
#BG_STAT_SCALE2 = 0.5
BG_STAT_SCALE1 = 0.0
BG_STAT_SCALE2 = 0.0

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

#BG_OPT_ACTIVE_PROFILE = "nominal_weighted"
BG_OPT_ACTIVE_PROFILE = "no_empirical_residual"

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
#BG_OPT_MM_SIMC_SCALE_MODE = "window"

PARTICLE_SUBTRACTION_MODE_SINGLE_SCALE = "single_scale"
PARTICLE_SUBTRACTION_MODE_COMPONENTS = "simc_shape_components"
PARTICLE_SUBTRACTION_MODE_DEFAULT = PARTICLE_SUBTRACTION_MODE_COMPONENTS
#PARTICLE_SUBTRACTION_MODE_DEFAULT = PARTICLE_SUBTRACTION_MODE_SINGLE_SCALE
PARTICLE_SUBTRACTION_MODES = (
    PARTICLE_SUBTRACTION_MODE_SINGLE_SCALE,
    PARTICLE_SUBTRACTION_MODE_COMPONENTS,
)
PARTICLE_SUBTRACTION_FALLBACK_MODE_DEFAULT = "error"
#PARTICLE_SUBTRACTION_FALLBACK_MODE_DEFAULT = PARTICLE_SUBTRACTION_MODE_SINGLE_SCALE
PARTICLE_SUBTRACTION_FALLBACK_MODES = (
    "error",
    PARTICLE_SUBTRACTION_MODE_SINGLE_SCALE,
    "zero",
    "skip_bin",
)
PARTICLE_SUBTRACTION_WEIGHT_DENOM_FLOOR_DEFAULT = 1e-12
PARTICLE_SUBTRACTION_WEIGHT_CLIP_MIN_DEFAULT = 0.0
PARTICLE_SUBTRACTION_WEIGHT_CLIP_MAX_DEFAULT = None
PARTICLE_SUBTRACTION_WEIGHT_WARN_MAX_DEFAULT = 10.0
SIMC_TREE_NAME_DEFAULT = "h10"
SIMC_PION_COMPONENT_BACKGROUND_MAP = {
    "pi_n": "neutron",
    "pi_delta": "delta",
    "pi_sidis": "sidis",
}

# Separate normalization windows used when scaling the pion-selected
# subtraction sample in kaon analyses.
PARTICLE_SUBTRACTION_WINDOW_CONFIG = {
    "kaon": {
        "pion": {
            "apply_mm_offset_data": True,
            "enabled_windows": {
                "pi_n": True,
                "pi_delta": False,
            },
            "windows": {
                "pi_n": (0.88, 0.94),
                "pi_delta": (1.35, 1.40),
            },
        },
    },
}

PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG = {
    "pion_control": {
        "apply_mm_offset_data": True,
        "staged_fit_passes": 1,
        "fit_order": ("pi_n", "pi_sidis", "pi_delta", "k_sigma0_signal"),
        "stage_amplitude_windows": {
            "pi_delta": (1.18, 1.23),
            "k_sigma0_signal": (1.18, 1.23),
        },
        "stage_amplitude_modes": {
            "pi_delta": "window_integral",
            "k_sigma0_signal": "window_integral",
        },
        "prior_scales": {
            "pi_n": 1.0,
            "pi_delta": 1.5,
            "pi_sidis": 2.0,
            "k_sigma0_signal": 1.0,
        },
        # Keep pion_control post-fit scales fixed at unity.
        # This fit defines the control-model denominator in w_pi(MM), so only
        # the kaon_nosub side should be manually tuned with post-fit scales.
        "postfit_component_scales": {
            "pi_n": 1.0,
            "pi_delta": 1.0,
            "pi_sidis": 1.0,
            "k_sigma0_signal": 1.0,
        },
        "joint_refinement_enabled": True,
        "particle_subtraction_max_fit_cycles": 50,
        "particle_subtraction_fit_tolerance": 1e-5,
        "oversub_sigma_tolerance": 2.0,
        "max_oversub_bin_count": None,
        "max_oversub_bin_fraction": None,
        #"max_oversub_bin_count": 8,
        #"max_oversub_bin_fraction": 0.20,
        "max_full_range_chi2_ndf": None,
        "enabled_windows": {
            "pi_n": True,
            "pi_delta": True,
            "pi_sidis": True,
            #"k_sigma0_signal": True,
            "k_sigma0_signal": False,
        },
        "windows": {
            "pi_n": (0.90, 0.95),
            "pi_delta": (1.18, 1.23),
            "pi_sidis": ((1.07, 1.10),(1.45, 1.50)),
            #"pi_sidis": ((1.05, 1.10), (1.25, 1.30), (1.45, 1.50)),
            "k_sigma0_signal": (1.17, 1.23),
        },
    },
    "kaon_nosub": {
        "apply_mm_offset_data": True,
        "staged_fit_passes": 1,
        "fit_order": ("pi_n", "pi_sidis", "pi_delta", "k_sigma0_signal"),
        "stage_amplitude_windows": {
            "pi_delta": (1.18, 1.23),
            "k_sigma0_signal": (1.17, 1.23),
        },
        "stage_amplitude_modes": {
            "pi_delta": "window_integral",
            "k_sigma0_signal": "window_integral",
        },
        "prior_scales": {
            "pi_n": 1.0,
            "pi_delta": 1.5,
            "pi_sidis": 2.0,
            "k_sigma0_signal": 1.0,
        },
        # Kaon-side post-fit scales may be tuned phenomenologically after the
        # raw component fit if the unscaled kaon pion-background model is too
        # large or shape-misaligned relative to the kaon no-sub spectrum.
        "postfit_component_scales": {
            "pi_n": 0.95,
            "pi_delta": 0.40,
            "pi_sidis": 0.35,
            "k_sigma0_signal": 1.0,
        },
        "include_kaon_signal_template": False,
        "joint_refinement_enabled": True,
        "particle_subtraction_prior_scale_k_lambda_signal": 1.0,
        "particle_subtraction_max_fit_cycles": 50,
        "particle_subtraction_fit_tolerance": 1e-5,
        "oversub_sigma_tolerance": 2.0,
        "max_oversub_bin_count": None,
        "max_oversub_bin_fraction": None,
        #"max_oversub_bin_count": 8,
        #"max_oversub_bin_fraction": 0.20,
        "max_full_range_chi2_ndf": None,
        "kaon_signal_tail_extension": 0.02,
        "enabled_windows": {
            "pi_n": True,
            "pi_delta": True,
            "pi_sidis": True,
            #"k_sigma0_signal": True,
            "k_sigma0_signal": False,
        },
        "windows": {
            "pi_n": (0.90, 0.95),
            "pi_delta": (1.18, 1.23),
            "pi_sidis": ((1.07, 1.10),(1.45, 1.50)),
            #"pi_sidis": ((1.05, 1.10), (1.25, 1.30), (1.45, 1.50)),
            "k_sigma0_signal": (1.17, 1.23),
        },
        "enabled_excluded_windows": {
            "sigma_peak": False,
        },
        "excluded_windows": {
            "sigma_peak": (1.17, 1.23),
        },
    },
}

PARTICLE_SUBTRACTION_WINDOW_CONFIG_OVERRIDES = {}

PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG_OVERRIDES = {}

PARTICLE_SUBTRACTION_CONFIG_MERGE_KEYS = frozenset(
    {
        "enabled_windows",
        "windows",
        "enabled_excluded_windows",
        "excluded_windows",
        "stage_amplitude_windows",
        "stage_amplitude_modes",
        "prior_scales",
        "postfit_component_scales",
    }
)

PARTICLE_SUBTRACTION_COMPONENT_PRIOR_SCALE_LEGACY_KEYS = {
    "pi_n": "particle_subtraction_prior_scale_pi_n",
    "pi_delta": "particle_subtraction_prior_scale_pi_delta",
    "pi_sidis": "particle_subtraction_prior_scale_pi_sidis",
    "k_lambda_signal": "particle_subtraction_prior_scale_k_lambda_signal",
    "k_sigma0_signal": "particle_subtraction_prior_scale_k_sigma0_signal",
}

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
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_lowe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Center_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Left_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
                    }
                },

                    # --- ** GOOD ** 2nd-order Chebyshev ---------------------------
                "fixquad_Right_highe": {
                    "func_expr": "cheb2", # 2nd-order Chebyshev polynomial
                    "n_par":      3,
                    "sidebands": {
                        #"left":  (1.08, 1.10), #Q4p4W2p74
                        #"right": (1.20, 1.25), #Q4p4W2p74
                        #"left":  (1.06, 1.10),
                        #"right": (1.18, 1.21),
                        "left":  (1.10, 1.15),
                        "right": (1.25, 1.30),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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
                        "right": (1.15, 1.20),
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


def _normalize_particle_subtraction_mode(value):
    mode = str(value or PARTICLE_SUBTRACTION_MODE_DEFAULT).strip().lower()
    if mode not in PARTICLE_SUBTRACTION_MODES:
        raise ValueError(
            "Invalid particle_subtraction_mode '{}'. Allowed values are {}.".format(
                value,
                ", ".join(PARTICLE_SUBTRACTION_MODES),
            )
        )
    return mode


def resolve_particle_subtraction_mode(inp_dict=None, mode=None):
    if mode is None and isinstance(inp_dict, dict):
        mode = inp_dict.get(
            "particle_subtraction_mode",
            PARTICLE_SUBTRACTION_MODE_DEFAULT,
        )
    elif mode is None:
        mode = PARTICLE_SUBTRACTION_MODE_DEFAULT
    return _normalize_particle_subtraction_mode(mode)


def _normalize_particle_subtraction_fallback_mode(value):
    mode = str(value or PARTICLE_SUBTRACTION_FALLBACK_MODE_DEFAULT).strip().lower()
    if mode not in PARTICLE_SUBTRACTION_FALLBACK_MODES:
        raise ValueError(
            "Invalid particle_subtraction_fallback_mode '{}'. Allowed values are {}.".format(
                value,
                ", ".join(PARTICLE_SUBTRACTION_FALLBACK_MODES),
            )
        )
    return mode


def resolve_particle_subtraction_fallback_mode(inp_dict=None, mode=None):
    if mode is None and isinstance(inp_dict, dict):
        mode = inp_dict.get(
            "particle_subtraction_fallback_mode",
            PARTICLE_SUBTRACTION_FALLBACK_MODE_DEFAULT,
        )
    elif mode is None:
        mode = PARTICLE_SUBTRACTION_FALLBACK_MODE_DEFAULT
    return _normalize_particle_subtraction_fallback_mode(mode)


def resolve_particle_subtraction_weight_denominator_floor(inp_dict=None):
    if isinstance(inp_dict, dict):
        return float(
            inp_dict.get(
                "particle_subtraction_weight_denom_floor",
                PARTICLE_SUBTRACTION_WEIGHT_DENOM_FLOOR_DEFAULT,
            )
        )
    return float(PARTICLE_SUBTRACTION_WEIGHT_DENOM_FLOOR_DEFAULT)


def resolve_particle_subtraction_weight_clip_bounds(inp_dict=None):
    clip_min = PARTICLE_SUBTRACTION_WEIGHT_CLIP_MIN_DEFAULT
    clip_max = PARTICLE_SUBTRACTION_WEIGHT_CLIP_MAX_DEFAULT
    if isinstance(inp_dict, dict):
        clip_min = inp_dict.get(
            "particle_subtraction_weight_clip_min",
            PARTICLE_SUBTRACTION_WEIGHT_CLIP_MIN_DEFAULT,
        )
        clip_max = inp_dict.get(
            "particle_subtraction_weight_clip_max",
            PARTICLE_SUBTRACTION_WEIGHT_CLIP_MAX_DEFAULT,
        )
    clip_min = None if clip_min is None else float(clip_min)
    clip_max = None if clip_max is None else float(clip_max)
    return clip_min, clip_max


def resolve_particle_subtraction_weight_warn_max(inp_dict=None):
    if isinstance(inp_dict, dict):
        return float(
            inp_dict.get(
                "particle_subtraction_weight_warn_max",
                PARTICLE_SUBTRACTION_WEIGHT_WARN_MAX_DEFAULT,
            )
        )
    return float(PARTICLE_SUBTRACTION_WEIGHT_WARN_MAX_DEFAULT)


def resolve_simc_tree_name(inp_dict=None):
    tree_name = ""
    if isinstance(inp_dict, dict):
        tree_name = str(inp_dict.get("simc_tree_name", SIMC_TREE_NAME_DEFAULT)).strip()
    return tree_name or SIMC_TREE_NAME_DEFAULT


def _build_default_simc_pion_component_file_map(background_samples):
    resolved = {phi_setting: {} for phi_setting in ("Center", "Left", "Right")}
    by_phi = {}
    if isinstance(background_samples, dict):
        by_phi = background_samples.get("by_phi", {}) or {}

    for phi_setting in resolved:
        phi_samples = by_phi.get(phi_setting, {}) or {}
        for component_name, background_name in SIMC_PION_COMPONENT_BACKGROUND_MAP.items():
            sample_entry = phi_samples.get(background_name, {}) or {}
            resolved[phi_setting][component_name] = sample_entry.get("root")

    return resolved


def resolve_simc_pion_component_files(inp_dict=None, phi_setting=None):
    raw_map = None
    background_samples = None
    if isinstance(inp_dict, dict):
        raw_map = inp_dict.get("simc_pion_component_files")
        background_samples = inp_dict.get("background_samples")

    resolved = _build_default_simc_pion_component_file_map(background_samples)
    if isinstance(raw_map, dict):
        for raw_phi_setting, raw_components in raw_map.items():
            if raw_phi_setting not in resolved:
                resolved[str(raw_phi_setting)] = {}
            if not isinstance(raw_components, dict):
                continue
            for component_name, filename in raw_components.items():
                resolved[str(raw_phi_setting)][str(component_name)] = filename

    for resolved_phi_setting, component_map in resolved.items():
        for component_name in SIMC_PION_COMPONENT_BACKGROUND_MAP:
            component_map.setdefault(component_name, None)

    if phi_setting is None:
        return deepcopy(resolved)
    return deepcopy(resolved.get(str(phi_setting), {}))


def _normalize_particle_subtraction_setting_key(q2=None, w=None):
    q2_text = str(q2 or "").strip()
    w_text = str(w or "").strip()
    if not q2_text or not w_text:
        return None
    if q2_text.lower().startswith("q"):
        q2_text = q2_text[1:]
    if w_text.lower().startswith("w"):
        w_text = w_text[1:]
    if not q2_text or not w_text:
        return None
    return "Q{}W{}".format(q2_text, w_text)


def get_particle_subtraction_setting_key(inp_dict=None, q2=None, w=None):
    if isinstance(inp_dict, dict):
        if q2 is None:
            q2 = inp_dict.get("Q2")
        if w is None:
            w = inp_dict.get("W")
    return _normalize_particle_subtraction_setting_key(q2=q2, w=w)


def _get_casefold_mapping_entry(mapping, key):
    if not isinstance(mapping, dict):
        return None
    normalized_key = str(key or "").strip().lower()
    if not normalized_key:
        return None
    for candidate_key, candidate_value in mapping.items():
        if str(candidate_key or "").strip().lower() == normalized_key:
            return candidate_value
    return None


def _deep_merge_particle_subtraction_config(base_config, override_config):
    if base_config is None:
        merged = {}
    else:
        merged = deepcopy(base_config)
    if not isinstance(override_config, dict):
        return merged
    for raw_key, raw_value in override_config.items():
        key = str(raw_key)
        if (
            key in PARTICLE_SUBTRACTION_CONFIG_MERGE_KEYS
            and isinstance(merged.get(key), dict)
            and isinstance(raw_value, dict)
        ):
            child = deepcopy(merged.get(key) or {})
            for child_key, child_value in raw_value.items():
                child[str(child_key)] = deepcopy(child_value)
            merged[key] = child
        else:
            merged[key] = deepcopy(raw_value)
    return merged


def _resolve_particle_subtraction_override_context(
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    resolved_setting_key = str(setting_key).strip() if setting_key is not None else ""
    if not resolved_setting_key:
        resolved_setting_key = get_particle_subtraction_setting_key(inp_dict)
    resolved_setting_key = resolved_setting_key or None
    resolved_phi_setting = phi_setting
    if resolved_phi_setting is None and isinstance(inp_dict, dict):
        resolved_phi_setting = inp_dict.get("phi_setting")
    resolved_phi_setting = str(resolved_phi_setting).strip() if resolved_phi_setting is not None else None
    if resolved_phi_setting == "":
        resolved_phi_setting = None
    return resolved_setting_key, resolved_phi_setting


def _resolve_particle_subtraction_override_layers(
    override_root,
    setting_key=None,
    phi_setting=None,
):
    if not isinstance(override_root, dict) or not setting_key:
        return []
    setting_entry = override_root.get(str(setting_key))
    if not isinstance(setting_entry, dict):
        return []

    layers = []
    default_override = setting_entry.get("default")
    if isinstance(default_override, dict) and default_override:
        layers.append(
            {
                "layer": "default",
                "path": "{}:default".format(setting_key),
                "payload": default_override,
            }
        )

    if phi_setting:
        phi_override = _get_casefold_mapping_entry(setting_entry.get("by_phi"), phi_setting)
        if isinstance(phi_override, dict) and phi_override:
            layers.append(
                {
                    "layer": "by_phi",
                    "path": "{}:by_phi:{}".format(setting_key, phi_setting),
                    "payload": phi_override,
                }
            )

    return layers


def _attach_particle_subtraction_resolution_metadata(
    config,
    setting_key=None,
    phi_setting=None,
    override_layers=None,
):
    if config is None:
        return None
    metadata_ready = deepcopy(config)
    metadata_ready["particle_subtraction_setting_key"] = setting_key
    metadata_ready["particle_subtraction_phi_setting"] = phi_setting
    metadata_ready["particle_subtraction_override_layers"] = list(override_layers or [])
    metadata_ready["particle_subtraction_override_applied"] = bool(override_layers)
    return metadata_ready


def get_particle_subtraction_window_config(
    particle_type,
    subtracted_particle,
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    particle_key = str(particle_type or "").strip().lower()
    subtracted_key = str(subtracted_particle or "").strip().lower()
    particle_config = PARTICLE_SUBTRACTION_WINDOW_CONFIG.get(particle_key, {})
    config = deepcopy(particle_config.get(subtracted_key)) if subtracted_key in particle_config else None

    resolved_setting_key, resolved_phi_setting = _resolve_particle_subtraction_override_context(
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        setting_key=setting_key,
    )
    override_layers = _resolve_particle_subtraction_override_layers(
        PARTICLE_SUBTRACTION_WINDOW_CONFIG_OVERRIDES,
        setting_key=resolved_setting_key,
        phi_setting=resolved_phi_setting,
    )

    for override_layer in override_layers:
        layer_payload = override_layer.get("payload")
        particle_override = _get_casefold_mapping_entry(layer_payload, particle_key)
        config_override = _get_casefold_mapping_entry(particle_override, subtracted_key)
        if isinstance(config_override, dict):
            config = _deep_merge_particle_subtraction_config(config, config_override)

    if config is None:
        return None
    return _attach_particle_subtraction_resolution_metadata(
        config,
        setting_key=resolved_setting_key,
        phi_setting=resolved_phi_setting,
        override_layers=[layer.get("path") for layer in override_layers],
    )


def resolve_particle_subtraction_windows(
    particle_type,
    subtracted_particle,
    mm_offset_data=0.0,
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    config = get_particle_subtraction_window_config(
        particle_type,
        subtracted_particle,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        setting_key=setting_key,
    )
    if not config:
        return {}

    offset = float(mm_offset_data) if bool(config.get("apply_mm_offset_data", False)) else 0.0
    enabled_windows = config.get("enabled_windows") or {}
    resolved = {}
    for window_name, bounds in (config.get("windows") or {}).items():
        if not bool(enabled_windows.get(window_name, True)):
            continue
        if len(bounds) != 2:
            raise ValueError(
                "Particle subtraction window '{}' for {} -> {} must contain exactly two bounds".format(
                    window_name,
                    particle_type,
                    subtracted_particle,
                )
            )
        low_edge, high_edge = bounds
        resolved[str(window_name)] = (float(low_edge) + offset, float(high_edge) + offset)
    return resolved


def get_particle_subtraction_component_fit_window_config(
    fit_target,
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    fit_target_key = str(fit_target or "").strip()
    config = (
        deepcopy(PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG.get(fit_target_key))
        if fit_target_key in PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG
        else None
    )

    resolved_setting_key, resolved_phi_setting = _resolve_particle_subtraction_override_context(
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        setting_key=setting_key,
    )
    override_layers = _resolve_particle_subtraction_override_layers(
        PARTICLE_SUBTRACTION_COMPONENT_FIT_WINDOW_CONFIG_OVERRIDES,
        setting_key=resolved_setting_key,
        phi_setting=resolved_phi_setting,
    )

    for override_layer in override_layers:
        config_override = _get_casefold_mapping_entry(
            override_layer.get("payload"),
            fit_target_key,
        )
        if isinstance(config_override, dict):
            config = _deep_merge_particle_subtraction_config(config, config_override)

    if config is None:
        return None
    return _attach_particle_subtraction_resolution_metadata(
        config,
        setting_key=resolved_setting_key,
        phi_setting=resolved_phi_setting,
        override_layers=[layer.get("path") for layer in override_layers],
    )


def _is_particle_subtraction_numeric_pair(bounds):
    if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
        return False
    try:
        float(bounds[0])
        float(bounds[1])
    except (TypeError, ValueError):
        return False
    return True


def _resolve_particle_subtraction_window_collection(bounds, offset, context):
    if bounds is None:
        return []

    if _is_particle_subtraction_numeric_pair(bounds):
        return [(float(bounds[0]) + offset, float(bounds[1]) + offset)]

    resolved = []
    if not isinstance(bounds, (list, tuple)):
        raise ValueError("{} must contain one window pair or an iterable of window pairs".format(context))

    for window_bounds in bounds:
        if not _is_particle_subtraction_numeric_pair(window_bounds):
            raise ValueError(
                "{} must contain only window pairs of the form (min, max)".format(context)
            )
        resolved.append(
            (float(window_bounds[0]) + offset, float(window_bounds[1]) + offset)
        )
    return resolved


def resolve_particle_subtraction_component_fit_windows(
    fit_target,
    mm_offset_data=0.0,
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    config = get_particle_subtraction_component_fit_window_config(
        fit_target,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        setting_key=setting_key,
    )
    if not config:
        return {}

    offset = float(mm_offset_data) if bool(config.get("apply_mm_offset_data", False)) else 0.0
    enabled_windows = config.get("enabled_windows") or {}
    resolved = {}
    for window_name, bounds in (config.get("windows") or {}).items():
        if not bool(enabled_windows.get(window_name, True)):
            continue
        resolved[str(window_name)] = _resolve_particle_subtraction_window_collection(
            bounds,
            offset,
            "Particle subtraction component-fit window '{}' for '{}'".format(
                window_name,
                fit_target,
            ),
        )
    return resolved


def resolve_particle_subtraction_component_fit_excluded_windows(
    fit_target,
    mm_offset_data=0.0,
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    config = get_particle_subtraction_component_fit_window_config(
        fit_target,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        setting_key=setting_key,
    )
    if not config:
        return []

    offset = float(mm_offset_data) if bool(config.get("apply_mm_offset_data", False)) else 0.0
    enabled_windows = config.get("enabled_excluded_windows") or {}
    resolved = []
    for window_name, bounds in (config.get("excluded_windows") or {}).items():
        if not bool(enabled_windows.get(window_name, True)):
            continue
        resolved.extend(
            _resolve_particle_subtraction_window_collection(
                bounds,
                offset,
                "Particle subtraction component-fit excluded window '{}' for '{}'".format(
                    window_name,
                    fit_target,
                ),
            )
        )
    return resolved


def resolve_particle_subtraction_component_stage_amplitude_windows(
    fit_target,
    mm_offset_data=0.0,
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    config = get_particle_subtraction_component_fit_window_config(
        fit_target,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        setting_key=setting_key,
    )
    if not config:
        return {}

    offset = float(mm_offset_data) if bool(config.get("apply_mm_offset_data", False)) else 0.0
    resolved = {}
    for component_name, bounds in (config.get("stage_amplitude_windows") or {}).items():
        resolved[str(component_name)] = _resolve_particle_subtraction_window_collection(
            bounds,
            offset,
            "Particle subtraction component stage-amplitude window '{}' for '{}'".format(
                component_name,
                fit_target,
            ),
        )
    return resolved


def resolve_particle_subtraction_component_stage_amplitude_modes(
    fit_target,
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    config = get_particle_subtraction_component_fit_window_config(
        fit_target,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        setting_key=setting_key,
    )
    if not config:
        return {}

    resolved = {}
    for component_name, mode_name in (config.get("stage_amplitude_modes") or {}).items():
        normalized_mode = str(mode_name or "").strip().lower() or "least_squares"
        if normalized_mode not in ("least_squares", "window_integral"):
            raise ValueError(
                "Particle subtraction component stage-amplitude mode '{}' for '{}' must be 'least_squares' or 'window_integral'".format(
                    component_name,
                    fit_target,
                )
            )
        resolved[str(component_name)] = normalized_mode
    return resolved


def resolve_particle_subtraction_component_prior_scales(
    fit_target,
    component_names=None,
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    config = get_particle_subtraction_component_fit_window_config(
        fit_target,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        setting_key=setting_key,
    )
    if not config:
        return {}

    configured_scales = config.get("prior_scales") or {}
    requested_names = [str(name) for name in (component_names or []) if str(name)]
    if not requested_names:
        requested_names = list(configured_scales.keys()) or ["pi_n", "pi_delta", "pi_sidis"]

    resolved = {}
    for component_name in requested_names:
        if component_name in configured_scales:
            scale_value = configured_scales.get(component_name, 1.0)
        else:
            legacy_key = PARTICLE_SUBTRACTION_COMPONENT_PRIOR_SCALE_LEGACY_KEYS.get(component_name)
            scale_value = config.get(legacy_key, 1.0) if legacy_key else 1.0
        scale_value = float(scale_value)
        if scale_value != scale_value or scale_value < 0.0:
            raise ValueError(
                "Particle subtraction prior scale '{}' for '{}' must be finite and non-negative".format(
                    component_name,
                    fit_target,
                )
            )
        resolved[component_name] = scale_value
    return resolved


def resolve_particle_subtraction_component_postfit_scales(
    fit_target,
    component_names=None,
    inp_dict=None,
    phi_setting=None,
    setting_key=None,
):
    config = get_particle_subtraction_component_fit_window_config(
        fit_target,
        inp_dict=inp_dict,
        phi_setting=phi_setting,
        setting_key=setting_key,
    )
    if not config:
        return {}

    configured_scales = config.get("postfit_component_scales") or {}
    requested_names = [str(name) for name in (component_names or []) if str(name)]
    if not requested_names:
        requested_names = ["pi_n", "pi_delta", "pi_sidis"]
        for component_name in configured_scales.keys():
            component_name = str(component_name)
            if component_name not in requested_names:
                requested_names.append(component_name)

    resolved = {}
    for component_name in requested_names:
        scale_value = float(configured_scales.get(component_name, 1.0))
        if scale_value != scale_value or scale_value < 0.0:
            raise ValueError(
                "Particle subtraction post-fit scale '{}' for '{}' must be finite and non-negative".format(
                    component_name,
                    fit_target,
                )
            )
        resolved[component_name] = scale_value
    return resolved


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
