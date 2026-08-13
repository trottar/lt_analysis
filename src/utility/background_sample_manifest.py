#! /usr/bin/python

"""Runtime source manifest for SIMC background templates."""

from __future__ import annotations

import os


BACKGROUND_SAMPLE_PHI_ENV = {
    "Right": "RIGHT",
    "Left": "LEFT",
    "Center": "CENTER",
}

BACKGROUND_SAMPLE_SUFFIXES = {
    "base": "",
    "root": ".root",
    "hist": ".hist",
    "gen": ".gen",
    "geni": ".geni",
    "random_state": "_start_random_state.dat",
}

GENERATED_BACKGROUND_SAMPLES = ("neutron", "delta", "sidis")
EXTERNAL_BACKGROUND_SAMPLES = ("sigma0",)


def _source_identity(q2, w, epsset, phi_setting):
    return {
        "Q2": str(q2),
        "W": str(w),
        "EPSSET": str(epsset),
        "phi_setting": str(phi_setting),
    }


def _explicit_sigma0_root(phi_env):
    value = str(os.environ.get("LT_BG_SIGMA0_{}_ROOT".format(phi_env), "")).strip()
    return value or None


def build_background_sample_manifest(ltanapath, q2, w, epsset):
    """Build generated and explicitly configured external SIMC sources.

    K-Sigma0 is deliberately excluded from the generated-background convention.
    It remains present in every per-phi manifest so consumers can distinguish an
    intentional unconfigured external source from a missing manifest entry.
    """
    env_base_dir = os.environ.get("LT_BG_SAMPLE_BASE")
    base_dir = env_base_dir or os.path.join(ltanapath, "background_samples", "OUTPUTS")
    sample_q2 = os.environ.get("LT_BG_SAMPLE_Q2", q2)
    sample_w = os.environ.get("LT_BG_SAMPLE_W", w)
    sample_eps = os.environ.get("LT_BG_SAMPLE_EPSILON", epsset)
    manifest = {
        "base_dir": base_dir,
        "q2": sample_q2,
        "w": sample_w,
        "epsilon": sample_eps,
        "by_background": {},
        "by_phi": {phi_label: {} for phi_label in BACKGROUND_SAMPLE_PHI_ENV},
    }

    for background in GENERATED_BACKGROUND_SAMPLES:
        env_background = background.upper()
        manifest["by_background"][background] = {}
        for phi_label, phi_env in BACKGROUND_SAMPLE_PHI_ENV.items():
            env_prefix = "LT_BG_{}_{}".format(env_background, phi_env)
            sample_name = "Prod_Coin_Q{}W{}{}_{}e".format(
                sample_q2,
                sample_w,
                phi_label.lower(),
                sample_eps,
            )
            default_base = os.path.join(base_dir, background, sample_name)
            sample_entry = {
                "source_strategy": "generated_background",
                "configured": True,
                "source_identity": _source_identity(
                    sample_q2,
                    sample_w,
                    sample_eps,
                    phi_label,
                ),
            }
            for key, suffix in BACKGROUND_SAMPLE_SUFFIXES.items():
                env_name = "{}_{}".format(env_prefix, key.upper())
                sample_entry[key] = os.environ.get(env_name, "{}{}".format(default_base, suffix))
            manifest["by_background"][background][phi_label] = sample_entry
            manifest["by_phi"][phi_label][background] = sample_entry

    manifest["by_background"]["sigma0"] = {}
    for phi_label, phi_env in BACKGROUND_SAMPLE_PHI_ENV.items():
        root_filename = _explicit_sigma0_root(phi_env)
        sample_entry = {
            "source_strategy": "external_required",
            "configured": bool(root_filename),
            "source_identity": _source_identity(
                sample_q2,
                sample_w,
                sample_eps,
                phi_label,
            ),
        }
        for key in BACKGROUND_SAMPLE_SUFFIXES:
            sample_entry[key] = root_filename if key == "root" else None
        manifest["by_background"]["sigma0"][phi_label] = sample_entry
        manifest["by_phi"][phi_label]["sigma0"] = sample_entry

    return manifest
