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
SIGMA0_EPSILON_TOKENS = ("low", "high")


def _source_identity(q2, w, epsset, phi_setting):
    return {
        "Q2": str(q2),
        "W": str(w),
        "EPSSET": str(epsset),
        "phi_setting": str(phi_setting),
    }


def normalize_sigma0_epsilon(epsset):
    """Return the only supported K-Sigma0 epsilon token.

    The external K-Sigma0 source is keyed by epsilon.  Accepting anything
    other than ``low`` or ``high`` would make it possible to select a source
    from a different analysis setting, so this normalizer is deliberately
    stricter than the generated-background naming convention.
    """
    normalized = str(epsset or "").strip().lower()
    if normalized not in SIGMA0_EPSILON_TOKENS:
        raise ValueError(
            "K-Sigma0 EPSSET must be one of {}; received {!r}".format(
                ", ".join(SIGMA0_EPSILON_TOKENS), epsset
            )
        )
    return normalized


def sigma0_environment_variable(phi_setting, epsset):
    """Build the authoritative K-Sigma0 external-source variable name."""
    try:
        phi_env = BACKGROUND_SAMPLE_PHI_ENV[phi_setting]
    except KeyError as exc:
        raise ValueError(
            "Unknown K-Sigma0 phi setting {!r}; expected one of {}".format(
                phi_setting, ", ".join(BACKGROUND_SAMPLE_PHI_ENV)
            )
        ) from exc
    epsilon = normalize_sigma0_epsilon(epsset)
    return "LT_BG_SIGMA0_{}_{}_ROOT".format(phi_env, epsilon.upper())


def _explicit_sigma0_root(phi_setting, epsset, environ=None):
    environment_variable = sigma0_environment_variable(phi_setting, epsset)
    source_environ = os.environ if environ is None else environ
    value = str(source_environ.get(environment_variable, "")).strip()
    return value or None


def build_background_sample_manifest(ltanapath, q2, w, epsset):
    """Build generated and explicitly configured external SIMC sources.

    K-Sigma0 is deliberately excluded from the generated-background convention.
    It remains present in every per-phi manifest so consumers can distinguish an
    intentional unconfigured external source from a missing manifest entry.
    """
    analysis_q2 = str(q2)
    analysis_w = str(w)
    analysis_eps = normalize_sigma0_epsilon(epsset)

    env_base_dir = os.environ.get("LT_BG_SAMPLE_BASE")
    base_dir = env_base_dir or os.path.join(ltanapath, "background_samples", "OUTPUTS")
    generated_q2 = str(os.environ.get("LT_BG_SAMPLE_Q2", q2))
    generated_w = str(os.environ.get("LT_BG_SAMPLE_W", w))
    generated_eps = normalize_sigma0_epsilon(
        os.environ.get("LT_BG_SAMPLE_EPSILON", epsset)
    )
    manifest = {
        "base_dir": base_dir,
        "q2": generated_q2,
        "w": generated_w,
        "epsilon": generated_eps,
        "by_background": {},
        "by_phi": {phi_label: {} for phi_label in BACKGROUND_SAMPLE_PHI_ENV},
    }

    for background in GENERATED_BACKGROUND_SAMPLES:
        env_background = background.upper()
        manifest["by_background"][background] = {}
        for phi_label, phi_env in BACKGROUND_SAMPLE_PHI_ENV.items():
            env_prefix = "LT_BG_{}_{}".format(env_background, phi_env)
            sample_name = "Prod_Coin_Q{}W{}{}_{}e".format(
                generated_q2,
                generated_w,
                phi_label.lower(),
                generated_eps,
            )
            default_base = os.path.join(base_dir, background, sample_name)
            sample_entry = {
                "source_strategy": "generated_background",
                "configured": True,
                "source_identity": _source_identity(
                    generated_q2,
                    generated_w,
                    generated_eps,
                    phi_label,
                ),
            }
            for key, suffix in BACKGROUND_SAMPLE_SUFFIXES.items():
                env_name = "{}_{}".format(env_prefix, key.upper())
                sample_entry[key] = os.environ.get(env_name, "{}{}".format(default_base, suffix))
            manifest["by_background"][background][phi_label] = sample_entry
            manifest["by_phi"][phi_label][background] = sample_entry

    manifest["by_background"]["sigma0"] = {}
    for phi_label in BACKGROUND_SAMPLE_PHI_ENV:
        environment_variable = sigma0_environment_variable(phi_label, analysis_eps)
        root_filename = _explicit_sigma0_root(phi_label, analysis_eps)
        sample_entry = {
            "source_strategy": "external_required",
            "configured": bool(root_filename),
            "environment_variable": environment_variable,
            "source_identity": _source_identity(
                analysis_q2,
                analysis_w,
                analysis_eps,
                phi_label,
            ),
        }
        for key in BACKGROUND_SAMPLE_SUFFIXES:
            sample_entry[key] = root_filename if key == "root" else None
        manifest["by_background"]["sigma0"][phi_label] = sample_entry
        manifest["by_phi"][phi_label]["sigma0"] = sample_entry

    return manifest
