#! /usr/bin/python
"""Create one low/high-epsilon canonical t/phi interval pair from raw support.

The production launcher first captures the selected pre-particle records for
both epsilon settings, then invokes this small orchestration helper before
either full analysis starts.  It intentionally knows nothing about pion fits
or normalized yields.
"""

import argparse
import os
import sys

import numpy as np


SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
SRC_DIR = os.path.dirname(SCRIPT_DIR)
for path in (SCRIPT_DIR, os.path.join(SRC_DIR, "binning")):
    if path not in sys.path:
        sys.path.insert(0, path)

from background_config import resolve_analysis_runtime_config  # noqa: E402
from find_bins import resolve_shared_canonical_phi_preflight  # noqa: E402


def _load_capture(path, epsilon):
    with np.load(path) as captured:
        t_values = np.asarray(captured["adj_t"], dtype=float)
        phi_values = np.asarray(captured["phi_value"], dtype=float)
        coefficients = np.asarray(captured["physical_coefficient"], dtype=float)
    if not (len(t_values) == len(phi_values) == len(coefficients)):
        raise ValueError("canonical preflight capture has inconsistent array lengths: {}".format(path))
    return [
        {
            "adj_t": float(t_value),
            "phi_value": float(phi_value),
            "physical_coefficient": float(coefficient),
            "source_label": "{}_captured_prepass".format(epsilon),
            "entry_index": int(index),
        }
        for index, (t_value, phi_value, coefficient) in enumerate(
            zip(t_values, phi_values, coefficients)
        )
    ]


def build_preflight_input(q2, w, particle_type):
    """Build the minimal canonical resolver contract from central config."""
    config = resolve_analysis_runtime_config(q2, w)
    return {
        "ParticleType": str(particle_type),
        "Q2": str(q2),
        "W": str(w),
        "EPSSET": "low",
        "NumtBins": int(config["num_t_bins"]),
        "NumPhiBins": int(config["num_phi_bins"]),
        "tmin": float(config["tmin"]),
        "tmax": float(config["tmax"]),
        "auto_rebin_phi": bool(config["auto_rebin_phi"]),
        "min_phi_bins": int(config["min_phi_bins"]),
        "t_phi_support_policy": str(config["t_phi_support_policy"]),
        "t_phi_support_min_events": int(config["t_phi_support_min_events"]),
        "analysis_runtime_config_hash": config["config_hash"],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("q2")
    parser.add_argument("w")
    parser.add_argument("particle_type")
    parser.add_argument("low_capture")
    parser.add_argument("high_capture")
    arguments = parser.parse_args(argv)
    inp_dict = build_preflight_input(arguments.q2, arguments.w, arguments.particle_type)
    result = resolve_shared_canonical_phi_preflight(
        {
            "low": [{"records": _load_capture(arguments.low_capture, "low")}],
            "high": [{"records": _load_capture(arguments.high_capture, "high")}],
        },
        inp_dict,
        write_interval_files=True,
    )
    metadata = result["metadata"]
    print(
        "Shared canonical preflight complete: Nt={} Nphi={} pair={}".format(
            len(result["t_bins"]) - 1,
            len(result["phi_bins"]) - 1,
            metadata["canonical_interval_pair_id"],
        )
    )


if __name__ == "__main__":
    main()
