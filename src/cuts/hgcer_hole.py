#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-09-26 17:13:00 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import csv
from functools import lru_cache
from pathlib import Path

import ROOT
from ROOT import TCutG


HGCER_HOLE_DIR = Path(__file__).resolve().parents[2] / "input" / "kaon" / "hgcer_hole"


def _format_kinematic(value, decimals):
    value_text = str(value).strip()
    if value_text.startswith(("Q", "W")):
        value_text = value_text[1:]
    if "p" in value_text:
        return value_text
    return f"{float(value_text):.{decimals}f}".replace(".", "p")


def _format_epsset(epsset):
    eps_text = str(epsset).strip().lower()
    if eps_text.endswith("e"):
        return eps_text
    if eps_text in ("low", "high"):
        return f"{eps_text}e"
    raise ValueError(f"Invalid EPSSET for HGCer hole cut: {epsset}")


def _format_phi_setting(phi_setting):
    if phi_setting is None:
        raise ValueError("HGCer hole cut requires phi_setting: Right, Left, or Center.")

    phi_text = str(phi_setting).strip().lower()
    if phi_text not in ("right", "left", "center"):
        raise ValueError(f"Invalid HGCer hole phi_setting: {phi_setting}")
    return phi_text


def _hgcer_hole_path(Q2, W, EPSSET, phi_setting, hole_dir=None):
    q2_tag = _format_kinematic(Q2, 1)
    w_tag = _format_kinematic(W, 2)
    eps_tag = _format_epsset(EPSSET)
    phi_tag = _format_phi_setting(phi_setting)
    base_dir = Path(hole_dir) if hole_dir is not None else HGCER_HOLE_DIR

    return base_dir / f"Q{q2_tag}_W{w_tag}_{phi_tag}_{eps_tag}.csv"


@lru_cache(maxsize=None)
def _load_hgcer_points(csv_path):
    csv_path = Path(csv_path)
    if not csv_path.is_file():
        available = sorted(path.name for path in csv_path.parent.glob("*.csv")) if csv_path.parent.is_dir() else []
        available_msg = ", ".join(available) if available else "none found"
        raise FileNotFoundError(
            f"HGCer hole cut file not found: {csv_path}. Available CSV files: {available_msg}"
        )

    points = []
    with csv_path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        if reader.fieldnames is None or "x" not in reader.fieldnames or "y" not in reader.fieldnames:
            raise ValueError(f"HGCer hole cut file must have x,y columns: {csv_path}")

        for row in reader:
            points.append((float(row["x"]), float(row["y"])))

    if len(points) < 3:
        raise ValueError(f"HGCer hole cut file needs at least three points: {csv_path}")

    if points[0] != points[-1]:
        points.append(points[0])

    return tuple(points)


def apply_HGCer_hole_cut(Q2, W, EPSSET, phi_setting=None, hole_dir=None):
    csv_path = _hgcer_hole_path(Q2, W, EPSSET, phi_setting, hole_dir)
    points = _load_hgcer_points(csv_path)

    cutg_name = f"cutg_hgcer_hole_{csv_path.stem}"
    cutg = TCutG(cutg_name, len(points))
    cutg.SetVarX("P_hgcer_xAtCer")
    cutg.SetVarY("P_hgcer_yAtCer")

    for i, (x_pos, y_pos) in enumerate(points):
        cutg.SetPoint(i, x_pos, y_pos)

    return cutg
