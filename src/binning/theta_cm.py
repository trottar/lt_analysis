#! /usr/bin/python

import math

MP = 0.93827231
MP2 = 0.88035493
MPI = 0.13956995
MPI2 = 0.01947977
MN = 0.93956563
MN2 = 0.88278357
MLAMBDA = 1.115683
MLAMBDA2 = 1.244749
MK = 0.493677
MK2 = 0.24387


def _select_particle_mass_sq(particle_type):
    if particle_type == "kaon":
        return MK, MK2
    if particle_type == "pion":
        return MPI, MPI2
    raise ValueError("Unsupported particle type: {}".format(particle_type))


def _select_target_recoil_mass_sq(pol):
    if float(pol) > 0.0:
        return MP, MP2, MLAMBDA, MLAMBDA2
    return MN, MN2, MP, MP2


def _calculate_tmin_and_denom(particle_type, pol, w, q2):
    m3, m32 = _select_particle_mass_sq(particle_type)
    m2, m22, m4, m42 = _select_target_recoil_mass_sq(pol)

    if w <= 0.0 or q2 < 0.0:
        return float("nan"), float("nan")

    s = w * w
    omega = (s + q2 - m22) / (2.0 * m2)
    q = math.sqrt(max(q2 + omega * omega, 0.0))
    m12 = -q2

    e1cm = (s + m12 - m22) / (2.0 * w)
    e3cm = (s + m32 - m42) / (2.0 * w)
    p1cm = q * m2 / w
    p3cm_sq = e3cm * e3cm - m32

    if p1cm <= 0.0 or p3cm_sq <= 0.0:
        return float("nan"), float("nan")

    p3cm = math.sqrt(p3cm_sq)
    denom = 4.0 * p1cm * p3cm
    if denom <= 0.0:
        return float("nan"), float("nan")

    tmin = -((e1cm - e3cm) ** 2 - (p1cm - p3cm) ** 2)
    if not math.isfinite(tmin):
        return float("nan"), float("nan")

    return tmin, denom


def calculate_tmin(particle_type, pol, w, q2):
    tmin, _ = _calculate_tmin_and_denom(particle_type, pol, w, q2)
    return tmin


def _calculate_half_angle_sin_sq(particle_type, pol, w, q2, minus_t):
    if minus_t < 0.0:
        return float("nan")

    tmin, denom = _calculate_tmin_and_denom(particle_type, pol, w, q2)
    if not math.isfinite(tmin) or not math.isfinite(denom):
        return float("nan")

    sin_half_sq = (minus_t - tmin) / denom

    if not math.isfinite(sin_half_sq):
        return float("nan")

    return min(max(sin_half_sq, 0.0), 1.0)


def calculate_sin_theta_cm(particle_type, pol, w, q2, minus_t):
    sin_half_sq = _calculate_half_angle_sin_sq(particle_type, pol, w, q2, minus_t)
    if not math.isfinite(sin_half_sq):
        return float("nan")

    return 2.0 * math.sqrt(sin_half_sq * max(1.0 - sin_half_sq, 0.0))


def calculate_theta_cm_rad(particle_type, pol, w, q2, minus_t):
    sin_half_sq = _calculate_half_angle_sin_sq(particle_type, pol, w, q2, minus_t)
    if not math.isfinite(sin_half_sq):
        return float("nan")

    return 2.0 * math.asin(math.sqrt(sin_half_sq))


def calculate_theta_cm_deg(particle_type, pol, w, q2, minus_t):
    theta_cm = calculate_theta_cm_rad(particle_type, pol, w, q2, minus_t)
    if not math.isfinite(theta_cm):
        return float("nan")
    return math.degrees(theta_cm)
