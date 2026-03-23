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

_TOL = 1.0e-10


def _select_particle_mass_sq(particle_type):
    if particle_type == "kaon":
        return MK, MK2
    if particle_type == "pion":
        return MPI, MPI2
    raise ValueError(f"Unsupported particle type: {particle_type}")


def _select_target_recoil_mass_sq(pol):
    if float(pol) > 0.0:
        return MP, MP2, MLAMBDA, MLAMBDA2
    return MN, MN2, MP, MP2


def _calculate_half_angle_sin_sq(particle_type, pol, w, q2, minus_t):
    m3, m32 = _select_particle_mass_sq(particle_type)
    m2, m22, m4, m42 = _select_target_recoil_mass_sq(pol)

    if w <= 0.0 or q2 < 0.0 or minus_t < 0.0:
        return float("nan")

    s = w * w
    omega = (s + q2 - m22) / (2.0 * m2)
    q = math.sqrt(max(q2 + omega * omega, 0.0))
    m12 = -q2  # virtual photon mass^2

    e1cm = (s + m12 - m22) / (2.0 * w)
    e3cm = (s + m32 - m42) / (2.0 * w)
    p1cm = q * m2 / w

    p3cm_sq = e3cm * e3cm - m32
    if p1cm <= 0.0 or p3cm_sq < -_TOL:
        return float("nan")
    p3cm = math.sqrt(max(p3cm_sq, 0.0))

    denom = 4.0 * p1cm * p3cm
    if denom <= _TOL:
        return float("nan")

    minus_t_min = -((e1cm - e3cm) ** 2 - (p1cm - p3cm) ** 2)
    raw = (minus_t - minus_t_min) / denom

    # reject truly unphysical points; only clip tiny floating-point leakage
    if raw < -_TOL or raw > 1.0 + _TOL:
        return float("nan")

    return min(max(raw, 0.0), 1.0)


def calculate_sin_theta_cm(particle_type, pol, w, q2, minus_t):
    sin_half_sq = _calculate_half_angle_sin_sq(particle_type, pol, w, q2, minus_t)
    if not math.isfinite(sin_half_sq):
        return float("nan")

    return 2.0 * math.sqrt(sin_half_sq * (1.0 - sin_half_sq))


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