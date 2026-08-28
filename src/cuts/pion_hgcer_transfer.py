"""Frozen, analysis-stage HGCer zero-photoelectron transfer diagnostics.

The module intentionally has no production subtraction entry point.  It reads
Part-1 records, audits their noRF/PID contract, freezes a response map, and may
construct detached proposed MM products for review by a later Part-3 decision.
"""

from __future__ import annotations

from array import array
from copy import deepcopy
import hashlib
import json
import math
import os

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import gammainc, gammaln

from background_config import (
    fingerprint_hgcer_pid_contract,
    hgcer_mask_accepts,
    normalize_hgcer_mask,
)

try:  # Keep pure-Python map tests usable without a PyROOT installation.
    import ROOT
except ImportError:  # pragma: no cover
    ROOT = None

try:
    from root_histogram_ownership import clone_root_histogram
except ImportError:  # pragma: no cover - only needed by PyROOT integration
    clone_root_histogram = None


ROOT_SAFE_TRANSFER_LABELS = {
    "part2": "PION HGCer ZERO-PE TRANSFER MAP - PART 2 - NON-AUTHORITATIVE",
    "unavailable": "PION HGCer PART 2 UNAVAILABLE - PRODUCTION UNCHANGED",
    "mask_audit": "Part 2 mask and noRF provenance audit",
    "family_audit": "Part 2 response-family audit",
    "p0": "Pion zero-photoelectron probability P0",
    "relative_p0": "Relative P0 statistical uncertainty",
    "transfer": "Proposed pion-to-kaon transfer R",
    "uncertainty": "Proposed transfer total uncertainty",
    "statistical_uncertainty": "Proposed transfer statistical uncertainty",
    "model_uncertainty": "Proposed transfer model uncertainty",
    "fallback_uncertainty": "Proposed transfer fallback uncertainty",
    "fit_status": "Part 2 fit status",
    "solution_source": "Part 2 solution source",
    "eligibility": "Part 3 review eligibility by cell",
    "response_fit": "Prompt-pion HGCer response fit",
    "simc_closure": "SIMC closure view - diagnostic only",
    "signal_closure": "K-Lambda and K-Sigma0 closure view - diagnostic only",
    "proposed_mm": "Proposed pion subtraction from proton-cleaned noRF host",
}


def pion_hgcer_transfer_display_text(kind):
    try:
        return ROOT_SAFE_TRANSFER_LABELS[str(kind)]
    except KeyError as exc:
        raise ValueError("unknown Part-2 display text '{}'".format(kind)) from exc


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(child) for key, child in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_ready(child) for child in value]
    if isinstance(value, np.generic):
        return _json_ready(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def _fingerprint(value):
    encoded = json.dumps(_json_ready(value), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _finite(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _mask_probability_discrete(mu, mask):
    """Poisson probability for the finite comparison language in the contract."""
    normalized = normalize_hgcer_mask(mask)
    operator, threshold = normalized["operator"], normalized["value"]
    if operator == "==":
        integer = int(round(threshold))
        return math.exp(-mu + integer * math.log(mu) - math.lgamma(integer + 1)) if (
            integer >= 0 and abs(threshold - integer) < 1.0e-12
        ) else 0.0
    floor = int(math.floor(threshold))
    cdf = sum(
        math.exp(-mu + number * math.log(mu) - math.lgamma(number + 1))
        for number in range(max(0, floor) + 1)
    )
    if operator == ">":
        return max(0.0, 1.0 - cdf)
    if operator == ">=":
        lower = int(math.ceil(threshold))
        return max(0.0, 1.0 - sum(
            math.exp(-mu + number * math.log(mu) - math.lgamma(number + 1))
            for number in range(max(0, lower))
        ))
    if operator == "<":
        return sum(
            math.exp(-mu + number * math.log(mu) - math.lgamma(number + 1))
            for number in range(max(0, int(math.ceil(threshold))))
        )
    return sum(
        math.exp(-mu + number * math.log(mu) - math.lgamma(number + 1))
        for number in range(max(0, int(math.floor(threshold)) + 1))
    )


def poisson_zero_probability(mu):
    return math.exp(-float(mu))


def zero_truncated_poisson_efficiency(mu, mask):
    """Unconditional detector efficiency; P0 remains structurally e^-mu."""
    return _mask_probability_discrete(float(mu), mask)


def zero_truncated_poisson_transfer(mu, kaon_mask, physical_control_mask):
    numerator = zero_truncated_poisson_efficiency(mu, kaon_mask)
    denominator = zero_truncated_poisson_efficiency(mu, physical_control_mask)
    return numerator / denominator if denominator > 0.0 else None


def validate_zero_photoelectron_response_family(name):
    """Reject positive-only families before they can define a zero-NPE leak."""
    family = str(name or "").strip().lower()
    if family in {"gamma", "gaussian", "positive_gamma", "positive_gaussian", "zero_inflated_gamma"}:
        raise ValueError("positive-only response family cannot determine P_hgcer_npeSum == 0 leakage")
    if family not in {
        "zero_truncated_poisson", "zero_truncated_negative_binomial",
        "zero_truncated_compound_poisson_gamma_gain",
        "zero_truncated_compound_poisson_exponential_gain", "auto",
    }:
        raise ValueError("unsupported zero-photoelectron response family '{}'".format(name))
    return family


def _poisson_logpmf(values, mu):
    values = np.asarray(values, dtype=float)
    return -mu + values * math.log(mu) - gammaln(values + 1.0)


def _poisson_fit_normalization(mu, response_mask, fit_range):
    lower, upper = [float(value) for value in fit_range]
    start = max(0, int(math.ceil(lower)))
    stop = int(math.floor(upper))
    return sum(
        math.exp(-mu + number * math.log(mu) - math.lgamma(number + 1))
        for number in range(start, stop + 1)
        if hgcer_mask_accepts(number, response_mask)
    )


def _finite_hessian_covariance(fun, parameters, step=1.0e-4):
    parameters = np.asarray(parameters, dtype=float)
    size = len(parameters)
    hessian = np.zeros((size, size), dtype=float)
    center = float(fun(parameters))
    for left in range(size):
        unit_left = np.zeros(size)
        unit_left[left] = step
        hessian[left, left] = (float(fun(parameters + unit_left)) - 2.0 * center + float(fun(parameters - unit_left))) / (step * step)
        for right in range(left + 1, size):
            unit_right = np.zeros(size)
            unit_right[right] = step
            value = (
                float(fun(parameters + unit_left + unit_right))
                - float(fun(parameters + unit_left - unit_right))
                - float(fun(parameters - unit_left + unit_right))
                + float(fun(parameters - unit_left - unit_right))
            ) / (4.0 * step * step)
            hessian[left, right] = hessian[right, left] = value
    try:
        covariance = np.linalg.inv(hessian)
        condition = float(np.linalg.cond(hessian))
    except np.linalg.LinAlgError:
        covariance, condition = None, float("inf")
    return covariance, condition


def _profile_log_mu(nll, best_parameters, config, profile_nll=None):
    best_parameters = np.asarray(best_parameters, dtype=float)
    best_nll = float(nll(best_parameters))
    center = float(best_parameters[0])
    grid = np.linspace(
        center - float(config["profile_log_mu_half_range"]),
        center + float(config["profile_log_mu_half_range"]),
        int(config["profile_points"]),
    )
    values = []
    for log_mu in grid:
        if profile_nll is None or len(best_parameters) == 1:
            value = float(nll(np.asarray([log_mu])))
        else:
            initial = best_parameters[1:]
            try:
                result = minimize(
                    lambda rest: profile_nll(float(log_mu), np.asarray(rest, dtype=float)),
                    initial, method="L-BFGS-B", bounds=[(-10.0, 8.0)] * len(initial),
                )
                value = float(result.fun) if result.success and math.isfinite(result.fun) else float("inf")
            except (ArithmeticError, FloatingPointError, ValueError):
                value = float("inf")
        values.append(value)
    target = best_nll + float(config["profile_delta_nll"])
    lower, upper = None, None
    for index in range(len(grid) - 1):
        left, right = values[index] - target, values[index + 1] - target
        if left == 0.0 or left * right <= 0.0:
            crossing = float(grid[index] + (grid[index + 1] - grid[index]) * (-left) / (right - left)) if right != left else float(grid[index])
            if crossing < center and lower is None:
                lower = crossing
            if crossing > center and upper is None:
                upper = crossing
    status = "two_sided" if lower is not None and upper is not None else "one_sided_or_unbounded"
    return {
        "status": status,
        "delta_nll": float(config["profile_delta_nll"]),
        "log_mu_interval": [lower, upper],
        "mu_interval": [math.exp(lower) if lower is not None else None, math.exp(upper) if upper is not None else None],
        "grid_log_mu": list(map(float, grid)),
        "grid_nll": [value if math.isfinite(value) else None for value in values],
    }


def _fit_poisson(values, response_mask, contract, config):
    values = np.asarray(values, dtype=float)
    fit_range = config["fit_range"]
    fit_values = values[(values >= fit_range[0]) & (values <= fit_range[1])]
    if len(fit_values) == 0 or np.any(np.abs(fit_values - np.rint(fit_values)) > float(config["integer_tolerance"])):
        return {"fit_status": "fit_invalid", "reason": "integer_poisson_fit_requires_integer_values"}
    rounded = np.rint(fit_values).astype(int)

    def nll(parameters):
        mu = math.exp(float(parameters[0]))
        normalization = _poisson_fit_normalization(mu, response_mask, fit_range)
        if normalization <= 0.0:
            return float("inf")
        return float(-np.sum(_poisson_logpmf(rounded, mu)) + len(rounded) * math.log(normalization))

    result = minimize_scalar(lambda log_mu: nll(np.asarray([log_mu])), bounds=(-12.0, 8.0), method="bounded")
    if not result.success or not math.isfinite(float(result.fun)):
        return {"fit_status": "fit_invalid", "reason": "poisson_optimizer_failed"}
    parameters = np.asarray([float(result.x)])
    covariance, condition = _finite_hessian_covariance(nll, parameters)
    return _complete_fit(
        family="zero_truncated_poisson", values=rounded, nll=nll, parameters=parameters,
        covariance=covariance, condition=condition, response_mask=response_mask,
        contract=contract, config=config,
    )


def _nbinom_logpmf(values, mean, shape):
    probability = shape / (shape + mean)
    values = np.asarray(values, dtype=float)
    return gammaln(values + shape) - gammaln(shape) - gammaln(values + 1.0) + shape * math.log(probability) + values * math.log(1.0 - probability)


def _fit_negative_binomial(values, response_mask, contract, config):
    values = np.asarray(values, dtype=float)
    fit_range = config["fit_range"]
    selected = values[(values >= fit_range[0]) & (values <= fit_range[1])]
    if len(selected) == 0 or np.any(np.abs(selected - np.rint(selected)) > float(config["integer_tolerance"])):
        return {"fit_status": "fit_invalid", "reason": "integer_negative_binomial_fit_requires_integer_values"}
    selected = np.rint(selected).astype(int)

    def nll(parameters):
        mean, shape = math.exp(float(parameters[0])), math.exp(float(parameters[1]))
        probability = shape / (shape + mean)
        normalization = sum(
            math.exp(float(_nbinom_logpmf([number], mean, shape)[0]))
            for number in range(max(0, int(math.ceil(fit_range[0]))), int(math.floor(fit_range[1])) + 1)
            if hgcer_mask_accepts(number, response_mask)
        )
        if normalization <= 0.0:
            return float("inf")
        return float(-np.sum(_nbinom_logpmf(selected, mean, shape)) + len(selected) * math.log(normalization))

    start = np.log([max(float(np.mean(selected)), 1.0e-3), 5.0])
    result = minimize(nll, start, method="L-BFGS-B", bounds=[(-12.0, 8.0), (-8.0, 12.0)])
    if not result.success or not math.isfinite(float(result.fun)):
        return {"fit_status": "fit_invalid", "reason": "negative_binomial_optimizer_failed"}
    covariance, condition = _finite_hessian_covariance(nll, result.x)
    # The positive NB component is an alternate shape diagnostic.  Its zero
    # atom remains the contract-linked e^-mu rather than a free inflation term.
    return _complete_fit(
        family="zero_truncated_negative_binomial", values=selected, nll=nll,
        parameters=np.asarray(result.x), covariance=covariance, condition=condition,
        response_mask=response_mask, contract=contract, config=config,
        positive_nbinom_shape=math.exp(float(result.x[1])),
    )


def _compound_density(values, mu, shape, scale):
    values = np.asarray(values, dtype=float)
    density = np.zeros_like(values, dtype=float)
    for count in range(1, 80):
        poisson = math.exp(-mu + count * math.log(mu) - math.lgamma(count + 1))
        density += poisson * np.exp((count * shape - 1.0) * np.log(np.maximum(values, 1.0e-300)) - values / scale - gammaln(count * shape) - count * shape * math.log(scale))
    return density


def _compound_interval_probability(mu, shape, scale, lower, upper):
    if not all(math.isfinite(value) and value > 0.0 for value in (mu, shape, scale)):
        return 0.0
    result = 0.0
    for count in range(1, 80):
        poisson = math.exp(-mu + count * math.log(mu) - math.lgamma(count + 1))
        result += poisson * (gammainc(count * shape, upper / scale) - gammainc(count * shape, lower / scale))
    return max(0.0, float(result))


def _fit_compound_poisson(values, response_mask, contract, config, *, exponential_gain=False):
    values = np.asarray(values, dtype=float)
    lower, upper = config["fit_range"]
    selected = values[(values >= lower) & (values <= upper)]
    if len(selected) == 0 or np.any(selected <= 0.0):
        return {"fit_status": "fit_invalid", "reason": "compound_poisson_requires_positive_continuous_values"}
    parameter_count = 2 if exponential_gain else 3

    def unpack(parameters):
        mu = math.exp(float(parameters[0]))
        shape = 1.0 if exponential_gain else math.exp(float(parameters[1]))
        scale = math.exp(float(parameters[-1]))
        return mu, shape, scale

    def nll(parameters):
        mu, shape, scale = unpack(parameters)
        normalization = _compound_interval_probability(mu, shape, scale, lower, upper)
        density = _compound_density(selected, mu, shape, scale)
        if normalization <= 0.0 or np.any(density <= 0.0) or not np.all(np.isfinite(density)):
            return float("inf")
        return float(-np.sum(np.log(density)) + len(selected) * math.log(normalization))

    observed_mean = max(float(np.mean(selected)), 0.1)
    starts = []
    for start_mu in (0.5, 1.0, 2.0, 4.0, 8.0):
        positive_fraction = 1.0 - math.exp(-start_mu)
        if exponential_gain:
            starts.append(np.log([start_mu, max(observed_mean * positive_fraction / start_mu, 0.05)]))
        else:
            for start_shape in (0.7, 1.5, 3.0, 6.0):
                starts.append(np.log([
                    start_mu, start_shape,
                    max(observed_mean * positive_fraction / (start_mu * start_shape), 0.05),
                ]))
    candidates = [
        minimize(nll, initial, method="L-BFGS-B", bounds=[(-10.0, 8.0)] * parameter_count)
        for initial in starts
    ]
    valid_candidates = [candidate for candidate in candidates if candidate.success and math.isfinite(float(candidate.fun))]
    result = min(valid_candidates, key=lambda candidate: float(candidate.fun)) if valid_candidates else None
    if result is None:
        return {"fit_status": "fit_invalid", "reason": "compound_poisson_optimizer_failed"}
    covariance, condition = _finite_hessian_covariance(nll, result.x)
    return _complete_fit(
        family="zero_truncated_compound_poisson_exponential_gain" if exponential_gain else "zero_truncated_compound_poisson_gamma_gain",
        values=selected, nll=nll, parameters=np.asarray(result.x), covariance=covariance,
        condition=condition, response_mask=response_mask, contract=contract, config=config,
        compound_parameters=unpack(result.x),
    )


def _complete_fit(*, family, values, nll, parameters, covariance, condition, response_mask, contract, config, positive_nbinom_shape=None, compound_parameters=None):
    mu = math.exp(float(parameters[0]))
    p0 = poisson_zero_probability(mu)
    kaon_mask = contract["masks"]["kaon_tree"]
    control_mask = contract["masks"]["physical_pion_control"]
    if compound_parameters is None and positive_nbinom_shape is None:
        kaon_efficiency = zero_truncated_poisson_efficiency(mu, kaon_mask)
        control_efficiency = zero_truncated_poisson_efficiency(mu, control_mask)
    elif compound_parameters is not None:
        _, shape, scale = compound_parameters
        kaon_efficiency = p0 if kaon_mask["operator"] == "==" and kaon_mask["value"] == 0.0 else None
        lower = 0.0
        if control_mask["operator"] == ">":
            lower = float(control_mask["value"])
        elif control_mask["operator"] == ">=":
            lower = float(control_mask["value"])
        else:
            return {"fit_status": "fit_invalid", "reason": "compound_poisson_control_mask_unsupported"}
        control_efficiency = _compound_interval_probability(mu, shape, scale, lower, 1.0e4)
    else:
        # Positive NB component with a zero atom fixed by mu.  Numerically sum
        # its conditional positive tail, then scale by (1-P0).
        mean, shape = mu, positive_nbinom_shape
        if kaon_mask["operator"] != "==" or kaon_mask["value"] != 0.0 or control_mask["operator"] != ">":
            return {"fit_status": "fit_invalid", "reason": "negative_binomial_mask_unsupported"}
        nbinom_zero = math.exp(float(_nbinom_logpmf([0], mean, shape)[0]))
        positive_tail = 1.0 - sum(math.exp(float(_nbinom_logpmf([number], mean, shape)[0])) for number in range(0, int(math.floor(control_mask["value"])) + 1))
        control_efficiency = (1.0 - p0) * positive_tail / max(1.0e-15, 1.0 - nbinom_zero)
        kaon_efficiency = p0
    if kaon_efficiency is None or control_efficiency is None or control_efficiency <= 0.0:
        return {"fit_status": "fit_invalid", "reason": "nonpositive_physical_control_efficiency"}
    transfer = float(kaon_efficiency / control_efficiency)
    covariance_status, max_correlation = "available", 0.0
    sigma_p0, sigma_transfer = None, None
    if (
        covariance is None or not np.all(np.isfinite(covariance))
        or np.any(np.diag(covariance) <= 0.0)
        or np.any(np.linalg.eigvalsh((covariance + covariance.T) / 2.0) <= 0.0)
    ):
        covariance_status = "unavailable"
    else:
        denominator = np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
        correlations = np.divide(covariance, denominator, out=np.zeros_like(covariance), where=denominator > 0.0)
        max_correlation = float(np.max(np.abs(correlations - np.eye(len(correlations))))) if len(correlations) > 1 else 0.0
        sigma_log_mu = math.sqrt(float(covariance[0, 0]))
        sigma_p0 = p0 * sigma_log_mu
    profile = _profile_log_mu(
        nll, parameters, config,
        profile_nll=(lambda log_mu, rest: nll(np.r_[log_mu, rest])) if len(parameters) > 1 else None,
    )
    rng_seed = int(_fingerprint({
        "contract": contract["fingerprint"], "family": family,
        "parameters": list(map(float, parameters)), "count": int(len(values)),
        "coordinate": config.get("_toy_context") or {},
    })[:16], 16) % (2 ** 32)
    rng = np.random.default_rng(rng_seed)
    toy_transfers, toy_p0 = [], []
    if covariance is not None and covariance_status == "available":
        for toy_parameters in rng.multivariate_normal(parameters, covariance, size=int(config["toy_count"])):
            if not np.all(np.isfinite(toy_parameters)):
                continue
            toy_mu = math.exp(float(toy_parameters[0]))
            if not math.isfinite(toy_mu) or toy_mu <= 0.0:
                continue
            numerator = poisson_zero_probability(toy_mu)
            if compound_parameters is not None:
                toy_shape = 1.0 if "exponential_gain" in family else math.exp(float(toy_parameters[1]))
                toy_scale = math.exp(float(toy_parameters[-1]))
                denominator = _compound_interval_probability(
                    toy_mu, toy_shape, toy_scale, float(control_mask["value"]), 1.0e4
                )
            elif positive_nbinom_shape is not None:
                toy_shape = math.exp(float(toy_parameters[1]))
                nbinom_zero = math.exp(float(_nbinom_logpmf([0], toy_mu, toy_shape)[0]))
                tail = 1.0 - sum(
                    math.exp(float(_nbinom_logpmf([number], toy_mu, toy_shape)[0]))
                    for number in range(0, int(math.floor(control_mask["value"])) + 1)
                )
                denominator = (1.0 - numerator) * tail / max(1.0e-15, 1.0 - nbinom_zero)
            else:
                denominator = zero_truncated_poisson_efficiency(toy_mu, control_mask)
            if denominator > 0.0 and math.isfinite(numerator / denominator):
                toy_p0.append(numerator)
                toy_transfers.append(numerator / denominator)
    toy_fraction = float(len(toy_transfers)) / float(max(1, int(config["toy_count"])))
    if toy_transfers:
        sigma_transfer = float(np.std(toy_transfers, ddof=1)) if len(toy_transfers) > 1 else 0.0
        sigma_p0 = float(np.std(toy_p0, ddof=1)) if len(toy_p0) > 1 else sigma_p0
    covariance_degenerate = covariance_status != "available" or condition > float(config["covariance_condition_number_max"]) or max_correlation >= float(config["pair_correlation_max"])
    physical = all(math.isfinite(value) and value >= 0.0 for value in (p0, kaon_efficiency, control_efficiency, transfer))
    if not physical or toy_fraction < float(config["minimum_accepted_toy_fraction"]):
        fit_status, reason = "fit_invalid", "physicality_or_toy_quality_gate_failed"
    elif covariance_degenerate or profile["status"] != "two_sided":
        fit_status, reason = "fit_valid_but_P0_weakly_constrained", "profile_or_covariance_identifiability_weak"
    else:
        fit_status, reason = "fit_valid", "all_direct_fit_gates_passed"
    return {
        "fit_status": fit_status,
        "reason": reason,
        "model_family": family,
        "fit_record_count": int(len(values)),
        "fit_parameters_log": list(map(float, parameters)),
        "mu": float(mu),
        "P0": float(p0),
        "P0_statistical_uncertainty": sigma_p0,
        "P0_relative_statistical_uncertainty": (float(sigma_p0 / p0) if sigma_p0 is not None and p0 > 0.0 else None),
        "P0_toy_interval_16_84": ([float(np.quantile(toy_p0, 0.16)), float(np.quantile(toy_p0, 0.84))] if toy_p0 else [None, None]),
        "kaon_efficiency": float(kaon_efficiency),
        "physical_control_efficiency": float(control_efficiency),
        "transfer": transfer,
        "transfer_statistical_uncertainty": sigma_transfer,
        "nll": float(nll(parameters)),
        "covariance_status": covariance_status,
        "covariance": covariance.tolist() if covariance is not None and np.all(np.isfinite(covariance)) else None,
        "covariance_condition_number": float(condition) if math.isfinite(condition) else None,
        "maximum_pair_correlation": max_correlation,
        "profile": profile,
        "accepted_toy_fraction": toy_fraction,
        "accepted_toy_count": int(len(toy_transfers)),
        "toy_count": int(config["toy_count"]),
        "toy_seed_context": deepcopy(config.get("_toy_context") or {}),
        "integer_like": compound_parameters is None,
        "response_mask": dict(response_mask),
    }


def fit_zero_photoelectron_response(values, *, response_mask, contract, config):
    """Fit a positive response sample with a model linked to a real P0 atom."""
    response_mask = normalize_hgcer_mask(response_mask, name="response sample mask")
    contract = dict(contract or fingerprint_hgcer_pid_contract())
    values = np.asarray([value for value in (_finite(item) for item in values) if value is not None], dtype=float)
    selected = values[np.asarray([hgcer_mask_accepts(value, response_mask) for value in values], dtype=bool)]
    if len(selected) < int(config["minimum_prompt_pion_records"]):
        return {"fit_status": "fit_invalid", "reason": "insufficient_prompt_pion_response_records", "fit_record_count": int(len(selected))}
    integer_like = bool(np.all(np.abs(selected - np.rint(selected)) <= float(config["integer_tolerance"])))
    if integer_like:
        primary = _fit_poisson(selected, response_mask, contract, config)
        alternate = _fit_negative_binomial(selected, response_mask, contract, config)
    else:
        primary = _fit_compound_poisson(selected, response_mask, contract, config)
        alternate = _fit_compound_poisson(selected, response_mask, contract, config, exponential_gain=True)
    primary["integer_like"] = integer_like
    primary["alternate_model_family"] = alternate.get("model_family")
    primary["alternate_fit_status"] = alternate.get("fit_status")
    primary["alternate_transfer"] = alternate.get("transfer")
    # The complete alternate result is diagnostic display data only.  Keeping
    # it beside the frozen primary fit lets the PDF draw the already-computed
    # alternate response curve without a second fit at render time.
    primary["alternate_fit"] = alternate
    if primary.get("transfer") is not None and alternate.get("fit_status") == "fit_valid" and alternate.get("transfer") is not None:
        primary["transfer_model_uncertainty"] = abs(float(primary["transfer"]) - float(alternate["transfer"]))
        primary["model_systematic_status"] = "available"
    else:
        primary["transfer_model_uncertainty"] = None
        primary["model_systematic_status"] = "unavailable"
    terms = [primary.get("transfer_statistical_uncertainty"), primary.get("transfer_model_uncertainty")]
    primary["transfer_total_uncertainty"] = math.sqrt(sum(float(term) ** 2 for term in terms if term is not None)) if any(term is not None for term in terms) else None
    return primary


def _audit_tree_records(records, mask):
    finite_values, violations = [], 0
    for record in records or ():
        value = _finite((record or {}).get("P_hgcer_npeSum"))
        if value is None:
            violations += 1
            continue
        finite_values.append(value)
        if not hgcer_mask_accepts(value, mask):
            violations += 1
    count = len(records or ())
    return {
        "record_count": count,
        "finite_record_count": len(finite_values),
        "violation_count": violations,
        "violation_fraction": float(violations) / float(count) if count else 1.0,
        "minimum_npe": min(finite_values) if finite_values else None,
        "maximum_npe": max(finite_values) if finite_values else None,
        "passed": bool(count > 0 and violations == 0),
    }


def audit_pion_hgcer_transfer_inputs(part1_payload, pion_control_cache, contract=None):
    """Fail closed for Part 2 only; this never filters or repairs records."""
    contract = fingerprint_hgcer_pid_contract() if contract is None else contract
    result = {
        "status": "available", "reason": None, "pid_contract": contract,
        "four_masks": {
            "S_K_tree": dict(contract["masks"]["kaon_tree"]),
            "S_pi_tree": dict(contract["masks"]["pion_tree"]),
            "S_pi_response_sample": {
                "mask": dict(contract["masks"]["pion_tree"]),
                "source_label": "prompt", "allcuts": True, "weighting": "unweighted",
                "rf_state": "noRF",
            },
            "S_pi_physical_control": dict(contract["masks"]["physical_pion_control"]),
        },
    }
    if not isinstance(part1_payload, dict) or part1_payload.get("status") != "available":
        result.update(status="unavailable", reason="part1_diagnostic_records_unavailable")
        return result
    provenance = part1_payload.get("source_provenance") or {}
    invalid_trees = []
    for side in ("kaon", "pion"):
        for label, entry in (provenance.get(side) or {}).items():
            name = str((entry or {}).get("tree_name") or "")
            if not name.endswith("_noRF"):
                invalid_trees.append("{}:{}={}".format(side, label, name))
    result["source_provenance"] = deepcopy(provenance)
    if invalid_trees:
        result.update(status="unavailable", reason="norf_provenance_failed", invalid_trees=invalid_trees)
        return result
    records = part1_payload.get("records") or {}
    result["event_content_audit"] = {
        "kaon": _audit_tree_records(records.get("kaon"), contract["masks"]["kaon_tree"]),
        "pion": _audit_tree_records(records.get("pion"), contract["masks"]["pion_tree"]),
    }
    if not all(entry["passed"] for entry in result["event_content_audit"].values()):
        result.update(status="unavailable", reason="pid_tree_event_content_audit_failed")
        return result
    cache_fingerprint = str((pion_control_cache or {}).get("physical_pion_control_mask_fingerprint") or "")
    result["cache_mask_fingerprint"] = cache_fingerprint
    if cache_fingerprint != contract["fingerprint"]:
        result.update(status="unavailable", reason="pion_control_cache_mask_contract_mismatch")
    return result


def _cell_key(cell):
    return int(cell["t_index"]), int(cell["delta_index"])


def _kinematic_summary(records):
    summary = {}
    for key in ("Q2", "W", "epsilon", "phi"):
        values = [_finite((record or {}).get(key)) for record in records or ()]
        values = [value for value in values if value is not None]
        summary[key] = {
            "record_count": len(values),
            "minimum": min(values) if values else None,
            "maximum": max(values) if values else None,
            "median": float(np.median(values)) if values else None,
        }
    return summary


def _unavailable_transfer(reason, config=None, audit=None):
    return {
        "status": "unavailable", "reason": str(reason),
        "diagnostic_label": pion_hgcer_transfer_display_text("unavailable"),
        "non_authoritative": True, "production_side_effect_free": True,
        "production_pion_subtraction_unchanged": True,
        "config": deepcopy(config or {}), "audit": deepcopy(audit or {}),
        "noRF_host_terminology": True, "rf_restoration_applied": False,
        "Part3_eligibility": "unavailable",
    }


def _response_display_curve(fit, response_mask, fit_range, sample_count):
    """Return a display-only expected response curve from an existing fit."""
    if not isinstance(fit, dict) or not fit.get("fit_parameters_log"):
        return None
    family = str(fit.get("model_family") or "")
    parameters = [float(value) for value in fit.get("fit_parameters_log") or ()]
    lower, upper = [float(value) for value in fit_range]
    mu = math.exp(parameters[0])
    if family in {"zero_truncated_poisson", "zero_truncated_negative_binomial"}:
        x_values = np.arange(max(0, int(math.ceil(lower))), int(math.floor(upper)) + 1, dtype=float)
        if family == "zero_truncated_poisson":
            probability = np.exp(_poisson_logpmf(x_values, mu))
            normalization = _poisson_fit_normalization(mu, response_mask, fit_range)
        else:
            shape = math.exp(parameters[1])
            probability = np.exp(_nbinom_logpmf(x_values, mu, shape))
            normalization = sum(
                math.exp(float(_nbinom_logpmf([number], mu, shape)[0]))
                for number in range(max(0, int(math.ceil(lower))), int(math.floor(upper)) + 1)
                if hgcer_mask_accepts(number, response_mask)
            )
        if normalization <= 0.0:
            return None
        return {
            "x": list(map(float, x_values)),
            "y": list(map(float, float(sample_count) * probability / normalization)),
            "kind": "expected_bin_count",
        }
    if family not in {
        "zero_truncated_compound_poisson_gamma_gain",
        "zero_truncated_compound_poisson_exponential_gain",
    }:
        return None
    shape = 1.0 if "exponential" in family else math.exp(parameters[1])
    scale = math.exp(parameters[-1])
    x_values = np.linspace(lower, upper, 240)
    normalization = _compound_interval_probability(mu, shape, scale, lower, upper)
    if normalization <= 0.0:
        return None
    # The fit domain is positive in the continuous case; the response mask
    # has already been enforced by the selected prompt-pion records.
    density = _compound_density(x_values, mu, shape, scale) / normalization
    return {
        "x": list(map(float, x_values)),
        "y": list(map(float, density * float(sample_count) * (upper - lower) / 40.0)),
        "kind": "expected_bin_count_density",
    }


def _make_response_display(values, fit, config):
    """Compact JSON-safe binned data and existing-fit curves for rendering."""
    lower, upper = [float(value) for value in config["fit_range"]]
    integer_like = bool((fit or {}).get("integer_like"))
    if integer_like:
        bin_edges = np.arange(max(-0.5, math.floor(lower) - 0.5), math.ceil(upper) + 1.5, 1.0)
    else:
        bin_edges = np.linspace(lower, upper, 41)
    selected = np.asarray([
        value for value in (_finite(item) for item in values or ())
        if value is not None and lower <= value <= upper
    ], dtype=float)
    counts, edges = np.histogram(selected, bins=bin_edges)
    response_mask = (fit or {}).get("response_mask") or config["pid_contract"]["masks"]["pion_tree"]
    display = {
        "fit_range": [lower, upper],
        "bin_edges": list(map(float, edges)),
        "bin_counts": list(map(float, counts)),
        "record_count": int(len(selected)),
        "primary_response": _response_display_curve(fit or {}, response_mask, (lower, upper), len(selected)),
        "alternate_response": None,
    }
    alternate = (fit or {}).get("alternate_fit") or {}
    if alternate.get("fit_status") in {"fit_valid", "fit_valid_but_P0_weakly_constrained"}:
        display["alternate_response"] = _response_display_curve(
            alternate, response_mask, (lower, upper), len(selected)
        )
    return _json_ready(display)


def _solution_display_fields(solution, fit):
    """Expose uncertainty components without modifying map construction."""
    result = dict(solution or {})
    source = result.get("solution_source")
    source_fit = result.get("pooled_fit") if source == "same_t_pooled" else fit
    if source == "direct" or source == "same_t_pooled":
        result["transfer_statistical_uncertainty"] = (source_fit or {}).get("transfer_statistical_uncertainty")
        result["transfer_model_uncertainty"] = (source_fit or {}).get("transfer_model_uncertainty")
        result["transfer_fallback_uncertainty"] = 0.0
    elif source in {"same_t_delta_bracket", "same_t_delta_edge"}:
        # Existing fallback totals are intentionally preserved.  Their extra
        # conservative term is shown separately rather than recomputed.
        result["transfer_statistical_uncertainty"] = None
        result["transfer_model_uncertainty"] = None
        result["transfer_fallback_uncertainty"] = result.get("transfer_total_uncertainty")
    else:
        result["transfer_statistical_uncertainty"] = None
        result["transfer_model_uncertainty"] = None
        result["transfer_fallback_uncertainty"] = None
    return result


def _cell_part3_eligibility(cell, fit, solution):
    if cell.get("support_class") != "supported":
        return {"status": "unsupported", "reason": "full_cell_not_supported"}
    if fit.get("fit_status") == "fit_valid_but_P0_weakly_constrained":
        return {"status": "weak_P0", "reason": fit.get("reason")}
    if solution.get("transfer") is None:
        return {"status": "unresolved", "reason": solution.get("solution_source")}
    if solution.get("solution_source") != "direct":
        return {"status": "review_required_fallback", "reason": solution.get("solution_source")}
    if fit.get("model_systematic_status") != "available":
        return {"status": "incomplete_model_uncertainty", "reason": "alternate_model_unavailable"}
    return {"status": "review_eligible", "reason": "direct_valid_model_checked"}


def _solve_fallbacks(cells, fits, response_samples, contract, config):
    solution = {}
    fallback = config["fallback"]
    supported = {key for key, cell in cells.items() if cell.get("support_class") == "supported"}
    for key, cell in cells.items():
        fit = fits.get(key) or {}
        if key not in supported:
            solution[key] = {"solution_source": "unsupported", "transfer": None, "transfer_total_uncertainty": None}
        elif fit.get("fit_status") == "fit_valid":
            solution[key] = {"solution_source": "direct", "transfer": fit.get("transfer"), "transfer_total_uncertainty": fit.get("transfer_total_uncertainty")}
        else:
            solution[key] = {"solution_source": "unresolved", "transfer": None, "transfer_total_uncertainty": None}
    for (t_index, delta_index), current in list(solution.items()):
        if current["solution_source"] != "unresolved":
            continue
        direct = sorted(
            (key, value) for key, value in solution.items()
            if key[0] == t_index and value["solution_source"] == "direct"
        )
        lower = [(key, value) for key, value in direct if key[1] < delta_index]
        upper = [(key, value) for key, value in direct if key[1] > delta_index]
        if lower and upper:
            left_key, left = lower[-1]
            right_key, right = upper[0]
            fraction = float(delta_index - left_key[1]) / float(right_key[1] - left_key[1])
            transfer = (1.0 - fraction) * float(left["transfer"]) + fraction * float(right["transfer"])
            base = max(float(left.get("transfer_total_uncertainty") or 0.0), float(right.get("transfer_total_uncertainty") or 0.0))
            solution[(t_index, delta_index)] = {"solution_source": "same_t_delta_bracket", "transfer": transfer, "transfer_total_uncertainty": base * (1.0 + fallback["bracketing_relative_uncertainty"])}
            continue
        # Edge only: a single direct neighbour cannot seed an interior gap.
        ordered_delta = sorted(key[1] for key in cells if key[0] == t_index)
        at_edge = delta_index in (ordered_delta[0], ordered_delta[-1]) if ordered_delta else False
        neighbours = lower[-1:] + upper[:1]
        if at_edge and len(neighbours) == 1:
            _, neighbour = neighbours[0]
            solution[(t_index, delta_index)] = {"solution_source": "same_t_delta_edge", "transfer": neighbour["transfer"], "transfer_total_uncertainty": float(neighbour.get("transfer_total_uncertainty") or 0.0) * (1.0 + fallback["edge_relative_uncertainty"])}
            continue
        pooled_values = []
        for candidate_key, values in response_samples.items():
            if candidate_key[0] == t_index and candidate_key in supported:
                pooled_values.extend(values)
        if len(pooled_values) >= int(fallback["pooled_minimum_prompt_pion_records"]):
            pooled_config = deepcopy(config)
            pooled_config["_toy_context"] = {"t_index": int(t_index), "delta_index": "pooled"}
            pooled_fit = fit_zero_photoelectron_response(
                pooled_values, response_mask=contract["masks"]["pion_tree"], contract=contract, config=pooled_config,
            )
            if pooled_fit.get("fit_status") == "fit_valid":
                pooled_fit["response_display"] = _make_response_display(
                    pooled_values, pooled_fit, config
                )
                solution[(t_index, delta_index)] = {"solution_source": "same_t_pooled", "transfer": pooled_fit.get("transfer"), "transfer_total_uncertainty": pooled_fit.get("transfer_total_uncertainty"), "pooled_fit": pooled_fit}
    return solution


def build_pion_hgcer_zerope_transfer_map(part1_payload, pion_control_cache, *, config):
    """Build and freeze the Part-2 map without touching production products."""
    config = deepcopy(config or {})
    contract = dict(config.get("pid_contract") or fingerprint_hgcer_pid_contract())
    audit = audit_pion_hgcer_transfer_inputs(part1_payload, pion_control_cache, contract)
    if audit.get("status") != "available":
        return _unavailable_transfer(audit.get("reason"), config, audit)
    cells = {_cell_key(cell): dict(cell) for cell in (part1_payload.get("cells") or ())}
    if not cells:
        return _unavailable_transfer("part1_cells_unavailable", config, audit)
    response_samples = {key: [] for key in cells}
    for record in (part1_payload.get("records") or {}).get("pion") or ():
        if str(record.get("source_label")) != "prompt" or not bool(record.get("allcuts")):
            continue
        key = (int(record.get("canonical_t_index")), int(record.get("delta_index")))
        value = _finite(record.get("P_hgcer_npeSum"))
        if key in response_samples and value is not None:
            response_samples[key].append(value)
    fits = {}
    for key, cell in cells.items():
        if cell.get("support_class") != "supported":
            fits[key] = {"fit_status": "fit_invalid", "reason": "full_cell_not_supported", "fit_record_count": len(response_samples[key])}
            continue
        direct_config = deepcopy(config)
        direct_config["_toy_context"] = {"t_index": int(key[0]), "delta_index": int(key[1])}
        fits[key] = fit_zero_photoelectron_response(
            response_samples[key], response_mask=contract["masks"]["pion_tree"], contract=contract, config=direct_config,
        )
        fits[key]["response_display"] = _make_response_display(
            response_samples[key], fits[key], config
        )
    solutions = _solve_fallbacks(cells, fits, response_samples, contract, config)
    serialized_cells = []
    for key in sorted(cells):
        entry = dict(cells[key])
        entry["response_prompt_record_count"] = int(len(response_samples[key]))
        entry["fit"] = fits[key]
        entry["solution"] = _solution_display_fields(solutions[key], fits[key])
        entry["part3_cell_eligibility"] = _cell_part3_eligibility(
            entry, fits[key], entry["solution"]
        )
        serialized_cells.append(entry)
    essentials = [{
        "t_index": cell["t_index"], "delta_index": cell["delta_index"],
        "support_class": cell.get("support_class"), "fit_status": cell["fit"].get("fit_status"),
        "transfer": cell["solution"].get("transfer"), "solution_source": cell["solution"].get("solution_source"),
    } for cell in serialized_cells]
    all_records = []
    for side in ("kaon", "pion"):
        all_records.extend((part1_payload.get("records") or {}).get(side) or ())
    kinematics = _kinematic_summary(all_records)
    map_fingerprint = _fingerprint({
        "t_edges": part1_payload.get("t_edges"), "delta_edges": part1_payload.get("delta_edges"),
        "coordinate_fingerprint": part1_payload.get("coordinate_fingerprint"),
        "pid_contract": contract, "response_config": config, "kinematics": kinematics,
        "cells": essentials,
    })
    supported_cells = [cell for cell in serialized_cells if cell.get("support_class") == "supported"]
    eligible = bool(supported_cells) and all(
        cell["solution"].get("transfer") is not None
        and cell["fit"].get("model_systematic_status") == "available"
        for cell in supported_cells
    )
    return {
        "status": "available", "diagnostic_label": pion_hgcer_transfer_display_text("part2"),
        "non_authoritative": True, "production_side_effect_free": True,
        "production_pion_subtraction_unchanged": True,
        "noRF_host_terminology": True, "rf_restoration_applied": False,
        "config": config, "pid_contract": contract, "audit": audit,
        "coordinate_fingerprint": part1_payload.get("coordinate_fingerprint"),
        "coordinate_contract": deepcopy(part1_payload.get("coordinate_contract") or {}),
        "kinematics": kinematics,
        "t_edges": list(part1_payload.get("t_edges") or ()),
        "delta_edges": list(part1_payload.get("delta_edges") or ()),
        "cells": tuple(serialized_cells), "map_fingerprint": map_fingerprint,
        "frozen": True,
        "Part3_eligibility": "review_eligible" if eligible else "review_required_incomplete_uncertainty_or_solution",
    }


_PART2_NUMERIC_MATRIX_PATHS = {
    "P0": ("fit", "P0"),
    "relative_P0_uncertainty": ("fit", "P0_relative_statistical_uncertainty"),
    "transfer": ("solution", "transfer"),
    "statistical_uncertainty": ("solution", "transfer_statistical_uncertainty"),
    "model_uncertainty": ("solution", "transfer_model_uncertainty"),
    "fallback_uncertainty": ("solution", "transfer_fallback_uncertainty"),
    "total_uncertainty": ("solution", "transfer_total_uncertainty"),
}


_PART2_CATEGORICAL_MATRIX_PATHS = {
    "response_family": ("fit", "model_family"),
    "fit_status": ("fit", "fit_status"),
    "solution_source": ("solution", "solution_source"),
    "part3_eligibility": ("part3_cell_eligibility", "status"),
}


def _cell_path_value(cell, path):
    value = cell
    for key in path:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def extract_pion_hgcer_transfer_tdelta_matrix(payload, field):
    """Return a masked t-delta display matrix; undefined values stay None."""
    if field not in _PART2_NUMERIC_MATRIX_PATHS:
        raise ValueError("unknown Part-2 numeric matrix '{}'".format(field))
    t_edges = list((payload or {}).get("t_edges") or ())
    delta_edges = list((payload or {}).get("delta_edges") or ())
    matrix = [[None for _ in range(max(0, len(delta_edges) - 1))] for _ in range(max(0, len(t_edges) - 1))]
    for cell in (payload or {}).get("cells") or ():
        try:
            t_index, delta_index = int(cell["t_index"]), int(cell["delta_index"])
        except (KeyError, TypeError, ValueError):
            continue
        if not (0 <= t_index < len(matrix) and 0 <= delta_index < len(matrix[t_index])):
            continue
        matrix[t_index][delta_index] = _finite(_cell_path_value(cell, _PART2_NUMERIC_MATRIX_PATHS[field]))
    return {"field": field, "t_edges": t_edges, "delta_edges": delta_edges, "values": matrix}


def extract_pion_hgcer_transfer_categorical_map(payload, field):
    """Return an exact categorical t-delta display matrix with masked cells."""
    if field not in _PART2_CATEGORICAL_MATRIX_PATHS:
        raise ValueError("unknown Part-2 categorical matrix '{}'".format(field))
    t_edges = list((payload or {}).get("t_edges") or ())
    delta_edges = list((payload or {}).get("delta_edges") or ())
    matrix = [[None for _ in range(max(0, len(delta_edges) - 1))] for _ in range(max(0, len(t_edges) - 1))]
    for cell in (payload or {}).get("cells") or ():
        try:
            t_index, delta_index = int(cell["t_index"]), int(cell["delta_index"])
        except (KeyError, TypeError, ValueError):
            continue
        if 0 <= t_index < len(matrix) and 0 <= delta_index < len(matrix[t_index]):
            value = _cell_path_value(cell, _PART2_CATEGORICAL_MATRIX_PATHS[field])
            matrix[t_index][delta_index] = str(value) if value is not None else None
    categories = sorted({value for row in matrix for value in row if value is not None})
    return {"field": field, "t_edges": t_edges, "delta_edges": delta_edges, "values": matrix, "categories": categories}


def extract_pion_hgcer_transfer_t_series(payload, t_index, field):
    """Return a masked per-t delta series from the frozen cell payload."""
    matrix = extract_pion_hgcer_transfer_tdelta_matrix(payload, field)
    index = int(t_index)
    if not 0 <= index < len(matrix["values"]):
        raise ValueError("Part-2 t index {} is outside the frozen map".format(index))
    edges = matrix["delta_edges"]
    return {
        "field": field, "t_index": index,
        "delta_centers": [0.5 * (float(edges[i]) + float(edges[i + 1])) for i in range(len(edges) - 1)],
        "delta_half_widths": [0.5 * (float(edges[i + 1]) - float(edges[i])) for i in range(len(edges) - 1)],
        "values": list(matrix["values"][index]),
    }


def extract_pion_hgcer_response_fit_pages(payload):
    """List direct and same-t pooled existing fits once, without fitting."""
    pages, pooled_t = [], set()
    for cell in sorted((payload or {}).get("cells") or (), key=lambda item: (item.get("t_index", -1), item.get("delta_index", -1))):
        fit = cell.get("fit") or {}
        display = fit.get("response_display") or {}
        if display.get("primary_response") is not None:
            pages.append({
                "fit_source": "direct", "t_index": int(cell["t_index"]),
                "delta_index": int(cell["delta_index"]), "fit": fit,
                "solution": cell.get("solution") or {},
                "eligibility": cell.get("part3_cell_eligibility") or {},
            })
        pooled = (cell.get("solution") or {}).get("pooled_fit") or {}
        pooled_display = pooled.get("response_display") or {}
        t_index = int(cell.get("t_index", -1))
        if t_index not in pooled_t and pooled_display.get("primary_response") is not None:
            pages.append({"fit_source": "same_t_pooled", "t_index": t_index, "delta_index": None, "fit": pooled, "solution": cell.get("solution") or {}, "eligibility": cell.get("part3_cell_eligibility") or {}})
            pooled_t.add(t_index)
    return pages


def _edge_index(value, edges):
    value = _finite(value)
    if value is None or len(edges or ()) < 2:
        return None
    numeric_edges = [float(edge) for edge in edges]
    if value < numeric_edges[0] or value > numeric_edges[-1]:
        return None
    if value == numeric_edges[-1]:
        return len(numeric_edges) - 2
    index = int(np.searchsorted(numeric_edges, value, side="right") - 1)
    return index if 0 <= index < len(numeric_edges) - 1 else None


def _histograms_match(left, right, tolerance=1.0e-10):
    if left is None or right is None:
        return False, "histogram_missing"
    if int(left.GetNbinsX()) != int(right.GetNbinsX()):
        return False, "bin_count_mismatch"
    for index in range(0, int(left.GetNbinsX()) + 2):
        if abs(float(left.GetBinContent(index)) - float(right.GetBinContent(index))) > tolerance:
            return False, "content_mismatch_bin{}".format(index)
        if abs(float(left.GetBinError(index)) - float(right.GetBinError(index))) > tolerance:
            return False, "error_mismatch_bin{}".format(index)
    return True, None


def _control_source_accounting_matches(section, source_summary):
    """Check the all-cuts source algebra used by a cached control projection."""
    accounting = (section or {}).get("source_accounting") or {}
    for label, observed in (source_summary or {}).items():
        expected = accounting.get(label) or {}
        if int(expected.get("allcuts_records", -1)) != int(observed.get("records", 0)):
            return False, "source_allcuts_record_count:{}".format(label)
        expected_sum = _finite(expected.get("signed_weight_sum"))
        if expected_sum is None or abs(expected_sum - float(observed.get("signed_sum", 0.0))) > 1.0e-10:
            return False, "source_signed_algebra:{}".format(label)
        expected_coefficient = _finite(expected.get("coefficient"))
        if expected_coefficient is None or expected_coefficient != _finite(observed.get("coefficient")):
            return False, "source_coefficient:{}".format(label)
    return True, None


def inspect_pion_hgcer_control_population(payload, pion_control_cache):
    """Pure-Python proof of the record contract behind H_pion_control_cut."""
    result = {"status": "available", "reason": None, "by_t": {}}
    contract, cache = (payload or {}).get("pid_contract") or {}, pion_control_cache or {}
    mask = ((contract.get("masks") or {}).get("physical_pion_control") or {})
    fingerprint = str(contract.get("fingerprint") or "")
    if not mask or not fingerprint or str(cache.get("physical_pion_control_mask_fingerprint") or "") != fingerprint:
        return {**result, "status": "unavailable", "reason": "physical_control_mask_fingerprint_mismatch"}
    if normalize_hgcer_mask(cache.get("physical_pion_control_mask") or {}, name="cache physical pion control") != normalize_hgcer_mask(mask, name="Part-2 physical pion control"):
        return {**result, "status": "unavailable", "reason": "physical_control_mask_expression_mismatch"}
    t_edges, delta_edges = list((payload or {}).get("t_edges") or ()), list((payload or {}).get("delta_edges") or ())
    coordinate_fingerprint = str((payload or {}).get("coordinate_fingerprint") or "")
    if not coordinate_fingerprint or str(cache.get("coordinate_fingerprint") or "") != coordinate_fingerprint:
        return {**result, "status": "unavailable", "reason": "control_coordinate_fingerprint_mismatch"}
    cached_t = list(cache.get("by_t") or ())
    if len(cached_t) != max(0, len(t_edges) - 1):
        return {**result, "status": "unavailable", "reason": "control_canonical_t_count_mismatch"}
    solution_lookup = {(int(cell["t_index"]), int(cell["delta_index"])): (cell.get("solution") or {}) for cell in (payload or {}).get("cells") or ()}
    for t_index, section in enumerate(cached_t):
        invalid, source_summary = [], {}
        for record in (section or {}).get("records") or ():
            if not bool(record.get("allcuts")):
                continue
            label, coefficient = str(record.get("source_label") or "unknown"), _finite(record.get("coefficient"))
            source = (cache.get("source_accounting") or {}).get(label) or {}
            if coefficient is None or coefficient != _finite(source.get("coefficient")):
                invalid.append("source_coefficient:{}".format(label)); continue
            source_summary.setdefault(label, {"records": 0, "signed_sum": 0.0, "coefficient": coefficient})
            source_summary[label]["records"] += 1; source_summary[label]["signed_sum"] += coefficient
            try:
                record_delta_index = int(record.get("delta_index"))
            except (TypeError, ValueError):
                record_delta_index = None
            delta_index = _edge_index(record.get("ssdelta"), delta_edges)
            try:
                record_t_index = int(record.get("t_index"))
            except (TypeError, ValueError):
                record_t_index = None
            if str(record.get("rf_state") or "") != "noRF": invalid.append("norf:{}".format(label))
            if str(record.get("coordinate_fingerprint") or "") != coordinate_fingerprint: invalid.append("coordinate_fingerprint:{}".format(label))
            if not hgcer_mask_accepts(record.get("P_hgcer_npeSum"), mask): invalid.append("physical_control_mask:{}".format(label))
            if record_t_index != t_index or _edge_index(record.get("adj_t"), t_edges) != t_index: invalid.append("canonical_t_assignment:{}".format(label))
            if delta_index is None or record_delta_index != delta_index: invalid.append("delta_assignment:{}".format(label))
            if ((solution_lookup.get((t_index, delta_index)) or {}).get("transfer")) is None: invalid.append("unresolved_transfer:{}".format(label))
        accounting_ok, accounting_reason = _control_source_accounting_matches(section, source_summary)
        if not accounting_ok:
            invalid.append(accounting_reason)
        result["by_t"][t_index] = {"status": "available" if not invalid else "unavailable", "reason": invalid[0] if invalid else None, "source_summary": source_summary}
    if any(entry.get("status") != "available" for entry in result["by_t"].values()):
        result.update(status="unavailable", reason="control_record_contract_failed")
    return result


def audit_pion_hgcer_control_population(payload, pion_control_cache):
    """Prove that the display control spectrum is the Part-2 source population.

    This is deliberately a read-only, diagnostic audit.  It reconstructs a
    detached signed projection from the exact cached all-cuts records and
    never changes cache records, map cells, or proposed weights.
    """
    inspection = inspect_pion_hgcer_control_population(payload, pion_control_cache)
    if inspection.get("status") != "available":
        return {**inspection, "global": None}
    result = {"status": "available", "reason": None, "by_t": {}, "global": None}
    contract = (payload or {}).get("pid_contract") or {}
    mask = ((contract.get("masks") or {}).get("physical_pion_control") or {})
    fingerprint = str(contract.get("fingerprint") or "")
    cache = pion_control_cache or {}
    if not mask or not fingerprint or str(cache.get("physical_pion_control_mask_fingerprint") or "") != fingerprint:
        result.update(status="unavailable", reason="physical_control_mask_fingerprint_mismatch")
        return result
    if normalize_hgcer_mask(cache.get("physical_pion_control_mask") or {}, name="cache physical pion control") != normalize_hgcer_mask(mask, name="Part-2 physical pion control"):
        result.update(status="unavailable", reason="physical_control_mask_expression_mismatch")
        return result
    t_edges = list((payload or {}).get("t_edges") or ())
    delta_edges = list((payload or {}).get("delta_edges") or ())
    coordinate_fingerprint = str((payload or {}).get("coordinate_fingerprint") or "")
    if not coordinate_fingerprint or str(cache.get("coordinate_fingerprint") or "") != coordinate_fingerprint:
        result.update(status="unavailable", reason="control_coordinate_fingerprint_mismatch")
        return result
    solution_lookup = {
        (int(cell["t_index"]), int(cell["delta_index"])): (cell.get("solution") or {})
        for cell in (payload or {}).get("cells") or ()
    }
    cached_t = list(cache.get("by_t") or ())
    if len(cached_t) != max(0, len(t_edges) - 1):
        result.update(status="unavailable", reason="control_canonical_t_count_mismatch")
        return result
    if ROOT is None or clone_root_histogram is None:
        result.update(status="unavailable", reason="PyROOT_unavailable_for_control_projection")
        return result
    global_before = None
    for t_index, section in enumerate(cached_t):
        reference = (section or {}).get("H_pion_control_cut")
        if reference is None:
            result.update(status="unavailable", reason="control_cut_histogram_missing:t{}".format(t_index + 1))
            return result
        projection = clone_root_histogram(
            reference, scope="pion_hgcer_part2", role="verified_control_before_transfer",
            name="H_MM_part2_control_verified_t{}".format(t_index + 1), reset=True, sumw2=True,
        )
        source_summary, invalid = {}, []
        for record in (section or {}).get("records") or ():
            if not bool(record.get("allcuts")):
                continue
            label = str(record.get("source_label") or "unknown")
            coefficient = _finite(record.get("coefficient"))
            source = (cache.get("source_accounting") or {}).get(label) or {}
            if coefficient is None or coefficient != _finite(source.get("coefficient")):
                invalid.append("source_coefficient:{}".format(label))
                continue
            source_summary.setdefault(label, {"records": 0, "signed_sum": 0.0, "coefficient": coefficient})
            source_summary[label]["records"] += 1
            source_summary[label]["signed_sum"] += coefficient
            if str(record.get("rf_state") or "") != "noRF":
                invalid.append("norf:{}".format(label))
            if str(record.get("coordinate_fingerprint") or "") != coordinate_fingerprint:
                invalid.append("coordinate_fingerprint:{}".format(label))
            if not hgcer_mask_accepts(record.get("P_hgcer_npeSum"), mask):
                invalid.append("physical_control_mask:{}".format(label))
            try:
                record_t_index = int(record.get("t_index"))
            except (TypeError, ValueError):
                record_t_index = None
            if record_t_index != t_index or _edge_index(record.get("adj_t"), t_edges) != t_index:
                invalid.append("canonical_t_assignment:{}".format(label))
            delta_index = _edge_index(record.get("ssdelta"), delta_edges)
            try:
                record_delta_index = int(record.get("delta_index"))
            except (TypeError, ValueError):
                record_delta_index = None
            if delta_index is None or record_delta_index != delta_index:
                invalid.append("delta_assignment:{}".format(label))
            solution = solution_lookup.get((t_index, delta_index)) or {}
            if solution.get("transfer") is None:
                invalid.append("unresolved_transfer:{}".format(label))
            projection.Fill(float(record["adj_MM"]), coefficient)
        accounting_ok, accounting_reason = _control_source_accounting_matches(section, source_summary)
        if not accounting_ok:
            invalid.append(accounting_reason)
        matched, mismatch_reason = _histograms_match(projection, reference)
        if invalid or not matched:
            result["by_t"][t_index] = {
                "status": "unavailable", "reason": (invalid[0] if invalid else mismatch_reason),
                "source_summary": source_summary,
            }
            continue
        before_integral = float(projection.Integral())
        result["by_t"][t_index] = {
            "status": "available", "before": projection,
            "source_summary": source_summary, "before_integral": before_integral,
            "physical_control_mask_fingerprint": fingerprint,
        }
        if global_before is None:
            global_before = clone_root_histogram(
                projection, scope="pion_hgcer_part2", role="verified_global_control_before_transfer",
                name="H_MM_part2_control_verified_global", reset=True, sumw2=True,
            )
        global_before.Add(projection)
    unavailable = [entry for entry in result["by_t"].values() if entry.get("status") != "available"]
    if unavailable:
        result.update(status="unavailable", reason="control_population_audit_failed")
        return result
    result["global"] = {
        "status": "available", "before": global_before,
        "before_integral": float(global_before.Integral()) if global_before is not None else None,
        "construction": "strict_sum_of_verified_per_t_control_histograms",
    }
    return result


def extract_pion_hgcer_transfer_series_detail(payload, t_index, field):
    """Return a series with uncertainty and source metadata, never zero-filled."""
    index = int(t_index)
    base = extract_pion_hgcer_transfer_t_series(payload, index, "transfer" if field == "transfer" else "P0")
    cells = {
        int(cell.get("delta_index", -1)): cell
        for cell in (payload or {}).get("cells") or () if int(cell.get("t_index", -1)) == index
    }
    values = []
    for delta_index, center in enumerate(base["delta_centers"]):
        cell = cells.get(delta_index) or {}
        fit, solution = cell.get("fit") or {}, cell.get("solution") or {}
        if field == "transfer":
            values.append({
                "delta_index": delta_index, "delta": center, "value": _finite(solution.get("transfer")),
                "statistical": _finite(solution.get("transfer_statistical_uncertainty")),
                "model": _finite(solution.get("transfer_model_uncertainty")),
                "fallback": _finite(solution.get("transfer_fallback_uncertainty")),
                "total": _finite(solution.get("transfer_total_uncertainty")),
                "source": solution.get("solution_source"),
            })
        else:
            direct_value = _finite(fit.get("P0")) if solution.get("solution_source") == "direct" else None
            values.append({
                "delta_index": delta_index, "delta": center, "value": direct_value,
                "statistical": _finite(fit.get("P0_statistical_uncertainty")) if direct_value is not None else None,
                "source": "direct" if direct_value is not None else None,
            })
    pooled = []
    seen = set()
    for cell in cells.values():
        fit = (cell.get("solution") or {}).get("pooled_fit") or {}
        key = _fingerprint({"mu": fit.get("mu"), "P0": fit.get("P0"), "source": fit.get("toy_seed_context")}) if fit else None
        if key and key not in seen and _finite(fit.get("P0")) is not None:
            pooled.append({"P0": _finite(fit.get("P0")), "statistical": _finite(fit.get("P0_statistical_uncertainty")), "source": "same_t_pooled"})
            seen.add(key)
    return {**base, "field": field, "points": values, "pooled": pooled}


def serialize_pion_hgcer_zerope_transfer(payload):
    payload = payload or {}
    result = {
        key: value for key, value in payload.items()
        if key not in {"histograms", "application_histograms", "application"}
    }
    application = payload.get("application")
    if isinstance(application, dict):
        result["application"] = {
            key: value for key, value in application.items()
            if key not in {"histograms", "application_histograms", "global_histograms"}
        }
    return _json_ready(result)


def write_pion_hgcer_zerope_transfer_json(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(serialize_pion_hgcer_zerope_transfer(payload), handle, sort_keys=True, indent=2, allow_nan=False)
    return path


def apply_frozen_pion_hgcer_transfer_map(payload, pion_control_cache, proton_t_products):
    """Build detached, proposed MM products from proton-cleaned noRF hosts."""
    if not isinstance(payload, dict) or payload.get("status") != "available" or not payload.get("frozen"):
        return {"status": "unavailable", "reason": "frozen_part2_map_unavailable"}
    if ROOT is None or clone_root_histogram is None:
        return {"status": "unavailable", "reason": "PyROOT_unavailable_for_proposed_MM"}
    solution_lookup = {
        (int(cell["t_index"]), int(cell["delta_index"])): cell.get("solution") or {}
        for cell in payload.get("cells") or ()
    }
    entries, detached = [], {}
    strict_host = strict_pion = strict_clean = 0.0
    for t_index, product in enumerate(proton_t_products or ()):
        no_rf_targets = (product or {}).get("cleaned_targets_pre_rf") or {}
        host = no_rf_targets.get("h_mm_nosub")
        if host is None:
            return {"status": "unavailable", "reason": "proton_cleaned_noRF_host_missing:t{}".format(t_index + 1)}
        pion_hist = clone_root_histogram(host, scope="pion_hgcer_part2", role="proposed_pion_mm", name="H_MM_part2_pion_noRF_t{}".format(t_index + 1), reset=True, sumw2=True)
        clean_hist = clone_root_histogram(host, scope="pion_hgcer_part2", role="proposed_clean_mm", name="H_MM_part2_kaon_clean_noRF_t{}".format(t_index + 1), reset=False, sumw2=True)
        applied, unresolved = 0, 0
        for record in ((pion_control_cache or {}).get("by_t") or ())[t_index].get("records") or ():
            if not bool(record.get("allcuts")):
                continue
            solution = solution_lookup.get((t_index, int(record.get("delta_index", -1)))) or {}
            transfer = solution.get("transfer")
            if transfer is None:
                unresolved += 1
                continue
            pion_hist.Fill(float(record["adj_MM"]), float(record["coefficient"]) * float(transfer))
            applied += 1
        clean_hist.Add(pion_hist, -1.0)
        host_integral, pion_integral, clean_integral = float(host.Integral()), float(pion_hist.Integral()), float(clean_hist.Integral())
        bin_count = int(host.GetNbinsX())
        affected_bins = 0
        pion_exceeds_host_bins = 0
        maximum_local_contamination = 0.0
        most_negative_clean_bin = None
        for bin_index in range(1, bin_count + 1):
            host_value = float(host.GetBinContent(bin_index))
            pion_value = float(pion_hist.GetBinContent(bin_index))
            clean_value = float(clean_hist.GetBinContent(bin_index))
            if pion_value != 0.0:
                affected_bins += 1
            if pion_value > host_value:
                pion_exceeds_host_bins += 1
            if host_value != 0.0:
                maximum_local_contamination = max(
                    maximum_local_contamination, abs(pion_value / host_value)
                )
            elif pion_value != 0.0:
                maximum_local_contamination = float("inf")
            if most_negative_clean_bin is None or clean_value < most_negative_clean_bin["content"]:
                most_negative_clean_bin = {
                    "bin": bin_index, "center": float(clean_hist.GetBinCenter(bin_index)),
                    "content": clean_value,
                }
        strict_host += host_integral; strict_pion += pion_integral; strict_clean += clean_integral
        entries.append({
            "t_index": t_index, "host_label": "proton-cleaned noRF host", "applied_records": applied,
            "unresolved_records": unresolved, "host_integral": host_integral,
            "pion_integral": pion_integral, "clean_integral": clean_integral,
            "pion_to_host": pion_integral / host_integral if host_integral else None,
            "pion_exceeds_host_bin_count": pion_exceeds_host_bins,
            "affected_bin_fraction": float(affected_bins) / float(max(1, bin_count)),
            "maximum_local_contamination": maximum_local_contamination,
            "most_negative_clean_bin": most_negative_clean_bin,
        })
        detached[t_index] = {"host": host, "pion": pion_hist, "clean": clean_hist}
    global_histograms = {}
    if detached:
        # These are sums of detached per-t products only.  They deliberately
        # never refill from records and therefore cannot change application
        # weights or hide negative proposed-clean bins.
        first = detached[sorted(detached)[0]]
        for role in ("host", "pion", "clean"):
            aggregate = clone_root_histogram(
                first[role], scope="pion_hgcer_part2", role="global_{}_mm".format(role),
                name="H_MM_part2_{}_noRF_global".format(role), reset=True, sumw2=True,
            )
            for product in detached.values():
                aggregate.Add(product[role])
            global_histograms[role] = aggregate
        global_closure = (
            float(global_histograms["host"].Integral())
            - float(global_histograms["pion"].Integral())
            - float(global_histograms["clean"].Integral())
        )
    else:
        global_closure = None
    return {
        "status": "available", "host_label": "proton-cleaned noRF host",
        "rf_restoration_applied": False, "t_products": entries, "histograms": detached,
        "global_histograms": global_histograms,
        "strict_global_sums": {
            "host": strict_host, "pion": strict_pion, "clean": strict_clean,
            "closure": strict_host - strict_pion - strict_clean,
            "detached_histogram_closure": global_closure,
        },
    }


PART2_PAGE_SPECS = {
    "unavailable": {"kind": "text", "renderer": "status", "required_roles": ()},
    "mask_audit": {"kind": "text", "renderer": "audit", "required_roles": ()},
    "response_family_audit": {"kind": "text", "renderer": "audit", "required_roles": ()},
    "response_family_map": {"kind": "graphical", "renderer": "categorical_map", "required_roles": ("response_family_map",)},
    "p0_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("P0_map",)},
    "relative_p0_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("relative_P0_map",)},
    "transfer_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("transfer_map",)},
    "statistical_uncertainty_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("statistical_uncertainty_map",)},
    "model_uncertainty_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("model_uncertainty_map",)},
    "fallback_uncertainty_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("fallback_uncertainty_map",)},
    "total_uncertainty_map": {"kind": "graphical", "renderer": "numeric_map", "required_roles": ("total_uncertainty_map",)},
    "fit_status_map": {"kind": "graphical", "renderer": "categorical_map", "required_roles": ("fit_status_map",)},
    "solution_source_map": {"kind": "graphical", "renderer": "categorical_map", "required_roles": ("solution_source_map",)},
    "part3_eligibility_map": {"kind": "graphical", "renderer": "categorical_map", "required_roles": ("part3_eligibility_map",)},
    "transfer_series": {"kind": "graphical", "renderer": "series", "required_roles": ("transfer_points", "statistical_errors", "total_uncertainty", "solution_source_legend")},
    "p0_series": {"kind": "graphical", "renderer": "series", "required_roles": ("P0_points", "P0_statistical_errors")},
    "response_fit": {"kind": "graphical", "renderer": "response_fit", "required_roles": ("prompt_pion_data", "primary_response", "fit_metadata", "profile_likelihood", "physical_control_threshold", "zero_pe_annotation")},
    "control_transfer": {"kind": "graphical", "renderer": "control_transfer", "required_roles": ("pion_control_mm", "transferred_pion_mm", "suppression_summary")},
    "proposed_mm": {"kind": "graphical", "renderer": "proposed_mm", "required_roles": ("host_mm", "pion_mm", "clean_mm", "pion_zoom", "summary_text")},
    "normalization_summary": {"kind": "graphical", "renderer": "normalization_summary", "required_roles": ("pion_control_integral_series", "transferred_pion_integral_series", "normalization_summary")},
    "simc_closure": {"kind": "graphical", "renderer": "simc_closure", "required_roles": ("proposed_pion_mm", "simc_reference")},
    "signal_closure": {"kind": "graphical", "renderer": "signal_closure", "required_roles": ("proposed_clean_mm", "KLambda_reference", "normalization_status")},
    "final_summary": {"kind": "graphical", "renderer": "final_summary", "required_roles": ("part2_summary_table",)},
    "closure_status": {"kind": "text", "renderer": "status", "required_roles": ()},
}


_PART2_RENDER_SERIAL = 0


def _part2_page(page_id, spec_key, *, t_index=None, status="ready", unavailable_reason=None, **detail):
    spec = PART2_PAGE_SPECS[spec_key]
    page = {
        "page_id": str(page_id), "spec_key": str(spec_key), "t_index": t_index,
        "page_kind": spec["kind"], "renderer": spec["renderer"],
        "required_roles": list(spec["required_roles"]), "status": status,
        "unavailable_reason": unavailable_reason,
    }
    page.update(detail)
    return page


def _page_has_defined_series(payload, t_index, field):
    return any(value is not None for value in extract_pion_hgcer_transfer_t_series(payload, t_index, field)["values"])


def _application_histograms_for_t(application, t_index):
    histograms = (application or {}).get("histograms") or {}
    return histograms.get(t_index) or histograms.get(str(t_index)) or {}


def _closure_group(renderer_inputs, t_index, kind):
    if t_index is None:
        group = (renderer_inputs or {}).get("global") or {}
        return (group.get(kind) or {}) if isinstance(group, dict) else {}
    by_t = (renderer_inputs or {}).get("by_t") or {}
    group = by_t.get(t_index) or by_t.get(str(t_index)) or {}
    return (group.get(kind) or {}) if isinstance(group, dict) else {}


def _closure_hist(value):
    return value.get("hist") if isinstance(value, dict) else value


def _closure_available(renderer_inputs, t_index, kind):
    group = _closure_group(renderer_inputs, t_index, kind)
    if kind == "signal":
        return _closure_hist(group.get("K-Lambda")) is not None
    return any(_closure_hist(value) is not None for value in group.values())


def expected_pion_hgcer_transfer_page_manifest(payload, *, renderer_inputs=None):
    """Build a runtime-sized Part-2 manifest without touching ROOT objects."""
    if not isinstance(payload, dict) or payload.get("status") != "available":
        return [_part2_page("part2_unavailable", "unavailable", status="unavailable", unavailable_reason=(payload or {}).get("reason", "unknown"))]
    pages = [
        _part2_page("part2_mask_audit", "mask_audit"),
        _part2_page("part2_response_family_audit", "response_family_audit"),
        _part2_page("part2_response_family_map", "response_family_map"),
        _part2_page("part2_p0_map", "p0_map"),
        _part2_page("part2_relative_p0_map", "relative_p0_map"),
        _part2_page("part2_transfer_map", "transfer_map"),
        _part2_page("part2_statistical_uncertainty_map", "statistical_uncertainty_map"),
        _part2_page("part2_model_uncertainty_map", "model_uncertainty_map"),
        _part2_page("part2_fallback_uncertainty_map", "fallback_uncertainty_map"),
        _part2_page("part2_uncertainty_map", "total_uncertainty_map"),
        _part2_page("part2_fit_status_map", "fit_status_map"),
        _part2_page("part2_solution_source_map", "solution_source_map"),
        _part2_page("part2_part3_eligibility_map", "part3_eligibility_map"),
    ]
    application = payload.get("application") or {}
    control = (renderer_inputs or {}).get("control") or {}
    if control.get("status") == "available" and application.get("status") == "available":
        pages.append(_part2_page("part2_transfer_normalization_summary", "normalization_summary"))
    else:
        pages.append(_part2_page("part2_transfer_normalization_summary_unavailable", "closure_status", status="unavailable", unavailable_reason=(control.get("reason") or "verified_control_population_unavailable")))
    t_count = max(0, len(payload.get("t_edges") or ()) - 1)
    response_pages = extract_pion_hgcer_response_fit_pages(payload)
    response_by_t = {}
    for response_page in response_pages:
        response_by_t.setdefault(response_page["t_index"], []).append(response_page)
    for t_index in range(t_count):
        number = t_index + 1
        for field, page_id, spec_key in (
            ("transfer", "part2_transfer_t{}".format(number), "transfer_series"),
            ("P0", "part2_p0_t{}".format(number), "p0_series"),
        ):
            if _page_has_defined_series(payload, t_index, field):
                pages.append(_part2_page(page_id, spec_key, t_index=t_index, field=field))
            else:
                pages.append(_part2_page(page_id + "_unavailable", "closure_status", t_index=t_index, status="unavailable", unavailable_reason="no_resolved_{}_values".format(field)))
        for response_page in response_by_t.get(t_index, ()):
            if response_page["fit_source"] == "direct":
                page_id = "part2_response_fit_t{}_delta{}".format(number, int(response_page["delta_index"]) + 1)
            else:
                page_id = "part2_response_fit_t{}_pooled".format(number)
            required = list(PART2_PAGE_SPECS["response_fit"]["required_roles"])
            if ((response_page["fit"].get("response_display") or {}).get("alternate_response")) is not None:
                required.append("alternate_response")
            else:
                required.append("alternate_response_unavailable")
            pages.append(_part2_page(page_id, "response_fit", t_index=t_index, fit_source=response_page["fit_source"], delta_index=response_page["delta_index"], fit=response_page["fit"], solution=response_page.get("solution") or {}, eligibility=response_page.get("eligibility") or {}, required_roles=required))
        control_entry = (control.get("by_t") or {}).get(t_index) or (control.get("by_t") or {}).get(str(t_index)) or {}
        if control_entry.get("status") == "available" and _application_histograms_for_t(application, t_index):
            pages.append(_part2_page("part2_control_to_transferred_t{}".format(number), "control_transfer", t_index=t_index))
        else:
            pages.append(_part2_page("part2_control_to_transferred_t{}_unavailable".format(number), "closure_status", t_index=t_index, status="unavailable", unavailable_reason=(control_entry.get("reason") or "verified_control_population_unavailable")))
        if _application_histograms_for_t(application, t_index):
            pages.append(_part2_page("part2_proposed_mm_t{}".format(number), "proposed_mm", t_index=t_index, application_scope="per_t"))
        else:
            pages.append(_part2_page("part2_proposed_mm_t{}_unavailable".format(number), "closure_status", t_index=t_index, status="unavailable", unavailable_reason="detached_noRF_host_MM_unavailable"))
        for closure_kind, label, spec_key in (
            ("simc", "simc", "simc_closure"), ("signal", "signal", "signal_closure"),
        ):
            page_id = "part2_{}_closure_t{}".format(label, number)
            if _closure_available(renderer_inputs, t_index, closure_kind) and _application_histograms_for_t(application, t_index):
                required = list(PART2_PAGE_SPECS[spec_key]["required_roles"])
                if closure_kind == "signal":
                    signal = _closure_group(renderer_inputs, t_index, "signal")
                    required.append("KSigma0_reference" if _closure_hist(signal.get("K-Sigma0")) is not None else "KSigma0_unavailable_status")
                pages.append(_part2_page(page_id, spec_key, t_index=t_index, closure_kind=closure_kind, required_roles=required))
            else:
                pages.append(_part2_page(page_id + "_unavailable", "closure_status", t_index=t_index, status="unavailable", unavailable_reason="{}_closure_inputs_unavailable".format(label)))
    if (application.get("global_histograms") or {}):
        global_control = (control.get("global") or {})
        if global_control.get("status") == "available":
            pages.append(_part2_page("part2_control_to_transferred_global", "control_transfer", application_scope="global"))
        else:
            pages.append(_part2_page("part2_control_to_transferred_global_unavailable", "closure_status", status="unavailable", unavailable_reason=(global_control.get("reason") or "verified_global_control_population_unavailable")))
        pages.append(_part2_page("part2_proposed_mm_global", "proposed_mm", application_scope="global", required_roles=list(PART2_PAGE_SPECS["proposed_mm"]["required_roles"]) + ["strict_sum_closure"]))
        for closure_kind, label, spec_key in (
            ("simc", "simc", "simc_closure"), ("signal", "signal", "signal_closure"),
        ):
            page_id = "part2_{}_closure_global".format(label)
            if _closure_available(renderer_inputs, None, closure_kind):
                required = list(PART2_PAGE_SPECS[spec_key]["required_roles"])
                if closure_kind == "signal":
                    signal = _closure_group(renderer_inputs, None, "signal")
                    required.append("KSigma0_reference" if _closure_hist(signal.get("K-Sigma0")) is not None else "KSigma0_unavailable_status")
                pages.append(_part2_page(page_id, spec_key, closure_kind=closure_kind, application_scope="global", required_roles=required))
            else:
                pages.append(_part2_page(page_id + "_unavailable", "closure_status", status="unavailable", unavailable_reason="{}_closure_inputs_unavailable".format(label)))
    else:
        pages.append(_part2_page("part2_proposed_mm_global_unavailable", "closure_status", status="unavailable", unavailable_reason="detached_global_noRF_host_MM_unavailable"))
    pages.append(_part2_page("part2_final_summary", "final_summary"))
    return pages


def _part2_render_name(stem):
    global _PART2_RENDER_SERIAL
    _PART2_RENDER_SERIAL += 1
    return "part2_{}_{}".format(stem, _PART2_RENDER_SERIAL)


def _render_text_page(canvas, page, payload, title_prefix):
    text = ROOT.TPaveText(0.08, 0.12, 0.92, 0.88, "NDC")
    text.SetFillStyle(0); text.SetBorderSize(0); text.SetTextAlign(12)
    if str(title_prefix).strip():
        text.AddText(str(title_prefix).strip())
    unavailable = page["status"] == "unavailable"
    text.AddText(pion_hgcer_transfer_display_text("unavailable" if unavailable else "part2"))
    text.AddText(str(page["page_id"]).replace("_", " "))
    if page["spec_key"] == "mask_audit":
        audit = (payload or {}).get("audit") or {}
        text.AddText("noRF audit: {}".format(audit.get("status", "unavailable")))
        text.AddText("Four-mask contract fingerprint: {}".format(str(((payload or {}).get("pid_contract") or {}).get("fingerprint", "unavailable"))[:16]))
    elif page["spec_key"] == "response_family_audit":
        families = sorted({str((cell.get("fit") or {}).get("model_family")) for cell in (payload or {}).get("cells") or () if (cell.get("fit") or {}).get("model_family")})
        text.AddText("Frozen response families: {}".format(", ".join(families) if families else "unavailable"))
        text.AddText("Positive-only Gamma/Gaussian leakage models are rejected.")
    elif unavailable:
        text.AddText("Reason: {}".format(page.get("unavailable_reason") or "unknown"))
    else:
        text.AddText("Frozen map fingerprint: {}".format(str((payload or {}).get("map_fingerprint") or "unavailable")[:16]))
        text.AddText("Part 3 review status: {}".format((payload or {}).get("Part3_eligibility")))
    text.Draw()
    return [("audit_text" if page["spec_key"] in {"mask_audit", "response_family_audit"} else "status_text", text)]


def _make_tdelta_hist(page, payload, matrix, title):
    delta_edges, t_edges = matrix["delta_edges"], matrix["t_edges"]
    if len(delta_edges) < 2 or len(t_edges) < 2:
        raise RuntimeError("Part-2 graphical map has no frozen t-delta geometry")
    hist = ROOT.TH2D(
        _part2_render_name(page["page_id"]), title,
        len(delta_edges) - 1, array("d", [float(value) for value in delta_edges]),
        len(t_edges) - 1, array("d", [float(value) for value in t_edges]),
    )
    hist.SetDirectory(0); hist.GetXaxis().SetTitle("SHMS delta (%)"); hist.GetYaxis().SetTitle("|t| (GeV^2)")
    return hist


def _render_numeric_map(canvas, page, payload, _renderer_inputs):
    field = page["spec_key"].replace("_map", "")
    field = {"p0": "P0", "relative_p0": "relative_P0_uncertainty", "statistical_uncertainty": "statistical_uncertainty", "model_uncertainty": "model_uncertainty", "fallback_uncertainty": "fallback_uncertainty", "total_uncertainty": "total_uncertainty"}.get(field, field)
    matrix = extract_pion_hgcer_transfer_tdelta_matrix(payload, field)
    role = page["required_roles"][0]
    hist = _make_tdelta_hist(page, payload, matrix, pion_hgcer_transfer_display_text({"P0": "p0", "relative_P0_uncertainty": "relative_p0", "transfer": "transfer", "statistical_uncertainty": "statistical_uncertainty", "model_uncertainty": "model_uncertainty", "fallback_uncertainty": "fallback_uncertainty", "total_uncertainty": "uncertainty"}[field]))
    masked = []
    for t_index, row in enumerate(matrix["values"]):
        for delta_index, value in enumerate(row):
            if value is not None:
                hist.SetBinContent(delta_index + 1, t_index + 1, float(value))
            else:
                masked.append((t_index, delta_index))
    hist.Draw("COLZ TEXT")
    delta_edges, t_edges = matrix["delta_edges"], matrix["t_edges"]
    mask_boxes = []
    for t_index, delta_index in masked:
        box = ROOT.TBox(
            float(delta_edges[delta_index]), float(t_edges[t_index]),
            float(delta_edges[delta_index + 1]), float(t_edges[t_index + 1]),
        )
        box.SetFillColor(ROOT.kGray + 1); box.SetFillStyle(3004); box.SetLineColor(ROOT.kGray + 2)
        box.Draw("same"); mask_boxes.append(box)
    note = ROOT.TLatex(); note.SetNDC(True); note.SetTextSize(0.025)
    note.DrawLatex(0.12, 0.02, "Hatched cells are undefined; they are not zero-filled.")
    return [(role, hist), ("undefined_cell_note", note)] + [("undefined_cell_mask", box) for box in mask_boxes]


def _render_categorical_map(canvas, page, payload, _renderer_inputs):
    field = page["spec_key"].replace("_map", "")
    matrix = extract_pion_hgcer_transfer_categorical_map(payload, field)
    categories = matrix["categories"]
    codes = {category: index + 1 for index, category in enumerate(categories)}
    role = page["required_roles"][0]
    hist = _make_tdelta_hist(page, payload, matrix, pion_hgcer_transfer_display_text({"response_family": "family_audit", "fit_status": "fit_status", "solution_source": "solution_source", "part3_eligibility": "eligibility"}[field]))
    masked = []
    for t_index, row in enumerate(matrix["values"]):
        for delta_index, value in enumerate(row):
            if value is not None:
                hist.SetBinContent(delta_index + 1, t_index + 1, float(codes[value]))
            else:
                masked.append((t_index, delta_index))
    hist.Draw("COLZ TEXT")
    mask_boxes = []
    for t_index, delta_index in masked:
        box = ROOT.TBox(
            float(matrix["delta_edges"][delta_index]), float(matrix["t_edges"][t_index]),
            float(matrix["delta_edges"][delta_index + 1]), float(matrix["t_edges"][t_index + 1]),
        )
        box.SetFillColor(ROOT.kGray + 1); box.SetFillStyle(3004); box.SetLineColor(ROOT.kGray + 2)
        box.Draw("same"); mask_boxes.append(box)
    legend = ROOT.TPaveText(0.70, 0.62, 0.97, 0.92, "NDC")
    legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextSize(0.018)
    legend.AddText("Category codes")
    for category in categories:
        legend.AddText("{} = {}".format(codes[category], category))
    legend.Draw()
    return [(role, hist), ("category_legend", legend)] + [("undefined_cell_mask", box) for box in mask_boxes]


def _runtime_t_label(payload, t_index):
    edges = list((payload or {}).get("t_edges") or ())
    if 0 <= int(t_index) < len(edges) - 1:
        return "|t| = [{:.3f}, {:.3f}] GeV^2".format(float(edges[t_index]), float(edges[t_index + 1]))
    return "|t| range unavailable"


def _runtime_delta_label(payload, delta_index):
    edges = list((payload or {}).get("delta_edges") or ())
    if delta_index is not None and 0 <= int(delta_index) < len(edges) - 1:
        return "delta = [{:.3g}, {:.3g}] %".format(float(edges[delta_index]), float(edges[delta_index + 1]))
    return "delta = same-t pooled"


def _graph_errors(name, rows, *, y_key="value", error_key=None, marker=20, color=None):
    rows = [row for row in rows if _finite(row.get(y_key)) is not None]
    if not rows:
        return None
    x = [float(row["delta"]) for row in rows]
    y = [float(row[y_key]) for row in rows]
    errors = [float(row.get(error_key) or 0.0) for row in rows] if error_key else [0.0] * len(rows)
    graph = ROOT.TGraphErrors(len(rows), array("d", x), array("d", y), array("d", [0.0] * len(rows)), array("d", errors))
    graph.SetName(_part2_render_name(name)); graph.SetMarkerStyle(marker); graph.SetMarkerSize(1.05); graph.SetLineWidth(2)
    if color is not None:
        graph.SetMarkerColor(color); graph.SetLineColor(color)
    return graph


def _render_series(canvas, page, payload, _renderer_inputs):
    field = page["field"]
    detail = extract_pion_hgcer_transfer_series_detail(payload, page["t_index"], field)
    if field == "P0":
        graph = _graph_errors("p0_points", detail["points"], error_key="statistical", marker=20, color=ROOT.kBlue + 1)
        if graph is None:
            raise RuntimeError("Part-2 P0 page selected without direct fitted P0 values")
        graph.SetTitle("Part 2 P0 versus delta - {};SHMS delta (%);P0".format(_runtime_t_label(payload, page["t_index"])))
        graph.Draw("APZ")
        primitives = [("P0_points", graph), ("P0_statistical_errors", graph)]
        for pooled_index, pooled in enumerate(detail.get("pooled") or ()):
            line = ROOT.TLine(float(detail["delta_centers"][0] - detail["delta_half_widths"][0]), float(pooled["P0"]), float(detail["delta_centers"][-1] + detail["delta_half_widths"][-1]), float(pooled["P0"]))
            line.SetLineStyle(2); line.SetLineColor(ROOT.kMagenta + 2); line.Draw("SAME")
            primitives.append(("same_t_pooled_P0_reference_{}".format(pooled_index), line))
        return primitives
    rows = [row for row in detail["points"] if row.get("value") is not None]
    if not rows:
        raise RuntimeError("Part-2 transfer page selected without resolved values")
    total = _graph_errors("transfer_total", [row for row in rows if row.get("total") is not None], error_key="total", marker=1, color=ROOT.kAzure + 2)
    if total is None:
        raise RuntimeError("Part-2 transfer page has no total uncertainty payload")
    total.SetFillColorAlpha(ROOT.kAzure + 1, 0.25); total.SetFillStyle(1001)
    total.SetTitle("Part 2 transfer versus delta - {};SHMS delta (%);R_pi_to_K".format(_runtime_t_label(payload, page["t_index"])))
    total.Draw("A3")
    source_styles = {"direct": (20, ROOT.kBlue + 1), "same_t_delta_bracket": (21, ROOT.kOrange + 7), "same_t_delta_edge": (22, ROOT.kGreen + 2), "same_t_pooled": (33, ROOT.kMagenta + 2)}
    legend = ROOT.TLegend(0.55, 0.67, 0.89, 0.89); legend.SetBorderSize(0)
    legend.AddEntry(total, "total uncertainty", "f")
    primitives = [("total_uncertainty", total)]
    statistical_drawn = False
    for source, (marker, color) in source_styles.items():
        source_rows = [row for row in rows if row.get("source") == source]
        graph = _graph_errors("transfer_{}".format(source), source_rows, error_key="statistical", marker=marker, color=color)
        if graph is None:
            continue
        graph.Draw("PZ SAME")
        legend.AddEntry(graph, source.replace("_", " "), "lep")
        primitives.append(("transfer_points", graph))
        if any(row.get("statistical") is not None for row in source_rows):
            primitives.append(("statistical_errors", graph)); statistical_drawn = True
    if not statistical_drawn:
        note = ROOT.TLatex(); note.SetNDC(True); note.SetTextSize(0.025)
        note.DrawLatex(0.13, 0.15, "Statistical uncertainty undefined for frozen fallback values.")
        primitives.append(("statistical_errors", note))
    legend.Draw(); primitives.append(("solution_source_legend", legend))
    return primitives


def _graph_from_series(name, series, color, style=1):
    if not isinstance(series, dict) or not series.get("x") or not series.get("y"):
        return None
    graph = ROOT.TGraph(len(series["x"]), array("d", [float(value) for value in series["x"]]), array("d", [float(value) for value in series["y"]]))
    graph.SetName(_part2_render_name(name)); graph.SetLineColor(color); graph.SetLineStyle(style); graph.SetLineWidth(2)
    return graph


def _render_response_fit(canvas, page, payload, _renderer_inputs):
    fit = page["fit"]
    solution, eligibility = page.get("solution") or {}, page.get("eligibility") or {}
    display = fit.get("response_display") or {}
    edges, counts = display.get("bin_edges") or (), display.get("bin_counts") or ()
    if len(edges) != len(counts) + 1 or not display.get("primary_response"):
        raise RuntimeError("Part-2 response-fit page has no binned existing-fit display payload")
    canvas.Divide(2, 1)
    canvas.cd(1)
    hist = ROOT.TH1D(_part2_render_name(page["page_id"]), "Part 2 pion HGCer response - {} - {};P_hgcer_npeSum;unweighted prompt-pion records".format(_runtime_t_label(payload, page["t_index"]), _runtime_delta_label(payload, page.get("delta_index"))), len(counts), array("d", [float(value) for value in edges]))
    hist.SetDirectory(0); hist.GetXaxis().SetTitle("P_hgcer_npeSum"); hist.GetYaxis().SetTitle("unweighted prompt-pion records")
    for index, value in enumerate(counts, start=1):
        hist.SetBinContent(index, float(value))
    hist.SetMarkerStyle(20); hist.Draw("E1")
    primary = _graph_from_series("primary_response", display.get("primary_response"), ROOT.kBlue + 1)
    if primary is None:
        raise RuntimeError("Part-2 response-fit primary curve unavailable")
    primary.Draw("L SAME")
    primitives = [("prompt_pion_data", hist), ("primary_response", primary)]
    alternate = _graph_from_series("alternate_response", display.get("alternate_response"), ROOT.kRed + 1, 2)
    if alternate is not None:
        alternate.Draw("L SAME")
        primitives.append(("alternate_response", alternate))
    threshold = float(((payload.get("pid_contract") or {}).get("masks") or {}).get("physical_pion_control", {}).get("value", 2.0))
    line = ROOT.TLine(threshold, 0.0, threshold, max(float(hist.GetMaximum()), 1.0))
    line.SetLineStyle(7); line.SetLineColor(ROOT.kGray + 2); line.Draw("SAME")
    primitives.append(("physical_control_threshold", line))
    zero_line = ROOT.TLine(0.0, 0.0, 0.0, max(float(hist.GetMaximum()), 1.0))
    zero_line.SetLineStyle(3); zero_line.SetLineColor(ROOT.kMagenta + 2); zero_line.Draw("SAME")
    zero_text = ROOT.TLatex(); zero_text.SetTextSize(0.026)
    zero_text.DrawLatex(0.02, max(float(hist.GetMaximum()) * 0.86, 0.8), "NPE = 0 inferred atom; P0 = e^{-mu}")
    primitives.append(("zero_pe_annotation", zero_text)); primitives.append(("zero_pe_atom", zero_line))
    legend = ROOT.TLegend(0.55, 0.68, 0.89, 0.88); legend.SetBorderSize(0)
    legend.AddEntry(hist, "prompt pion data", "lep"); legend.AddEntry(primary, "primary response", "l")
    if alternate is not None:
        legend.AddEntry(alternate, "alternate response", "l")
    else:
        legend.AddEntry(0, "alternate response unavailable", "")
        primitives.append(("alternate_response_unavailable", legend))
    legend.Draw()
    canvas.cd(2)
    profile = fit.get("profile") or {}
    profile_x = [value for value in profile.get("grid_log_mu") or () if value is not None]
    profile_y = list(profile.get("grid_nll") or ())
    profile_pairs = [(x, y) for x, y in zip(profile_x, profile_y) if y is not None]
    if not profile_pairs:
        raise RuntimeError("Part-2 response-fit profile payload unavailable")
    px, py = zip(*profile_pairs)
    profile_graph = ROOT.TGraph(len(px), array("d", px), array("d", py))
    profile_graph.SetName(_part2_render_name("profile")); profile_graph.SetLineWidth(2); profile_graph.SetTitle("mu profile;log(mu);NLL")
    profile_graph.Draw("AL")
    primitives.append(("profile_likelihood", profile_graph))
    summary = ROOT.TPaveText(0.48, 0.53, 0.90, 0.89, "NDC")
    summary.SetFillStyle(0); summary.SetBorderSize(1); summary.SetTextSize(0.022); summary.SetTextAlign(12)
    summary.AddText("{}; {}".format(_runtime_t_label(payload, page["t_index"]), _runtime_delta_label(payload, page.get("delta_index"))))
    summary.AddText("source={}  family={}".format(page.get("fit_source"), fit.get("model_family")))
    summary.AddText("alternate={} ({})".format(fit.get("alternate_model_family", "unavailable"), fit.get("alternate_fit_status", "unavailable")))
    summary.AddText("prompt records={}  mu={:.5g}".format(fit.get("fit_record_count", "n/a"), float(fit.get("mu") or 0.0)))
    summary.AddText("P0={:.5g} +/- {:.3g}  rel={}".format(float(fit.get("P0") or 0.0), float(fit.get("P0_statistical_uncertainty") or 0.0), fit.get("P0_relative_statistical_uncertainty")))
    summary.AddText("R={}; stat={}; model={}; fallback={}; total={}".format(solution.get("transfer"), solution.get("transfer_statistical_uncertainty"), solution.get("transfer_model_uncertainty"), solution.get("transfer_fallback_uncertainty"), solution.get("transfer_total_uncertainty")))
    summary.AddText("fit={} profile={} toys={:.3g}".format(fit.get("fit_status"), (fit.get("profile") or {}).get("status"), float(fit.get("accepted_toy_fraction") or 0.0)))
    summary.AddText("covariance cond={} eligibility={}".format(fit.get("covariance_condition_number"), eligibility.get("status", "unavailable")))
    summary.Draw(); primitives.append(("fit_metadata", summary))
    return primitives


def _display_clone(hist, name, color, *, style=1):
    clone = hist.Clone(_part2_render_name(name))
    if hasattr(clone, "SetDirectory"):
        clone.SetDirectory(0)
    clone.SetLineColor(color); clone.SetLineStyle(style); clone.SetLineWidth(2); clone.SetFillStyle(0)
    return clone


def _render_hist_overlay(canvas, title, sources):
    clones = [(role, _display_clone(hist, role, color, style=style), label) for role, hist, label, color, style in sources if hist is not None]
    if not clones:
        raise RuntimeError("Part-2 graphical overlay has no existing detached histogram")
    minimum = min(float(clone.GetMinimum()) for _, clone, _ in clones)
    maximum = max(float(clone.GetMaximum()) for _, clone, _ in clones)
    first_role, first, _ = clones[0]
    first.SetTitle(title); first.GetXaxis().SetTitle("MM (GeV)"); first.GetYaxis().SetTitle("weighted counts")
    first.SetMinimum(minimum * 1.10 if minimum < 0.0 else 0.0); first.SetMaximum(maximum * 1.18 if maximum > 0.0 else 1.0)
    first.Draw("HIST")
    for _, clone, _ in clones[1:]:
        clone.Draw("HIST SAME")
    legend = ROOT.TLegend(0.58, 0.68, 0.89, 0.89); legend.SetBorderSize(0)
    for _, clone, label in clones:
        legend.AddEntry(clone, label, "l")
    legend.Draw()
    return [(role, clone) for role, clone, _ in clones] + [("overlay_legend", legend)]


def _render_proposed_mm(canvas, page, payload, _renderer_inputs):
    application = payload.get("application") or {}
    if page.get("application_scope") == "global":
        histograms = application.get("global_histograms") or {}
        title = "Detached global proposed pion subtraction from noRF host"
        metrics = dict(application.get("strict_global_sums") or {})
        metrics.update({
            "host_integral": _finite(metrics.get("host")), "pion_integral": _finite(metrics.get("pion")),
            "clean_integral": _finite(metrics.get("clean")),
            "pion_to_host": (_finite(metrics.get("pion")) / _finite(metrics.get("host")) if _finite(metrics.get("host")) not in (None, 0.0) and _finite(metrics.get("pion")) is not None else None),
        })
    else:
        histograms = _application_histograms_for_t(application, page["t_index"])
        title = "Part 2 proposed pion MM - {}".format(_runtime_t_label(payload, page["t_index"]))
        metrics = next((entry for entry in (application.get("t_products") or ()) if int(entry.get("t_index", -1)) == int(page["t_index"])), {})
    canvas.Divide(1, 2)
    canvas.cd(1)
    primitives = _render_hist_overlay(canvas, title, [
        ("host_mm", histograms.get("host"), "proton-cleaned noRF host", ROOT.kBlack, 1),
        ("pion_mm", histograms.get("pion"), "proposed pion MM", ROOT.kRed + 1, 1),
        ("clean_mm", histograms.get("clean"), "proposed clean MM", ROOT.kBlue + 1, 1),
    ])
    summary = ROOT.TPaveText(0.55, 0.55, 0.90, 0.89, "NDC")
    summary.SetFillStyle(0); summary.SetBorderSize(1); summary.SetTextSize(0.020); summary.SetTextAlign(12)
    summary.AddText("host={}; pion={}; clean={}".format(metrics.get("host_integral", metrics.get("host")), metrics.get("pion_integral", metrics.get("pion")), metrics.get("clean_integral", metrics.get("clean"))))
    summary.AddText("pion/host={}".format(metrics.get("pion_to_host")))
    summary.AddText("pion > host bins={}".format(metrics.get("pion_exceeds_host_bin_count", "n/a")))
    summary.AddText("largest |pion/host|={}".format(metrics.get("maximum_local_contamination", "n/a")))
    summary.AddText("most negative clean={}".format(metrics.get("most_negative_clean_bin", "n/a")))
    summary.AddText("unresolved transferred records={}".format(metrics.get("unresolved_records", "n/a")))
    summary.AddText("strict closure={}".format(metrics.get("closure", "n/a")))
    summary.Draw(); primitives.append(("summary_text", summary))
    if page.get("application_scope") == "global":
        primitives.append(("strict_sum_closure", summary))
    canvas.cd(2)
    pion_zoom = _display_clone(histograms.get("pion"), "pion_zoom", ROOT.kRed + 1) if histograms.get("pion") is not None else None
    if pion_zoom is None:
        raise RuntimeError("Part-2 proposed MM pion zoom histogram missing")
    pion_zoom.SetTitle("same absolute normalization - zoomed vertical scale;MM (GeV);weighted counts")
    pion_zoom.SetMinimum(min(float(pion_zoom.GetMinimum()), 0.0) * 1.10 if float(pion_zoom.GetMinimum()) < 0.0 else 0.0)
    pion_zoom.SetMaximum(max(float(pion_zoom.GetMaximum()) * 1.18, 1.0e-12)); pion_zoom.Draw("HIST")
    primitives.append(("pion_zoom", pion_zoom))
    return primitives


def _render_control_transfer(canvas, page, payload, renderer_inputs):
    application = payload.get("application") or {}
    control = (renderer_inputs or {}).get("control") or {}
    if page.get("application_scope") == "global":
        section, histograms = control.get("global") or {}, application.get("global_histograms") or {}
        title_scope = "strict sum of canonical t bins"
    else:
        section = (control.get("by_t") or {}).get(page["t_index"]) or (control.get("by_t") or {}).get(str(page["t_index"])) or {}
        histograms = _application_histograms_for_t(application, page["t_index"])
        title_scope = _runtime_t_label(payload, page["t_index"])
    before, after = section.get("before"), histograms.get("pion")
    primitives = _render_hist_overlay(canvas, "Part 2 signed pion control to transferred pion - {}".format(title_scope), [
        ("pion_control_mm", before, "signed pion control before R", ROOT.kBlack, 1),
        ("transferred_pion_mm", after, "HGCer transferred pion", ROOT.kRed + 1, 1),
    ])
    before_integral, after_integral = float(before.Integral()), float(after.Integral())
    summary = ROOT.TPaveText(0.55, 0.68, 0.90, 0.88, "NDC")
    summary.SetFillStyle(0); summary.SetBorderSize(1); summary.SetTextSize(0.024); summary.SetTextAlign(12)
    summary.AddText("control integral={:.6g}".format(before_integral))
    summary.AddText("transferred integral={:.6g}".format(after_integral))
    summary.AddText("transferred/control={}".format(after_integral / before_integral if before_integral != 0.0 else "undefined"))
    summary.Draw(); primitives.append(("suppression_summary", summary))
    return primitives


def _normalization_rows(payload, renderer_inputs):
    application, control = (payload or {}).get("application") or {}, (renderer_inputs or {}).get("control") or {}
    rows = []
    for entry in application.get("t_products") or ():
        t_index = int(entry.get("t_index", -1))
        section = (control.get("by_t") or {}).get(t_index) or (control.get("by_t") or {}).get(str(t_index)) or {}
        before = _finite(section.get("before_integral"))
        after = _finite(entry.get("pion_integral"))
        host = _finite(entry.get("host_integral"))
        rows.append({"t_index": t_index, "before": before, "after": after, "effective_transfer": after / before if before not in (None, 0.0) and after is not None else None, "host": host, "pion_to_host": after / host if host not in (None, 0.0) and after is not None else None})
    return rows


def _render_normalization_summary(canvas, page, payload, renderer_inputs):
    rows = _normalization_rows(payload, renderer_inputs)
    if not rows:
        raise RuntimeError("Part-2 transfer normalization summary has no verified rows")
    x = [float(row["t_index"] + 1) for row in rows]
    before = ROOT.TGraph(len(rows), array("d", x), array("d", [float(row["before"] or 0.0) for row in rows]))
    after = ROOT.TGraph(len(rows), array("d", x), array("d", [float(row["after"] or 0.0) for row in rows]))
    before.SetName(_part2_render_name("control_integrals")); before.SetMarkerStyle(20); before.SetLineColor(ROOT.kBlack); before.SetTitle("Part 2 transfer normalization;canonical t index;signed MM integral")
    after.SetName(_part2_render_name("transferred_integrals")); after.SetMarkerStyle(21); after.SetLineColor(ROOT.kRed + 1)
    before.Draw("ALP"); after.Draw("LP SAME")
    legend = ROOT.TLegend(0.55, 0.72, 0.89, 0.88); legend.SetBorderSize(0)
    legend.AddEntry(before, "signed pion control before R", "lp"); legend.AddEntry(after, "transferred pion after R", "lp"); legend.Draw()
    summary = ROOT.TPaveText(0.12, 0.58, 0.50, 0.89, "NDC")
    summary.SetFillStyle(0); summary.SetBorderSize(1); summary.SetTextSize(0.021); summary.SetTextAlign(12)
    for row in rows:
        summary.AddText("t{}: before={} after={} R_eff={} pion/host={}".format(row["t_index"] + 1, row["before"], row["after"], row["effective_transfer"], row["pion_to_host"]))
    summary.Draw()
    return [("pion_control_integral_series", before), ("transferred_pion_integral_series", after), ("normalization_summary", summary), ("normalization_legend", legend)]


def _render_final_summary(canvas, _page, payload, renderer_inputs):
    application = (payload or {}).get("application") or {}
    rows = _normalization_rows(payload, renderer_inputs)
    text = ROOT.TPaveText(0.08, 0.12, 0.92, 0.88, "NDC")
    text.SetFillStyle(0); text.SetBorderSize(1); text.SetTextAlign(12); text.SetTextSize(0.026)
    text.AddText("PART 2 FINAL DIAGNOSTIC SUMMARY - NON-AUTHORITATIVE")
    text.AddText("Frozen map fingerprint: {}".format(str((payload or {}).get("map_fingerprint") or "unavailable")[:16]))
    text.AddText("strict host/pion/clean/closure: {} / {} / {} / {}".format((application.get("strict_global_sums") or {}).get("host"), (application.get("strict_global_sums") or {}).get("pion"), (application.get("strict_global_sums") or {}).get("clean"), (application.get("strict_global_sums") or {}).get("closure")))
    control_global = ((renderer_inputs or {}).get("control") or {}).get("global") or {}
    text.AddText("strict global pion control before R: {} ({})".format(control_global.get("before_integral"), control_global.get("construction", "unavailable")))
    for row in rows:
        text.AddText("t{}: control={} transferred={} ratio={}".format(row["t_index"] + 1, row["before"], row["after"], row["effective_transfer"]))
    text.AddText("Part 3 remains review-only; production subtraction unchanged.")
    text.Draw()
    return [("part2_summary_table", text)]


def _render_closure(canvas, page, payload, renderer_inputs):
    application = payload.get("application") or {}
    if page.get("application_scope") == "global":
        histograms = application.get("global_histograms") or {}
        closure_index = None
    else:
        histograms = _application_histograms_for_t(application, page["t_index"])
        closure_index = page["t_index"]
    group = _closure_group(renderer_inputs, closure_index, page["closure_kind"])
    if page["closure_kind"] == "simc":
        reference_role, proposed_role = "simc_reference", "proposed_pion_mm"
        title, proposed = pion_hgcer_transfer_display_text("simc_closure"), histograms.get("pion")
        references = [(key, _closure_hist(value)) for key, value in group.items() if _closure_hist(value) is not None]
        sources = [(proposed_role, proposed, "proposed pion MM", ROOT.kRed + 1, 1)]
        sources.extend((reference_role if index == 0 else "simc_component_{}".format(index), hist, str(key), ROOT.kAzure + index + 1, 2) for index, (key, hist) in enumerate(references))
        primitives = _render_hist_overlay(canvas, title, sources)
        annotation = ROOT.TPaveText(0.12, 0.76, 0.50, 0.89, "NDC")
        annotation.SetFillStyle(0); annotation.SetBorderSize(1); annotation.SetTextSize(0.022)
        annotation.AddText("Existing fitted SIMC components")
        annotation.AddText("Diagnostic closure only")
        annotation.AddText("Not used to determine R_pi_to_K")
        annotation.Draw(); primitives.append(("simc_closure_annotation", annotation))
        return primitives
    else:
        proposed_role = "proposed_clean_mm"
        title, proposed = pion_hgcer_transfer_display_text("signal_closure"), histograms.get("clean")
        lambda_entry, sigma_entry = group.get("K-Lambda") or {}, group.get("K-Sigma0")
        lambda_hist = _closure_hist(lambda_entry)
        if lambda_hist is None:
            raise RuntimeError("Part-2 signal closure K-Lambda display reference unavailable")
        sigma_hist = _closure_hist(sigma_entry)
        sources = [
            (proposed_role, proposed, "proposed clean MM", ROOT.kBlue + 1, 1),
            ("KLambda_reference", lambda_hist, "K-Lambda established display reference", ROOT.kGreen + 2, 2),
        ]
        if sigma_hist is not None:
            sources.append(("KSigma0_reference", sigma_hist, "K-Sigma0 reference", ROOT.kCyan + 2, 2))
        primitives = _render_hist_overlay(canvas, title, sources)
        status = ROOT.TPaveText(0.12, 0.73, 0.53, 0.89, "NDC")
        status.SetFillStyle(0); status.SetBorderSize(1); status.SetTextSize(0.021); status.SetTextAlign(12)
        if isinstance(lambda_entry, dict):
            status.AddText("K-Lambda display source={}".format(lambda_entry.get("source", "unavailable")))
            status.AddText("normalization={}; scale={}".format(lambda_entry.get("normalization", "unavailable"), lambda_entry.get("scale", "unavailable")))
        else:
            status.AddText("K-Lambda display normalization metadata unavailable")
        if sigma_hist is None:
            status.AddText("K-Sigma0 unavailable for this scope")
            primitives.append(("KSigma0_unavailable_status", status))
        else:
            status.AddText("K-Sigma0 available for this scope")
        status.Draw(); primitives.append(("normalization_status", status))
        return primitives


def _validate_part2_graphical_page(page, primitives):
    if page["page_kind"] != "graphical":
        return
    roles = {role for role, _ in primitives}
    missing = set(page["required_roles"]) - roles
    if missing:
        raise RuntimeError("Part-2 graphical page '{}' is missing semantic primitives: {}".format(page["page_id"], ", ".join(sorted(missing))))


def _render_part2_page(canvas, page, payload, title_prefix, renderer_inputs):
    if page["page_kind"] == "text":
        return _render_text_page(canvas, page, payload, title_prefix)
    if page["renderer"] == "numeric_map":
        return _render_numeric_map(canvas, page, payload, renderer_inputs)
    if page["renderer"] == "categorical_map":
        return _render_categorical_map(canvas, page, payload, renderer_inputs)
    if page["renderer"] == "series":
        return _render_series(canvas, page, payload, renderer_inputs)
    if page["renderer"] == "response_fit":
        return _render_response_fit(canvas, page, payload, renderer_inputs)
    if page["renderer"] == "control_transfer":
        return _render_control_transfer(canvas, page, payload, renderer_inputs)
    if page["renderer"] == "proposed_mm":
        return _render_proposed_mm(canvas, page, payload, renderer_inputs)
    if page["renderer"] == "normalization_summary":
        return _render_normalization_summary(canvas, page, payload, renderer_inputs)
    if page["renderer"] in {"simc_closure", "signal_closure"}:
        return _render_closure(canvas, page, payload, renderer_inputs)
    if page["renderer"] == "final_summary":
        return _render_final_summary(canvas, page, payload, renderer_inputs)
    raise RuntimeError("unknown Part-2 page renderer '{}'".format(page["renderer"]))


def render_pion_hgcer_zerope_transfer_pages(pdf_name, payload, *, title_prefix="", page_manifest=None, close_pdf=True, renderer_inputs=None):
    """Render frozen Part-2 diagnostics and prove each graphical primitive role."""
    if ROOT is None:
        return []
    manifest = expected_pion_hgcer_transfer_page_manifest(payload, renderer_inputs=renderer_inputs)
    emitted = []
    for index, page in enumerate(manifest):
        canvas = ROOT.TCanvas(_part2_render_name("canvas"), page["page_id"], 1100, 760)
        primitives = _render_part2_page(canvas, page, payload or {}, title_prefix, renderer_inputs or {})
        _validate_part2_graphical_page(page, primitives)
        record = {
            "page_id": page["page_id"], "t_index": page.get("t_index"),
            "page_kind": page["page_kind"],
            "status": "unavailable" if page["status"] == "unavailable" else "rendered",
            "semantic_primitives": [{"role": role, "type": primitive.ClassName()} for role, primitive in primitives],
            "unavailable_reason": page.get("unavailable_reason"),
        }
        if page["spec_key"] == "response_fit":
            record.update({"delta_index": page.get("delta_index"), "fit_source": page.get("fit_source")})
        if page["spec_key"] in {"proposed_mm", "control_transfer"}:
            application = (payload or {}).get("application") or {}
            if page.get("application_scope") == "global":
                summary = application.get("strict_global_sums") or {}
                record.update({"host_integral": summary.get("host"), "pion_integral": summary.get("pion"), "clean_integral": summary.get("clean")})
            else:
                entry = next((item for item in (application.get("t_products") or ()) if int(item.get("t_index", -1)) == int(page.get("t_index", -1))), {})
                record.update({"host_integral": entry.get("host_integral"), "pion_integral": entry.get("pion_integral"), "clean_integral": entry.get("clean_integral")})
        if page["spec_key"] == "control_transfer":
            control = (renderer_inputs or {}).get("control") or {}
            control_entry = (control.get("global") or {}) if page.get("application_scope") == "global" else ((control.get("by_t") or {}).get(page.get("t_index")) or (control.get("by_t") or {}).get(str(page.get("t_index"))) or {})
            before_integral = _finite(control_entry.get("before_integral"))
            after_integral = _finite(record.get("pion_integral"))
            record.update({
                "control_before_integral": before_integral,
                "transferred_pion_integral": after_integral,
                "control_to_transferred_ratio": after_integral / before_integral if before_integral not in (None, 0.0) and after_integral is not None else None,
                "control_construction": control_entry.get("construction") if page.get("application_scope") == "global" else "verified_per_t_control_projection",
            })
        suffix = ")" if close_pdf and index == len(manifest) - 1 else ""
        canvas.Print(str(pdf_name) + suffix)
        emitted.append(record)
    if isinstance(page_manifest, list):
        page_manifest.extend(_json_ready(emitted))
    return _json_ready(emitted)


def render_pion_hgcer_zerope_transfer_failure_page(pdf_name, reason, *, title_prefix="", page_manifest=None, close_pdf=True):
    """Close a diagnostic PDF after a Part-2 presentation-only failure."""
    if ROOT is None:
        return []
    page = _part2_page("part2_rendering_failure", "unavailable", status="unavailable", unavailable_reason=str(reason))
    canvas = ROOT.TCanvas(_part2_render_name("failure"), page["page_id"], 1100, 760)
    _render_text_page(canvas, page, {"reason": str(reason)}, title_prefix)
    canvas.Print(str(pdf_name) + (")" if close_pdf else ""))
    record = {"page_id": page["page_id"], "t_index": None, "page_kind": "text", "status": "unavailable", "semantic_primitives": [], "unavailable_reason": str(reason)}
    if isinstance(page_manifest, list):
        page_manifest.append(record)
    return [record]
