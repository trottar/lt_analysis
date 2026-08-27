"""Detached Part-2 HGCer zero-photoelectron pion transfer diagnostics.

The module is intentionally pure-Python first: response fitting, masks,
fallbacks, and map freezing do not require ROOT.  ROOT is used only for the
optional PDF products and proposed MM histogram clones.  Nothing here mutates
the authoritative pion/proton subtraction objects.
"""

from __future__ import annotations

from array import array
from collections.abc import Mapping
from copy import deepcopy
import hashlib
import json
import math
import re
from types import MappingProxyType

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
from scipy.stats import gamma as gamma_distribution
from scipy.stats import nbinom as negative_binomial_distribution
from scipy.stats import poisson as poisson_distribution

try:  # Pure-Python contracts deliberately run without PyROOT.
    import ROOT
except ImportError:  # pragma: no cover - depends on the local analysis host
    ROOT = None

from canonical_binning import find_canonical_bin


TRANSFER_LABEL = "PION HGCer ZERO-PE TRANSFER - PART 2 - NON-AUTHORITATIVE"
_MASK_FIELD = "P_hgcer_npeSum"
_MASK_PATTERN = re.compile(
    r"\bP_hgcer_npeSum\b\s*(==|>=|<=|>|<)\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
)


class PionHGCerTransferFailure(RuntimeError):
    """Recoverable Part-2 build error with an explicit diagnostic stage."""

    def __init__(self, diagnostic_stage, original_exception, provenance=None):
        self.diagnostic_stage = str(diagnostic_stage)
        self.original_exception = original_exception
        self.provenance = provenance
        super().__init__(str(original_exception))


def _finite(value, default=None):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def _json_ready(value):
    if isinstance(value, MappingProxyType):
        value = dict(value)
    if isinstance(value, dict):
        return {str(key): _json_ready(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(child) for child in value]
    if isinstance(value, np.generic):
        return _json_ready(value.item())
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, np.ndarray):
        return _json_ready(value.tolist())
    return value


def _fingerprint(value):
    encoded = json.dumps(_json_ready(value), sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _freeze(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(child) for key, child in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(child) for child in value)
    if isinstance(value, tuple):
        return tuple(_freeze(child) for child in value)
    return value


def _mask_operator(operator):
    text = str(operator or "").strip()
    if text not in {"==", ">", ">=", "<", "<="}:
        raise ValueError("unsupported HGCer mask operator '{}'".format(operator))
    return text


def normalize_hgcer_mask(mask, *, source=""):
    """Return one canonical scalar HGCer mask, never an inferred threshold."""
    if not isinstance(mask, dict):
        raise ValueError("HGCer mask from {} must be a mapping".format(source or "unknown source"))
    field = str(mask.get("field") or _MASK_FIELD).strip()
    if field != _MASK_FIELD:
        raise ValueError("HGCer mask from {} must constrain {}".format(source or "unknown source", _MASK_FIELD))
    value = _finite(mask.get("value"))
    if value is None:
        raise ValueError("HGCer mask from {} has non-finite threshold".format(source or "unknown source"))
    operator = _mask_operator(mask.get("operator"))
    expression = str(mask.get("expression") or "{} {} {}".format(field, operator, value)).strip()
    return {
        "field": field,
        "operator": operator,
        "value": float(value),
        "expression": expression,
        "source": str(mask.get("source") or source or "unknown"),
    }


def parse_hgcer_mask(expression, *, source=""):
    """Extract exactly one direct NPE predicate from replay-cut provenance."""
    text = str(expression or "").strip()
    matches = list(_MASK_PATTERN.finditer(text))
    if len(matches) != 1:
        raise ValueError(
            "HGCer provenance from {} must contain exactly one direct {} predicate".format(
                source or "unknown source", _MASK_FIELD
            )
        )
    match = matches[0]
    return normalize_hgcer_mask(
        {
            "field": _MASK_FIELD,
            "operator": match.group(1),
            "value": float(match.group(2)),
            "expression": text,
            "source": source,
        },
        source=source,
    )


def hgcer_mask_accepts(mask, value):
    mask = normalize_hgcer_mask(mask)
    scalar = _finite(value)
    if scalar is None:
        return False
    threshold = float(mask["value"])
    operator = mask["operator"]
    if operator == "==":
        return scalar == threshold
    if operator == ">":
        return scalar > threshold
    if operator == ">=":
        return scalar >= threshold
    if operator == "<":
        return scalar < threshold
    return scalar <= threshold


def _same_mask(left, right):
    return (
        left["field"] == right["field"]
        and left["operator"] == right["operator"]
        and abs(float(left["value"]) - float(right["value"])) <= 1.0e-12
    )


def _manifest_tree_entry(manifest, tree_name):
    manifest = manifest or {}
    entries = manifest.get("trees") if isinstance(manifest, dict) else None
    if not isinstance(entries, dict):
        entries = manifest if isinstance(manifest, dict) else {}
    entry = entries.get(str(tree_name))
    if not isinstance(entry, dict):
        raise ValueError("missing persisted PID provenance for tree '{}'".format(tree_name))
    expression = entry.get("expression") or entry.get("cut_expression")
    if not expression:
        raise ValueError("persisted PID provenance has no expression for tree '{}'".format(tree_name))
    return entry, str(expression)


def resolve_hgcer_transfer_masks(source_entries, physical_control_mask):
    """Resolve/fingerprint the four Part-2 masks from persisted source truth.

    ``source_entries`` is ``{"kaon": [{tree_name, manifest, ...}],
    "pion": [...]}``.  Every source of a side must resolve to the same scalar
    NPE predicate; otherwise a map cannot claim a single detector response.
    """
    resolved = {}
    provenance = {}
    for side in ("kaon", "pion"):
        masks = []
        rows = []
        for entry in list((source_entries or {}).get(side) or []):
            tree_name = str((entry or {}).get("tree_name") or "")
            if not tree_name.endswith("_noRF") or str((entry or {}).get("rf_state") or "") != "noRF":
                raise ValueError("Part-2 PID source '{}' is not explicitly noRF".format(tree_name))
            expected_role = "kaon_pid" if side == "kaon" else "pion_pid"
            if str((entry or {}).get("pid_role") or "") != expected_role:
                raise ValueError("Part-2 PID source '{}' has an invalid PID role".format(tree_name))
            manifest_entry, expression = _manifest_tree_entry((entry or {}).get("manifest"), tree_name)
            if str(manifest_entry.get("rf_state") or "") != "noRF":
                raise ValueError("persisted PID provenance for '{}' is not explicitly noRF".format(tree_name))
            if str(manifest_entry.get("pid_role") or "") != expected_role:
                raise ValueError("persisted PID provenance for '{}' has an invalid PID role".format(tree_name))
            source = "{}:{}".format(side, tree_name)
            mask = parse_hgcer_mask(expression, source=source)
            masks.append(mask)
            rows.append(
                {
                    "tree_name": tree_name,
                    "source_label": (entry or {}).get("source_label"),
                    "expression": expression,
                    "cut_name": manifest_entry.get("cut_name"),
                    "cut_file": manifest_entry.get("cut_file"),
                    "manifest_fingerprint": manifest_entry.get("fingerprint"),
                    "rf_state": entry.get("rf_state"),
                    "pid_role": entry.get("pid_role"),
                    "signed_coefficient": _finite(entry.get("signed_coefficient")),
                    "coordinate_fingerprint": entry.get("coordinate_fingerprint"),
                    "proton_factor_scope": entry.get("proton_factor_scope"),
                    "mask": mask,
                }
            )
        if not masks:
            raise ValueError("no persisted {} PID-tree provenance is available".format(side))
        if any(not _same_mask(masks[0], candidate) for candidate in masks[1:]):
            raise ValueError("{} PID-tree HGCer predicates disagree across sources".format(side))
        resolved[side] = masks[0]
        provenance[side] = rows

    kaon_mask = resolved["kaon"]
    if not (kaon_mask["operator"] == "==" and abs(float(kaon_mask["value"])) <= 1.0e-12):
        raise ValueError("Part-2 zero-photoelectron transfer requires kaon PID mask '{} == 0'".format(_MASK_FIELD))
    pion_mask = resolved["pion"]
    if not (pion_mask["operator"] == ">" and abs(float(pion_mask["value"])) <= 1.0e-12):
        raise ValueError("Part-2 pion PID response sample requires '{} > 0'".format(_MASK_FIELD))
    control_mask = normalize_hgcer_mask(physical_control_mask, source="pion_control_cache")
    excludes_zero = (
        (control_mask["operator"] == ">" and float(control_mask["value"]) >= 0.0)
        or (control_mask["operator"] == ">=" and float(control_mask["value"]) > 0.0)
    )
    if not excludes_zero:
        raise ValueError("physical pion-control HGCer mask must select positive NPE response")
    response_sample_mask = deepcopy(pion_mask)
    response_sample_mask["source"] = "prompt_pion_tree_allcuts_unweighted_response_sample"
    masks = {
        "S_K_tree": kaon_mask,
        "S_pi_tree": pion_mask,
        "S_pi_response_sample": response_sample_mask,
        "S_pi_physical_control": control_mask,
    }
    return {
        "masks": masks,
        "response_sample_provenance": {
            "source_label": "prompt",
            "allcuts": True,
            "weighting": "unweighted",
            "tree_mask": deepcopy(pion_mask),
        },
        "source_provenance": provenance,
        "pid_selection_fingerprint": _fingerprint({
            "masks": masks,
            "response_sample_provenance": {
                "source_label": "prompt", "allcuts": True,
                "weighting": "unweighted", "tree_mask": pion_mask,
            },
            "sources": provenance,
        }),
    }


def _is_integer_like(values, tolerance):
    if len(values) == 0:
        return False
    array_values = np.asarray(values, dtype=float)
    return bool(np.all(np.abs(array_values - np.rint(array_values)) <= float(tolerance)))


def choose_pion_hgcer_response_family(values, config):
    requested = str((config or {}).get("response_family") or "auto").strip().lower()
    if requested not in {"auto", "poisson", "compound_poisson"}:
        raise ValueError(
            "response_family '{}' is not a zero-atom photoelectron response".format(requested)
        )
    integer_like = _is_integer_like(values, (config or {}).get("integer_tolerance", 1.0e-6))
    if requested == "poisson":
        if not integer_like:
            raise ValueError("configured Poisson response requires integer-like prompt pion NPE")
        return "zero_truncated_poisson", "zero_truncated_negative_binomial", True
    if requested == "compound_poisson":
        return "zero_truncated_compound_poisson_gamma", "zero_truncated_compound_poisson_exponential", False
    if integer_like:
        return "zero_truncated_poisson", "zero_truncated_negative_binomial", True
    return "zero_truncated_compound_poisson_gamma", "zero_truncated_compound_poisson_exponential", False


def _poisson_p0(mu):
    return float(math.exp(-float(mu)))


def _model_kind_parameters(theta, family):
    parameters = np.exp(np.asarray(theta, dtype=float))
    if not np.all(np.isfinite(parameters)) or np.any(parameters <= 0.0):
        raise ValueError("nonphysical response parameters")
    if family == "zero_truncated_poisson":
        return {"mu": float(parameters[0])}
    if family == "zero_truncated_negative_binomial":
        # Parameterize the NB alternate by its zero-photoelectron rate rather
        # than its ordinary mean.  With m = r(exp(mu/r)-1), a negative
        # binomial with shape r and mean m has P(N=0)=exp(-mu), so the zero
        # atom remains structurally tied to the same mu used by every family.
        mu, dispersion = float(parameters[0]), float(parameters[1])
        exponent = mu / dispersion
        if exponent > 700.0:
            raise ValueError("negative-binomial mean overflow")
        return {
            "mu": mu,
            "dispersion": dispersion,
            "mean": float(dispersion * math.expm1(exponent)),
        }
    if family == "zero_truncated_compound_poisson_gamma":
        return {"mu": float(parameters[0]), "gain_shape": float(parameters[1]), "gain_scale": float(parameters[2])}
    if family == "zero_truncated_compound_poisson_exponential":
        return {"mu": float(parameters[0]), "gain_scale": float(parameters[1])}
    raise ValueError("unknown response family '{}'".format(family))


def _poisson_weights(mu, *, tolerance=1.0e-12, maximum_terms=1000):
    weight = math.exp(-float(mu))
    weights = []
    cumulative = weight
    for count in range(1, int(maximum_terms) + 1):
        weight *= float(mu) / float(count)
        weights.append(weight)
        cumulative += weight
        if count >= int(math.ceil(mu)) and 1.0 - cumulative <= tolerance:
            break
    return weights


def _model_p0(parameters, family):
    return _poisson_p0(parameters["mu"])


def _model_log_probability(values, parameters, family):
    values = np.asarray(values, dtype=float)
    if np.any(values < 0.0):
        return np.full(values.shape, -np.inf, dtype=float)
    if family == "zero_truncated_poisson":
        rounded = np.rint(values).astype(int)
        if not np.all(np.abs(values - rounded) <= 1.0e-6):
            return np.full(values.shape, -np.inf, dtype=float)
        return poisson_distribution.logpmf(rounded, parameters["mu"])
    if family == "zero_truncated_negative_binomial":
        rounded = np.rint(values).astype(int)
        if not np.all(np.abs(values - rounded) <= 1.0e-6):
            return np.full(values.shape, -np.inf, dtype=float)
        dispersion = parameters["dispersion"]
        probability = dispersion / (dispersion + parameters["mean"])
        return negative_binomial_distribution.logpmf(rounded, dispersion, probability)
    mu = float(parameters["mu"])
    gain_shape = 1.0 if family.endswith("exponential") else float(parameters["gain_shape"])
    gain_scale = float(parameters["gain_scale"])
    positive = values > 0.0
    output = np.full(values.shape, -np.inf, dtype=float)
    if not np.any(positive):
        return output
    x_values = values[positive]
    density = np.zeros(x_values.shape, dtype=float)
    for count, weight in enumerate(_poisson_weights(mu), start=1):
        density += weight * gamma_distribution.pdf(
            x_values, a=float(count) * gain_shape, scale=gain_scale
        )
    output[positive] = np.log(np.maximum(density, np.finfo(float).tiny))
    return output


def _model_cdf(value, parameters, family):
    scalar = _finite(value)
    if scalar is None:
        return float("nan")
    if scalar < 0.0:
        return 0.0
    if family == "zero_truncated_poisson":
        return float(poisson_distribution.cdf(math.floor(scalar), parameters["mu"]))
    if family == "zero_truncated_negative_binomial":
        dispersion = parameters["dispersion"]
        probability = dispersion / (dispersion + parameters["mean"])
        return float(negative_binomial_distribution.cdf(math.floor(scalar), dispersion, probability))
    mu = float(parameters["mu"])
    gain_shape = 1.0 if family.endswith("exponential") else float(parameters["gain_shape"])
    gain_scale = float(parameters["gain_scale"])
    total = _poisson_p0(mu)
    for count, weight in enumerate(_poisson_weights(mu), start=1):
        total += weight * float(gamma_distribution.cdf(scalar, a=float(count) * gain_shape, scale=gain_scale))
    return min(max(float(total), 0.0), 1.0)


def _discrete_family(family):
    return family in {"zero_truncated_poisson", "zero_truncated_negative_binomial"}


def model_mask_probability(mask, parameters, family):
    """Integrate a zero-atom response exactly over one scalar HGCer mask."""
    mask = normalize_hgcer_mask(mask)
    threshold = float(mask["value"])
    operator = mask["operator"]
    discrete = _discrete_family(family)
    if operator == "==":
        if abs(threshold) <= 1.0e-12:
            return _model_p0(parameters, family)
        if discrete and abs(threshold - round(threshold)) <= 1.0e-12 and threshold >= 0.0:
            lower = _model_cdf(threshold, parameters, family)
            upper = _model_cdf(threshold - 1.0, parameters, family)
            return max(0.0, lower - upper)
        return 0.0
    if discrete:
        if operator == ">":
            return max(0.0, 1.0 - _model_cdf(math.floor(threshold), parameters, family))
        if operator == ">=":
            return max(0.0, 1.0 - _model_cdf(math.ceil(threshold) - 1.0, parameters, family))
        if operator == "<":
            return max(0.0, _model_cdf(math.ceil(threshold) - 1.0, parameters, family))
        return max(0.0, _model_cdf(math.floor(threshold), parameters, family))
    if operator in {">", ">="}:
        if threshold <= 0.0 and operator == ">=":
            return 1.0
        return max(0.0, 1.0 - _model_cdf(threshold, parameters, family))
    if operator in {"<", "<="}:
        return max(0.0, _model_cdf(threshold, parameters, family))
    return 0.0


def _response_fit_probability(mask, fit_range, parameters, family):
    """Probability for the observed mask intersected with the configured range."""
    lower, upper = [float(value) for value in fit_range]
    normalized = normalize_hgcer_mask(mask)
    if normalized["operator"] not in {">", ">="}:
        raise ValueError("response-sample mask must be a positive lower-tail selection")
    threshold = max(lower, float(normalized["value"]))
    if _discrete_family(family):
        minimum = math.floor(threshold) + 1 if normalized["operator"] == ">" else math.ceil(threshold)
        maximum = math.floor(upper)
        if maximum < minimum:
            return 0.0
        return max(0.0, _model_cdf(maximum, parameters, family) - _model_cdf(minimum - 1.0, parameters, family))
    if upper <= threshold:
        return 0.0
    return max(0.0, _model_cdf(upper, parameters, family) - _model_cdf(threshold, parameters, family))


def _initial_theta(values, family):
    mean = max(float(np.mean(values)), 0.05)
    if family == "zero_truncated_poisson":
        return np.log([max(mean, 0.2)])
    if family == "zero_truncated_negative_binomial":
        return np.log([max(mean, 0.2), 3.0])
    if family == "zero_truncated_compound_poisson_gamma":
        return np.log([1.0, 2.0, max(mean / 2.0, 0.05)])
    return np.log([1.0, max(mean, 0.05)])


def _parameter_bounds(family):
    count = len(_initial_theta(np.asarray([1.0]), family))
    return [(-9.0, 8.0)] * count


def _negative_log_likelihood(theta, values, response_mask, fit_range, family):
    try:
        parameters = _model_kind_parameters(theta, family)
        normalization = _response_fit_probability(response_mask, fit_range, parameters, family)
        if not math.isfinite(normalization) or normalization <= 0.0:
            return 1.0e100
        log_density = _model_log_probability(values, parameters, family)
        if not np.all(np.isfinite(log_density)):
            return 1.0e100
        return float(-np.sum(log_density) + len(values) * math.log(normalization))
    except Exception:
        return 1.0e100


def _numerical_hessian(function, point, *, step=1.0e-4):
    point = np.asarray(point, dtype=float)
    size = len(point)
    hessian = np.zeros((size, size), dtype=float)
    origin = float(function(point))
    for row in range(size):
        h_row = step * max(1.0, abs(point[row]))
        plus = point.copy(); plus[row] += h_row
        minus = point.copy(); minus[row] -= h_row
        hessian[row, row] = (function(plus) - 2.0 * origin + function(minus)) / (h_row * h_row)
        for column in range(row + 1, size):
            h_column = step * max(1.0, abs(point[column]))
            pp = point.copy(); pp[row] += h_row; pp[column] += h_column
            pm = point.copy(); pm[row] += h_row; pm[column] -= h_column
            mp = point.copy(); mp[row] -= h_row; mp[column] += h_column
            mm = point.copy(); mm[row] -= h_row; mm[column] -= h_column
            value = (function(pp) - function(pm) - function(mp) + function(mm)) / (4.0 * h_row * h_column)
            hessian[row, column] = value
            hessian[column, row] = value
    return hessian


def _covariance_diagnostics(function, optimum, config):
    hessian = _numerical_hessian(function, optimum)
    try:
        covariance = np.linalg.inv(hessian)
        eigenvalues = np.linalg.eigvalsh(covariance)
    except np.linalg.LinAlgError:
        return None, {"status": "singular", "near_degenerate": True}
    diagonal = np.diag(covariance)
    if np.any(~np.isfinite(covariance)) or np.any(diagonal <= 0.0) or np.any(eigenvalues <= 0.0):
        return None, {"status": "non_positive", "near_degenerate": True}
    standard = np.sqrt(diagonal)
    correlation = covariance / np.outer(standard, standard)
    maximum_pair = 0.0
    if len(correlation) > 1:
        maximum_pair = float(np.max(np.abs(correlation - np.eye(len(correlation)))))
    condition = float(np.linalg.cond(correlation))
    near_degenerate = bool(
        not math.isfinite(condition)
        or condition > float(config.get("correlation_condition_number_max", 1.0e8))
        or maximum_pair >= float(config.get("pair_correlation_abs_max", 0.995))
    )
    return covariance, {
        "status": "finite",
        "hessian": hessian.tolist(),
        "covariance": covariance.tolist(),
        "correlation": correlation.tolist(),
        "correlation_condition_number": condition,
        "maximum_pair_correlation": maximum_pair,
        "near_degenerate": near_degenerate,
    }


def _profile_mu(function, optimum, bounds, config):
    optimum = np.asarray(optimum, dtype=float)
    span = float(config.get("profile_scan_log_mu_half_width", 3.0))
    points = int(config.get("profile_scan_points", 31))
    grid = np.linspace(optimum[0] - span, optimum[0] + span, points)
    profile = []
    for log_mu in grid:
        if len(optimum) == 1:
            nll = float(function(np.asarray([log_mu])))
        else:
            def reduced(rest):
                full = optimum.copy()
                full[0] = log_mu
                full[1:] = rest
                return function(full)
            fitted = minimize(reduced, optimum[1:], method="L-BFGS-B", bounds=bounds[1:])
            nll = float(fitted.fun) if fitted.success and math.isfinite(float(fitted.fun)) else float("inf")
        profile.append((float(log_mu), nll))
    nll_values = np.asarray([entry[1] for entry in profile], dtype=float)
    if not np.any(np.isfinite(nll_values)):
        return {"status": "failed", "points": [{"log_mu": x, "nll": y} for x, y in profile]}
    baseline = float(function(optimum))
    target = baseline + float(config.get("profile_delta_nll", 0.5))
    best_index = int(np.nanargmin(nll_values))
    def crossing(indices):
        previous = best_index
        for index in indices:
            left_nll = nll_values[previous] - target
            right_nll = nll_values[index] - target
            if math.isfinite(left_nll) and math.isfinite(right_nll) and left_nll * right_nll <= 0.0:
                left_x, right_x = grid[previous], grid[index]
                if abs(right_nll - left_nll) <= 1.0e-15:
                    return float((left_x + right_x) / 2.0)
                return float(left_x + (right_x - left_x) * (-left_nll) / (right_nll - left_nll))
            previous = index
        return None
    left = crossing(range(best_index - 1, -1, -1))
    right = crossing(range(best_index + 1, len(grid)))
    status = "two_sided_68pct" if left is not None and right is not None else ("one_sided" if left is not None or right is not None else "unbounded")
    return {
        "status": status,
        "log_mu_low": left,
        "log_mu_high": right,
        "delta_nll": float(config.get("profile_delta_nll", 0.5)),
        "points": [{"log_mu": float(x), "nll": None if not math.isfinite(y) else float(y)} for x, y in profile],
    }


def _response_toys(theta, covariance, family, kaon_mask, control_mask, count, seed):
    rng = np.random.default_rng(int(seed) % (2 ** 63 - 1))
    draws = rng.multivariate_normal(np.asarray(theta, dtype=float), np.asarray(covariance, dtype=float), size=int(count))
    rows = []
    for draw in draws:
        try:
            parameters = _model_kind_parameters(draw, family)
            p0 = _model_p0(parameters, family)
            epsilon_k = model_mask_probability(kaon_mask, parameters, family)
            epsilon_pi = model_mask_probability(control_mask, parameters, family)
            transfer = epsilon_k / epsilon_pi
            if all(math.isfinite(value) and value >= 0.0 for value in (p0, epsilon_k, epsilon_pi, transfer)) and epsilon_pi > 0.0:
                rows.append((p0, epsilon_k, epsilon_pi, transfer))
        except Exception:
            continue
    return np.asarray(rows, dtype=float)


def _interval_summary(values):
    if values is None or len(values) == 0:
        return {"central": None, "p16": None, "p84": None, "sigma": None}
    p16, p84 = np.percentile(values, [16.0, 84.0])
    return {
        "central": float(np.median(values)),
        "p16": float(p16),
        "p84": float(p84),
        "sigma": float((p84 - p16) / 2.0),
    }


def _invalid_fit(reason, *, family=None, n_records=0, excluded_records=0):
    return {
        "fit_status": "fit_invalid",
        "reason": str(reason),
        "response_family": family,
        "prompt_pion_records": int(n_records),
        "excluded_fit_range_records": int(excluded_records),
        "R_pi_to_K": None,
    }


def fit_zero_truncated_pion_response(
    values,
    *,
    response_mask,
    kaon_mask,
    physical_control_mask,
    config,
    fingerprint_context="",
):
    """Fit a zero-atom response and infer transfer probabilities from toys."""
    config = dict(config or {})
    finite_values = np.asarray([float(value) for value in values if _finite(value) is not None], dtype=float)
    lower, upper = [float(value) for value in config.get("fit_range", (0.0, 20.0))]
    fit_values = finite_values[(finite_values >= lower) & (finite_values <= upper)]
    if len(fit_values) < int(config.get("minimum_prompt_pion_records", 20)):
        return _invalid_fit("insufficient_prompt_pion_records", n_records=len(fit_values), excluded_records=len(finite_values) - len(fit_values))
    try:
        primary_family, alternate_family, integer_like = choose_pion_hgcer_response_family(fit_values, config)
    except Exception as exc:
        return _invalid_fit(exc, n_records=len(fit_values), excluded_records=len(finite_values) - len(fit_values))

    def fit_family(family):
        initial = _initial_theta(fit_values, family)
        bounds = _parameter_bounds(family)
        function = lambda theta: _negative_log_likelihood(theta, fit_values, response_mask, (lower, upper), family)
        result = minimize(function, initial, method="L-BFGS-B", bounds=bounds)
        if not result.success or not math.isfinite(float(result.fun)):
            return _invalid_fit("optimizer_failed:{}".format(result.message), family=family, n_records=len(fit_values), excluded_records=len(finite_values) - len(fit_values))
        theta = np.asarray(result.x, dtype=float)
        bound_margin = float(config.get("fit_parameter_bound_margin", 1.0e-5))
        if any(
            theta[index] <= float(bound[0]) + bound_margin
            or theta[index] >= float(bound[1]) - bound_margin
            for index, bound in enumerate(bounds)
        ):
            return _invalid_fit(
                "response_parameter_at_bound", family=family,
                n_records=len(fit_values),
                excluded_records=len(finite_values) - len(fit_values),
            )
        try:
            parameters = _model_kind_parameters(theta, family)
            epsilon_k = model_mask_probability(kaon_mask, parameters, family)
            epsilon_pi = model_mask_probability(physical_control_mask, parameters, family)
            p0 = _model_p0(parameters, family)
        except Exception as exc:
            return _invalid_fit(exc, family=family, n_records=len(fit_values), excluded_records=len(finite_values) - len(fit_values))
        if not all(math.isfinite(value) and value >= 0.0 for value in (p0, epsilon_k, epsilon_pi)) or epsilon_pi <= 0.0:
            return _invalid_fit("nonphysical_efficiency", family=family, n_records=len(fit_values), excluded_records=len(finite_values) - len(fit_values))
        covariance, covariance_diagnostics = _covariance_diagnostics(function, theta, config)
        if covariance is None:
            return _invalid_fit("covariance_{}".format(covariance_diagnostics.get("status")), family=family, n_records=len(fit_values), excluded_records=len(finite_values) - len(fit_values))
        profile = _profile_mu(function, theta, bounds, config)
        seed = int(_fingerprint({"context": fingerprint_context, "family": family, "theta": theta.tolist()})[:16], 16)
        toys = _response_toys(theta, covariance, family, kaon_mask, physical_control_mask, config.get("toy_count", 2000), seed)
        accepted_fraction = float(len(toys)) / float(config.get("toy_count", 2000))
        p0_summary = _interval_summary(toys[:, 0] if len(toys) else [])
        transfer_summary = _interval_summary(toys[:, 3] if len(toys) else [])
        p0_relative_uncertainty = (
            float(p0_summary["sigma"] / p0) if p0 > 0.0 and p0_summary["sigma"] is not None else None
        )
        toy_quality_failed = accepted_fraction < float(
            config.get("minimum_accepted_toy_fraction", 0.50)
        )
        weak = bool(
            covariance_diagnostics.get("near_degenerate")
            or profile.get("status") != "two_sided_68pct"
        )
        return {
            "fit_status": (
                "fit_invalid" if toy_quality_failed
                else "fit_valid_but_P0_weakly_constrained" if weak else "fit_valid"
            ),
            "reason": (
                "accepted_toy_fraction_below_minimum" if toy_quality_failed
                else "P0_identifiability_weak" if weak else ""
            ),
            "response_family": family,
            "integer_like": bool(integer_like),
            "prompt_pion_records": int(len(fit_values)),
            "excluded_fit_range_records": int(len(finite_values) - len(fit_values)),
            "fit_range": [lower, upper],
            "nll": float(result.fun),
            "parameters": parameters,
            "log_parameters": theta.tolist(),
            "covariance_diagnostics": covariance_diagnostics,
            "profile_likelihood": profile,
            "P0": float(p0),
            "P0_statistical_uncertainty": p0_summary["sigma"],
            "P0_toy_interval_16_84": [p0_summary["p16"], p0_summary["p84"]],
            "P0_relative_uncertainty": p0_relative_uncertainty,
            "accepted_toy_fraction": accepted_fraction,
            "epsilon_pi_to_K": float(epsilon_k),
            "epsilon_pi_to_pi_control": float(epsilon_pi),
            "R_pi_to_K": float(epsilon_k / epsilon_pi),
            "sigma_stat": transfer_summary["sigma"],
            "transfer_toy_interval_16_84": [transfer_summary["p16"], transfer_summary["p84"]],
            "toy_count": int(config.get("toy_count", 2000)),
            "_theta": theta,
            "_covariance": covariance,
        }

    primary = fit_family(primary_family)
    alternate = fit_family(alternate_family)
    primary["alternate_fit"] = {
        key: value for key, value in alternate.items() if not str(key).startswith("_")
    }
    if primary.get("fit_status") in {"fit_valid", "fit_valid_but_P0_weakly_constrained"} and alternate.get("fit_status") in {"fit_valid", "fit_valid_but_P0_weakly_constrained"}:
        primary["sigma_model"] = abs(float(primary["R_pi_to_K"]) - float(alternate["R_pi_to_K"]))
    else:
        primary["sigma_model"] = None
    primary["sigma_total"] = (
        math.sqrt(float(primary["sigma_stat"]) ** 2 + float(primary["sigma_model"]) ** 2)
        if primary.get("sigma_stat") is not None and primary.get("sigma_model") is not None
        else None
    )
    return primary


def _cell_key(cell):
    return int(cell["t_index"]), int(cell["delta_index"])


def _source_cells(part1_payload):
    return {_cell_key(cell): dict(cell) for cell in list((part1_payload or {}).get("cells") or ())}


def _unavailable(reason, stage, *, config=None, provenance=None, part1_payload=None):
    return {
        "status": "unavailable",
        "diagnostic_label": TRANSFER_LABEL,
        "non_authoritative": True,
        "production_side_effect_free": True,
        "reason": str(reason),
        "diagnostic_stage": str(stage),
        "config": deepcopy(config or {}),
        "source_provenance": deepcopy(provenance or {}),
        "coordinate_fingerprint": (part1_payload or {}).get("coordinate_fingerprint"),
        "t_edges": list((part1_payload or {}).get("t_edges") or ()),
        "delta_edges": list((part1_payload or {}).get("delta_edges") or ()),
    }


def _prompt_values_by_cell(part1_payload, response_mask):
    grouped = {}
    for record in list(((part1_payload or {}).get("records") or {}).get("pion") or ()):
        if str(record.get("source_label")) != "prompt" or not bool(record.get("allcuts", False)):
            continue
        value = _finite(record.get(_MASK_FIELD))
        if value is None or not hgcer_mask_accepts(response_mask, value):
            continue
        key = (int(record.get("canonical_t_index")), int(record.get("delta_index")))
        grouped.setdefault(key, []).append(float(value))
    return grouped


def _fallback_sigma(*values):
    finite = [float(value) for value in values if _finite(value) is not None]
    return math.sqrt(sum(value * value for value in finite)) if finite else None


def _apply_delta_interpolation(cells, t_edges, delta_edges, config):
    fallback = dict(config.get("fallback") or {})
    n_t = max(0, len(t_edges) - 1)
    n_delta = max(0, len(delta_edges) - 1)
    for t_index in range(n_t):
        for delta_index in range(n_delta):
            cell = cells[(t_index, delta_index)]
            if cell.get("solution_source") != "unresolved" or cell.get("support_class") != "supported":
                continue
            left = cells.get((t_index, delta_index - 1))
            right = cells.get((t_index, delta_index + 1))
            valid_left = left if left and left.get("solution_source") == "direct" else None
            valid_right = right if right and right.get("solution_source") == "direct" else None
            if valid_left and valid_right:
                center = (float(delta_edges[delta_index]) + float(delta_edges[delta_index + 1])) / 2.0
                left_center = (float(delta_edges[delta_index - 1]) + float(delta_edges[delta_index])) / 2.0
                right_center = (float(delta_edges[delta_index + 1]) + float(delta_edges[delta_index + 2])) / 2.0
                fraction = (center - left_center) / (right_center - left_center)
                transfer = (1.0 - fraction) * float(valid_left["R_pi_to_K"]) + fraction * float(valid_right["R_pi_to_K"])
                statistical = _fallback_sigma(
                    (1.0 - fraction) * float(valid_left.get("sigma_stat") or 0.0),
                    fraction * float(valid_right.get("sigma_stat") or 0.0),
                )
                model = _fallback_sigma(
                    (1.0 - fraction) * float(valid_left.get("sigma_model") or 0.0),
                    fraction * float(valid_right.get("sigma_model") or 0.0),
                )
                fallback_sigma = abs(transfer) * float(fallback.get("delta_interpolation_relative_uncertainty", 0.20))
                cell.update({"R_pi_to_K": transfer, "sigma_stat": statistical, "sigma_model": model, "sigma_fallback": fallback_sigma, "solution_source": "delta_interpolation", "fallback_source": ["t{}d{}".format(t_index + 1, delta_index), "t{}d{}".format(t_index + 1, delta_index + 2)]})
            elif (delta_index == 0 and valid_right) or (delta_index == n_delta - 1 and valid_left):
                neighbor = valid_right or valid_left
                transfer = float(neighbor["R_pi_to_K"])
                cell.update({"R_pi_to_K": transfer, "sigma_stat": neighbor.get("sigma_stat"), "sigma_model": neighbor.get("sigma_model"), "sigma_fallback": abs(transfer) * float(fallback.get("edge_neighbor_relative_uncertainty", 0.30)), "solution_source": "delta_interpolation", "fallback_source": ["edge_neighbor"]})


def _finish_uncertainties(cell):
    pieces = [cell.get("sigma_stat"), cell.get("sigma_model"), cell.get("sigma_fallback")]
    finite = [float(value) for value in pieces if _finite(value) is not None]
    cell["sigma_total"] = math.sqrt(sum(value * value for value in finite)) if finite else None


def build_pion_hgcer_transfer_map(
    *,
    part1_payload,
    pion_control_cache,
    source_entries,
    config,
):
    """Build/freeze the Part-2 transfer map without touching production state."""
    config = deepcopy(config or {})
    if str((part1_payload or {}).get("status") or "unavailable") != "available":
        return _unavailable("Part-1 HGCer diagnostic is unavailable", "part1_input", config=config, part1_payload=part1_payload)
    physical_control_mask = (pion_control_cache or {}).get("physical_control_hgcer_mask")
    if physical_control_mask is None:
        return _unavailable("pion-control cache has no explicit physical HGCer mask", "physical_control_mask", config=config, part1_payload=part1_payload)
    try:
        selection = resolve_hgcer_transfer_masks(source_entries, physical_control_mask)
    except Exception as exc:
        return _unavailable(exc, "pid_mask_provenance", config=config, part1_payload=part1_payload)
    t_edges = [float(value) for value in list(part1_payload.get("t_edges") or ())]
    delta_edges = [float(value) for value in list(part1_payload.get("delta_edges") or ())]
    if len(t_edges) < 2 or len(delta_edges) < 2:
        return _unavailable("Part-1 geometry is incomplete", "geometry", config=config, provenance=selection, part1_payload=part1_payload)
    source_cells = _source_cells(part1_payload)
    prompt_values = _prompt_values_by_cell(
        part1_payload, selection["masks"]["S_pi_response_sample"]
    )
    cells = {}
    for key, source_cell in source_cells.items():
        t_index, delta_index = key
        cell = {
            "t_index": t_index,
            "t_low": float(source_cell["t_low"]),
            "t_high": float(source_cell["t_high"]),
            "delta_index": delta_index,
            "delta_low": float(source_cell["delta_low"]),
            "delta_high": float(source_cell["delta_high"]),
            "full_cell_support": {"kaon": deepcopy(source_cell.get("kaon") or {}), "pion": deepcopy(source_cell.get("pion") or {})},
            "support_class": str(source_cell.get("support_class") or "unsupported"),
            "prompt_pion_records": int(len(prompt_values.get(key, ()))),
            "fit_status": "not_attempted",
            "solution_source": "unresolved",
            "fallback_source": [],
            "R_pi_to_K": None,
            "sigma_stat": None,
            "sigma_model": None,
            "sigma_fallback": None,
            "sigma_total": None,
        }
        if cell["support_class"] == "supported":
            fit = fit_zero_truncated_pion_response(
                prompt_values.get(key, ()),
                response_mask=selection["masks"]["S_pi_response_sample"],
                kaon_mask=selection["masks"]["S_K_tree"],
                physical_control_mask=selection["masks"]["S_pi_physical_control"],
                config=config,
                fingerprint_context="{}:{}:{}".format(selection["pid_selection_fingerprint"], t_index, delta_index),
            )
            cell.update({key: value for key, value in fit.items() if not str(key).startswith("_")})
            if fit.get("fit_status") == "fit_valid":
                cell["solution_source"] = "direct"
            elif fit.get("fit_status") == "fit_valid_but_P0_weakly_constrained":
                cell["solution_source"] = "unresolved"
        cells[key] = cell
    _apply_delta_interpolation(cells, t_edges, delta_edges, config)
    fallback = dict(config.get("fallback") or {})
    for t_index in range(len(t_edges) - 1):
        unresolved = [cell for (index, _), cell in cells.items() if index == t_index and cell.get("solution_source") == "unresolved" and cell.get("support_class") == "supported"]
        if not unresolved:
            continue
        pooled_values = [value for (index, _), values in prompt_values.items() if index == t_index for value in values]
        if len(pooled_values) < int(fallback.get("minimum_pooled_prompt_pion_records", 20)):
            continue
        pooled = fit_zero_truncated_pion_response(
            pooled_values,
            response_mask=selection["masks"]["S_pi_response_sample"],
            kaon_mask=selection["masks"]["S_K_tree"],
            physical_control_mask=selection["masks"]["S_pi_physical_control"],
            config=config,
            fingerprint_context="{}:t{}:pooled".format(selection["pid_selection_fingerprint"], t_index),
        )
        if pooled.get("fit_status") != "fit_valid":
            continue
        for cell in unresolved:
            cell.update({
                "R_pi_to_K": pooled.get("R_pi_to_K"),
                "sigma_stat": pooled.get("sigma_stat"),
                "sigma_model": pooled.get("sigma_model"),
                "sigma_fallback": abs(float(pooled.get("R_pi_to_K") or 0.0)) * float(fallback.get("edge_neighbor_relative_uncertainty", 0.30)),
                "solution_source": "t_pooled",
                "fallback_source": ["t{} pooled".format(t_index + 1)],
                "pooled_fit": {key: value for key, value in pooled.items() if not str(key).startswith("_")},
            })
    for cell in cells.values():
        if cell.get("solution_source") == "unresolved":
            cell["solution_source"] = "unsupported"
        _finish_uncertainties(cell)
    serial_cells = [cells[key] for key in sorted(cells)]
    map_essentials = [{key: value for key, value in cell.items() if key not in {"full_cell_support", "pooled_fit"}} for cell in serial_cells]
    transfer_fingerprint = _fingerprint({
        "t_edges": t_edges,
        "delta_edges": delta_edges,
        "coordinate_fingerprint": part1_payload.get("coordinate_fingerprint"),
        "pid_selection_fingerprint": selection["pid_selection_fingerprint"],
        "config": config,
        "cells": map_essentials,
    })
    payload = {
        "status": "available",
        "diagnostic_label": TRANSFER_LABEL,
        "non_authoritative": True,
        "production_side_effect_free": True,
        "coordinate_fingerprint": part1_payload.get("coordinate_fingerprint"),
        "t_edges": t_edges,
        "delta_edges": delta_edges,
        "config": config,
        "config_fingerprint": _fingerprint(config),
        "masks": selection["masks"],
        "response_sample_provenance": selection["response_sample_provenance"],
        "source_provenance": selection["source_provenance"],
        "pid_selection_fingerprint": selection["pid_selection_fingerprint"],
        "transfer_map_fingerprint": transfer_fingerprint,
        "cells": serial_cells,
        "closure": {},
    }
    return _freeze(payload)


def _detached_root_histogram(source, name, *, reset=False):
    """Clone a ROOT histogram into diagnostic ownership, if ROOT is present."""
    if ROOT is None or source is None:
        return None
    clone = source.Clone(str(name))
    clone.SetDirectory(0)
    if reset:
        clone.Reset("ICES")
    clone.Sumw2()
    return clone


def _root_histogram_integral(histogram):
    if histogram is None:
        return None
    return float(histogram.Integral(0, histogram.GetNbinsX() + 1))


def _root_histogram_max_difference(left, right):
    if left is None or right is None:
        return None
    return max(
        abs(float(left.GetBinContent(index)) - float(right.GetBinContent(index)))
        for index in range(0, left.GetNbinsX() + 2)
    )


def _build_detached_proposed_mm_products(application, pion_control_cache, proton_cleaning_application):
    """Build proposal-only pion and kaon MM products from already frozen inputs.

    The helper deliberately clones both the host and the signed pion estimate;
    no histogram in the authoritative random/dummy, proton-cleaning, or pion
    subtraction chain is ever modified.
    """
    if ROOT is None or not isinstance(proton_cleaning_application, Mapping):
        return {}, {}
    def index_t_products(products):
        if isinstance(products, Mapping):
            return {int(key): value for key, value in products.items()}
        return {
            int(product.get("t_index")): product
            for product in (products or ())
            if isinstance(product, Mapping) and product.get("t_index") is not None
        }

    canonical_products = index_t_products(
        proton_cleaning_application.get("canonical_t_products")
    )
    cache_by_t = index_t_products((pion_control_cache or {}).get("by_t"))
    suffix = str(application.get("transfer_map_fingerprint") or "detached")[:12]
    products = {}
    summaries = {}
    for t_index, summary in (application.get("per_t") or {}).items():
        key = int(t_index)
        source_cache = cache_by_t.get(key, {}) or {}
        source_product = canonical_products.get(key, {}) or {}
        pion_template = source_cache.get("H_pion_control_cut")
        host_template = ((source_product.get("final_targets") or {}).get("h_mm"))
        if pion_template is None or host_template is None:
            summaries[key] = {
                "status": "unavailable",
                "reason": "missing_detached_t_histogram",
            }
            continue
        pion_estimate = _detached_root_histogram(
            pion_template, "H_pion_hgcer_transfer_proposed_t{}_{}".format(key + 1, suffix), reset=True
        )
        host = _detached_root_histogram(
            host_template, "H_kaon_hgcer_transfer_host_t{}_{}".format(key + 1, suffix)
        )
        cleaned = _detached_root_histogram(
            host_template, "H_kaon_hgcer_transfer_cleaned_t{}_{}".format(key + 1, suffix)
        )
        for record in application.get("records") or ():
            if int(record.get("t_index", -1)) != key:
                continue
            mm_value = _finite(record.get("adj_MM"))
            if mm_value is not None:
                pion_estimate.Fill(mm_value, float(record["weight"]))
        cleaned.Add(pion_estimate, -1.0)
        oversubtraction_bins = 0
        min_cleaned_bin = None
        for bin_index in range(1, cleaned.GetNbinsX() + 1):
            cleaned_value = float(cleaned.GetBinContent(bin_index))
            pion_value = float(pion_estimate.GetBinContent(bin_index))
            host_value = float(host.GetBinContent(bin_index))
            if pion_value > host_value:
                oversubtraction_bins += 1
            min_cleaned_bin = cleaned_value if min_cleaned_bin is None else min(min_cleaned_bin, cleaned_value)
        products[key] = {
            "pion_estimate_mm": pion_estimate,
            "kaon_host_mm": host,
            "kaon_cleaned_mm": cleaned,
        }
        summaries[key] = {
            "status": "available",
            "pion_estimate_integral": _root_histogram_integral(pion_estimate),
            "host_integral": _root_histogram_integral(host),
            "cleaned_integral": _root_histogram_integral(cleaned),
            "oversubtraction_bin_count": int(oversubtraction_bins),
            "minimum_cleaned_bin": min_cleaned_bin,
            "record_count": int(summary.get("record_count", 0)),
        }
    if not products:
        return products, summaries
    first_product = next(iter(products.values()))
    global_pion = _detached_root_histogram(first_product["pion_estimate_mm"], "H_pion_hgcer_transfer_proposed_global_{}".format(suffix), reset=True)
    global_host = _detached_root_histogram(first_product["kaon_host_mm"], "H_kaon_hgcer_transfer_host_global_{}".format(suffix), reset=True)
    global_cleaned = _detached_root_histogram(first_product["kaon_cleaned_mm"], "H_kaon_hgcer_transfer_cleaned_global_{}".format(suffix), reset=True)
    for product in products.values():
        global_pion.Add(product["pion_estimate_mm"])
        global_host.Add(product["kaon_host_mm"])
        global_cleaned.Add(product["kaon_cleaned_mm"])
    setting_host = ((proton_cleaning_application.get("final_targets") or {}).get("h_mm"))
    host_difference = _root_histogram_max_difference(global_host, setting_host)
    host_scale = max(abs(_root_histogram_integral(global_host) or 0.0), 1.0)
    strict_global_sums = host_difference is not None and host_difference <= 1.0e-9 * host_scale
    products["global"] = {
        "pion_estimate_mm": global_pion,
        "kaon_host_mm": global_host,
        "kaon_cleaned_mm": global_cleaned,
    }
    summaries["global"] = {
        "status": "available",
        "pion_estimate_integral": _root_histogram_integral(global_pion),
        "host_integral": _root_histogram_integral(global_host),
        "cleaned_integral": _root_histogram_integral(global_cleaned),
        "host_global_max_bin_difference": host_difference,
        "strict_global_sums": bool(strict_global_sums),
    }
    return products, summaries


def apply_pion_hgcer_transfer_diagnostic(
    transfer_map, pion_control_cache, proton_cleaning_application=None
):
    """Apply only the frozen map to detached signed pion-control accounting."""
    if str((transfer_map or {}).get("status") or "unavailable") != "available":
        return {"status": "unavailable", "reason": "transfer map unavailable", "non_authoritative": True}
    lookup = {
        (int(cell["t_index"]), int(cell["delta_index"])): cell
        for cell in list(transfer_map.get("cells") or ())
    }
    delta_edges = [float(value) for value in list(transfer_map.get("delta_edges") or ())]
    per_t = {}
    applied_records = []
    for record in list((pion_control_cache or {}).get("records") or ()):
        if not bool(record.get("allcuts", False)):
            continue
        delta_value = _finite(record.get("ssdelta"))
        if delta_value is None:
            continue
        delta_index = find_canonical_bin(delta_value, tuple(delta_edges))
        t_index = int(record.get("t_index", -1))
        cell = lookup.get((t_index, int(delta_index))) if delta_index >= 0 else None
        if cell is None or cell.get("solution_source") == "unsupported" or _finite(cell.get("R_pi_to_K")) is None:
            continue
        coefficient = float(record.get("coefficient", 0.0) or 0.0)
        weight = coefficient * float(cell["R_pi_to_K"])
        summary = per_t.setdefault(t_index, {"pion_estimate_integral": 0.0, "record_count": 0})
        summary["pion_estimate_integral"] += weight
        summary["record_count"] += 1
        applied_records.append({
            "t_index": t_index,
            "delta_index": int(delta_index),
            "adj_MM": _finite(record.get("adj_MM")),
            "coefficient": coefficient,
            "transfer": float(cell["R_pi_to_K"]),
            "weight": weight,
            "proton_cleaning_factor": None,
        })
    result = {
        "status": "available",
        "non_authoritative": True,
        "transfer_map_fingerprint": transfer_map.get("transfer_map_fingerprint"),
        "per_t": per_t,
        "records": applied_records,
        "source_weight_contract": "signed pion coefficient times frozen R(t,delta); no proton factor; no SIMC MM weight",
    }
    root_products, root_summary = _build_detached_proposed_mm_products(
        result, pion_control_cache, proton_cleaning_application
    )
    if root_products:
        result["_root_products"] = root_products
    if root_summary:
        result["mm_products"] = root_summary
        global_summary = root_summary.get("global") or {}
        if global_summary and not bool(global_summary.get("strict_global_sums", False)):
            result["status"] = "unavailable"
            result["reason"] = "strict_global_kaon_host_sum_failed"
    return _freeze(result)


def summarize_pion_hgcer_transfer_closure(application, closure_inputs):
    """Summarize closure-only reference shapes without feeding any fit or weight.

    The returned integrals are deliberately descriptive.  SIMC and protected
    K-Lambda/K-Sigma0 shapes remain validation references, never PID inputs or
    normalizers for the response transfer.
    """
    root_products = (application or {}).get("_root_products") or {}
    proposed = (root_products.get("global") or {}).get("kaon_cleaned_mm")
    if ROOT is None or proposed is None:
        return {
            "status": "unavailable",
            "reason": "detached_global_MM_product_unavailable",
            "authoritative": False,
        }
    rows = {}
    for label, histogram in (closure_inputs or {}).items():
        if histogram is None:
            rows[str(label)] = {"status": "unavailable"}
            continue
        reference_integral = _root_histogram_integral(histogram)
        proposed_integral = _root_histogram_integral(proposed)
        rows[str(label)] = {
            "status": "available",
            "reference_integral": reference_integral,
            "proposed_cleaned_kaon_integral": proposed_integral,
            "integral_ratio_to_proposed": (
                proposed_integral / reference_integral
                if reference_integral not in (None, 0.0) else None
            ),
        }
    return {"status": "available", "authoritative": False, "references": rows}


def serialize_pion_hgcer_transfer_map(
    payload, include_records=False, *, application=None, page_manifest=None, closure=None
):
    payload = payload or {}
    result = {
        "status": payload.get("status", "unavailable"),
        "reason": payload.get("reason"),
        "diagnostic_stage": payload.get("diagnostic_stage"),
        "diagnostic_label": payload.get("diagnostic_label", TRANSFER_LABEL),
        "non_authoritative": True,
        "production_side_effect_free": True,
        "coordinate_fingerprint": payload.get("coordinate_fingerprint"),
        "t_edges": list(payload.get("t_edges") or ()),
        "delta_edges": list(payload.get("delta_edges") or ()),
        "config": payload.get("config") or {},
        "config_fingerprint": payload.get("config_fingerprint"),
        "masks": payload.get("masks") or {},
        "response_sample_provenance": payload.get("response_sample_provenance") or {},
        "source_provenance": payload.get("source_provenance") or {},
        "pid_selection_fingerprint": payload.get("pid_selection_fingerprint"),
        "transfer_map_fingerprint": payload.get("transfer_map_fingerprint"),
        "cells": list(payload.get("cells") or ()),
        "closure": closure if closure is not None else (payload.get("closure") or {}),
        "page_manifest": list(page_manifest if page_manifest is not None else (payload.get("page_manifest") or ())),
    }
    if application is not None:
        result["proposed_application"] = {
            "status": application.get("status", "unavailable"),
            "non_authoritative": True,
            "transfer_map_fingerprint": application.get("transfer_map_fingerprint"),
            "per_t": application.get("per_t") or {},
            "mm_products": application.get("mm_products") or {},
            "source_weight_contract": application.get("source_weight_contract"),
        }
    if include_records:
        result["records"] = list(payload.get("records") or ())
    return _json_ready(result)


def write_pion_hgcer_transfer_json(
    path, payload, *, application=None, page_manifest=None, closure=None
):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(
            serialize_pion_hgcer_transfer_map(
                payload, application=application, page_manifest=page_manifest,
                closure=closure,
            ),
            handle, sort_keys=True, indent=2, allow_nan=False,
        )
    return path


def expected_pion_hgcer_transfer_page_manifest(payload):
    if str((payload or {}).get("status") or "unavailable") != "available":
        return [{"page_id": "pion_hgcer.part2.status", "scope": "setting-wide", "authoritative": False}]
    pages = [
        {"page_id": "pion_hgcer.part2.provenance", "scope": "setting-wide", "authoritative": False},
        {"page_id": "pion_hgcer.part2.response_family", "scope": "setting-wide", "authoritative": False},
        {"page_id": "pion_hgcer.part2.p0", "scope": "setting-wide", "authoritative": False},
        {"page_id": "pion_hgcer.part2.p0_relative_uncertainty", "scope": "setting-wide", "authoritative": False},
        {"page_id": "pion_hgcer.part2.transfer", "scope": "setting-wide", "authoritative": False},
        {"page_id": "pion_hgcer.part2.uncertainty", "scope": "setting-wide", "authoritative": False},
        {"page_id": "pion_hgcer.part2.status_map", "scope": "setting-wide", "authoritative": False},
        {"page_id": "pion_hgcer.part2.solution_source", "scope": "setting-wide", "authoritative": False},
    ]
    if bool((payload.get("config") or {}).get("render_direct_fit_cells", True)):
        for cell in sorted(
            payload.get("cells") or (),
            key=lambda item: (int(item.get("t_index", -1)), int(item.get("delta_index", -1))),
        ):
            if str(cell.get("fit_status")) in {
                "fit_valid", "fit_valid_but_P0_weakly_constrained"
            }:
                pages.append({
                    "page_id": "pion_hgcer.part2.response.t{}.d{}".format(
                        int(cell["t_index"]) + 1, int(cell["delta_index"]) + 1
                    ),
                    "scope": "t{}d{}".format(
                        int(cell["t_index"]) + 1, int(cell["delta_index"]) + 1
                    ),
                    "authoritative": False,
                })
    for t_index in range(max(0, len(payload.get("t_edges") or ()) - 1)):
        pages.append({"page_id": "pion_hgcer.part2.t{}".format(t_index + 1), "scope": "t{}".format(t_index + 1), "authoritative": False})
        pages.append({"page_id": "pion_hgcer.part2.mm.t{}".format(t_index + 1), "scope": "t{}".format(t_index + 1), "authoritative": False})
    pages.append({"page_id": "pion_hgcer.part2.closure", "scope": "setting-wide", "authoritative": False})
    pages.append({"page_id": "pion_hgcer.part2.mm.summary", "scope": "setting-wide", "authoritative": False})
    return pages


def _root_text_page(pdf_name, title, lines, *, close_pdf=False):
    canvas = ROOT.TCanvas("C_pion_hgcer_transfer_text", title, 1000, 800)
    box = ROOT.TPaveText(0.06, 0.08, 0.94, 0.92, "NDC")
    box.SetFillStyle(0)
    box.SetBorderSize(0)
    box.AddText(title)
    for line in lines:
        box.AddText(str(line))
    box.Draw()
    canvas.Print(str(pdf_name) + (")" if close_pdf else ""))
    return canvas


def render_pion_hgcer_transfer_pages(
    pdf_name, payload, *, title_prefix="", page_manifest=None, close_pdf=False,
    application=None, closure=None, closure_inputs=None,
):
    """Append compact Part-2 PDF pages without changing any ROOT histogram owner."""
    if ROOT is None or not isinstance(payload, Mapping):
        return []
    manifest = page_manifest if page_manifest is not None else []
    expected = expected_pion_hgcer_transfer_page_manifest(payload)
    if str(payload.get("status") or "unavailable") != "available":
        _root_text_page(pdf_name, TRANSFER_LABEL, [title_prefix, "status: unavailable", "reason: {}".format(payload.get("reason") or "not built")], close_pdf=close_pdf)
        manifest.extend(expected)
        return expected
    _root_text_page(pdf_name, TRANSFER_LABEL, [title_prefix, "S_K_tree: {}".format((payload.get("masks") or {}).get("S_K_tree", {}).get("expression")), "S_pi_physical_control: {}".format((payload.get("masks") or {}).get("S_pi_physical_control", {}).get("expression")), "transfer fingerprint: {}".format(payload.get("transfer_map_fingerprint"))])
    t_edges = [float(value) for value in payload.get("t_edges") or ()]
    delta_edges = [float(value) for value in payload.get("delta_edges") or ()]
    maps = (
        ("response_family", "response_family", "response family"),
        ("p0", "P0", "P0"),
        ("p0_relative_uncertainty", "P0_relative_uncertainty", "relative P0 uncertainty"),
        ("transfer", "R_pi_to_K", "R_pi_to_K"),
        ("uncertainty", "sigma_total", "total transfer uncertainty"),
        ("status_map", "fit_status", "fit status"),
        ("solution_source", "solution_source", "solution source"),
    )
    status_values = {"fit_invalid": 0.0, "not_attempted": 0.0, "fit_valid_but_P0_weakly_constrained": 1.0, "fit_valid": 2.0}
    response_values = {
        "zero_truncated_poisson": 1.0,
        "zero_truncated_negative_binomial": 2.0,
        "zero_truncated_compound_poisson_gamma": 3.0,
        "zero_truncated_compound_poisson_exponential": 4.0,
    }
    solution_values = {"unsupported": 0.0, "direct": 1.0, "delta_interpolation": 2.0, "t_pooled": 3.0}
    for index, (page_key, field, label) in enumerate(maps):
        histogram = ROOT.TH2D("H_pion_hgcer_part2_{}".format(page_key), "{};SHMS delta [%];analysis |t| [GeV^2]".format(label), len(delta_edges) - 1, array("d", delta_edges), len(t_edges) - 1, array("d", t_edges))
        histogram.SetDirectory(0)
        for cell in payload.get("cells") or ():
            if field == "fit_status":
                value = status_values.get(str(cell.get(field)), 0.0)
            elif field == "response_family":
                value = response_values.get(str(cell.get(field)), 0.0)
            elif field == "solution_source":
                value = solution_values.get(str(cell.get(field)), 0.0)
            else:
                value = _finite(cell.get(field), 0.0)
            histogram.SetBinContent(int(cell["delta_index"]) + 1, int(cell["t_index"]) + 1, float(value or 0.0))
        canvas = ROOT.TCanvas("C_pion_hgcer_part2_{}".format(page_key), "{} {}".format(TRANSFER_LABEL, label), 1000, 800)
        histogram.Draw("COLZ TEXT")
        canvas.Print(str(pdf_name))
    if bool((payload.get("config") or {}).get("render_direct_fit_cells", True)):
        lower, upper = [float(value) for value in (payload.get("config") or {}).get("fit_range", (0.0, 20.0))]
        for cell in sorted(
            payload.get("cells") or (),
            key=lambda item: (int(item.get("t_index", -1)), int(item.get("delta_index", -1))),
        ):
            if str(cell.get("fit_status")) not in {
                "fit_valid", "fit_valid_but_P0_weakly_constrained"
            }:
                continue
            label = "t{} delta{}".format(int(cell["t_index"]) + 1, int(cell["delta_index"]) + 1)
            canvas = ROOT.TCanvas(
                "C_pion_hgcer_part2_response_t{}_d{}".format(
                    int(cell["t_index"]) + 1, int(cell["delta_index"]) + 1
                ),
                "{} response/profile {}".format(TRANSFER_LABEL, label), 1200, 500,
            )
            canvas.Divide(2, 1)
            canvas.cd(1)
            density = ROOT.TH1D(
                "H_pion_hgcer_part2_response_t{}_d{}".format(
                    int(cell["t_index"]) + 1, int(cell["delta_index"]) + 1
                ),
                "positive response model {};NPE;relative density".format(label),
                80, lower, upper,
            )
            density.SetDirectory(0)
            parameters = cell.get("parameters") or {}
            family = cell.get("response_family")
            if parameters and family:
                try:
                    values = np.asarray([
                        density.GetBinCenter(bin_index)
                        for bin_index in range(1, density.GetNbinsX() + 1)
                    ])
                    log_probability = _model_log_probability(values, parameters, family)
                    normalization = _response_fit_probability(
                        (payload.get("masks") or {}).get("S_pi_response_sample"),
                        (lower, upper), parameters, family,
                    )
                    for bin_index, value in enumerate(log_probability, start=1):
                        density.SetBinContent(
                            bin_index,
                            math.exp(float(value)) / normalization
                            if math.isfinite(float(value)) and normalization > 0.0 else 0.0,
                        )
                except Exception:
                    pass
            density.Draw("HIST")
            canvas.cd(2)
            profile = (cell.get("profile_likelihood") or {}).get("points") or ()
            finite_profile = [
                (float(row["log_mu"]), float(row["nll"]))
                for row in profile if _finite(row.get("log_mu")) is not None and _finite(row.get("nll")) is not None
            ]
            if finite_profile:
                graph = ROOT.TGraph(
                    len(finite_profile),
                    array("d", [row[0] for row in finite_profile]),
                    array("d", [row[1] for row in finite_profile]),
                )
                graph.SetTitle("mu profile {};log(mu);NLL".format(label))
                graph.Draw("ALP")
            else:
                ROOT.TLatex(0.12, 0.55, "profile unavailable").Draw()
            canvas.Print(str(pdf_name))
    for t_index in range(max(0, len(t_edges) - 1)):
        rows = [cell for cell in payload.get("cells") or () if int(cell["t_index"]) == t_index]
        lines = ["{} canonical t{} [{:.4f}, {:.4f}]".format(title_prefix, t_index + 1, t_edges[t_index], t_edges[t_index + 1])]
        for cell in rows:
            lines.append("delta{} {} R={} P0={} status={}".format(int(cell["delta_index"]) + 1, cell.get("solution_source"), cell.get("R_pi_to_K"), cell.get("P0"), cell.get("fit_status")))
        _root_text_page(pdf_name, TRANSFER_LABEL, lines)
    root_products = (application or {}).get("_root_products") or {}
    for t_index in range(max(0, len(t_edges) - 1)):
        product = root_products.get(t_index, root_products.get(str(t_index))) or {}
        if not product:
            _root_text_page(
                pdf_name,
                TRANSFER_LABEL,
                ["proposed MM t{} unavailable".format(t_index + 1)],
            )
            continue
        canvas = ROOT.TCanvas(
            "C_pion_hgcer_part2_mm_t{}".format(t_index + 1),
            "{} proposed MM t{}".format(TRANSFER_LABEL, t_index + 1), 1350, 500,
        )
        canvas.Divide(3, 1)
        for pad, (name, label) in enumerate((
            ("kaon_host_mm", "proton-cleaned kaon host"),
            ("pion_estimate_mm", "proposed signed pion estimate"),
            ("kaon_cleaned_mm", "proposed pion-cleaned kaon"),
        ), start=1):
            canvas.cd(pad)
            histogram = product.get(name)
            if histogram is not None:
                histogram.SetTitle("{};MM [GeV];weighted events".format(label))
                histogram.Draw("E")
        canvas.Print(str(pdf_name))
    global_mm = root_products.get("global") or {}
    closure_rows = ((closure or {}).get("references") or {})
    proposed_cleaned = global_mm.get("kaon_cleaned_mm")
    closure_canvas = ROOT.TCanvas(
        "C_pion_hgcer_part2_closure", "{} closure references".format(TRANSFER_LABEL), 1350, 500
    )
    closure_canvas.Divide(3, 1)
    closure_specs = (
        ("pion_SIMC", "SIMC pion closure"),
        ("K_Lambda", "K-Lambda retention closure"),
        ("K_Sigma0", "K-Sigma0 retention closure"),
    )
    for pad, (key, label) in enumerate(closure_specs, start=1):
        closure_canvas.cd(pad)
        reference = (closure_inputs or {}).get(key)
        if proposed_cleaned is not None:
            proposed_draw = _detached_root_histogram(
                proposed_cleaned, "H_pion_hgcer_part2_closure_proposed_{}_{}".format(key, str(payload.get("transfer_map_fingerprint") or "")[:10])
            )
            proposed_draw.SetLineColor(ROOT.kBlack)
            proposed_draw.SetTitle("{};MM [GeV];weighted events".format(label))
            proposed_draw.Draw("HIST")
        if reference is not None:
            reference_draw = _detached_root_histogram(
                reference, "H_pion_hgcer_part2_closure_reference_{}_{}".format(key, str(payload.get("transfer_map_fingerprint") or "")[:10])
            )
            reference_draw.SetLineColor(ROOT.kBlue + 1)
            reference_draw.Draw("HIST SAME")
        elif proposed_cleaned is None:
            ROOT.TLatex(0.12, 0.55, "closure unavailable").Draw()
        row = closure_rows.get(key) or {}
        if row:
            ROOT.TLatex(0.12, 0.82, "status: {}".format(row.get("status"))).Draw()
    closure_canvas.Print(str(pdf_name))
    summary_lines = [
        "Part 2 map frozen; proposed MM products are detached.",
        "SIMC and K-Lambda/K-Sigma0 remain closure-only references.",
    ]
    global_summary = ((application or {}).get("mm_products") or {}).get("global") or {}
    if global_summary:
        summary_lines.append("strict global host sum: {}".format(global_summary.get("strict_global_sums")))
        summary_lines.append("global proposed pion integral: {}".format(global_summary.get("pion_estimate_integral")))
    if global_mm:
        summary_lines.append("global MM proposal is available for closure inspection")
    _root_text_page(pdf_name, TRANSFER_LABEL, summary_lines, close_pdf=close_pdf)
    manifest.extend(expected)
    return expected
