#!/usr/bin/env python3
"""Read retained pion-chain JSON; never fit, correct, or accept production.

Usage: python testing/analyze_pion_transfer_validity.py audit.json [--report new.json]
       [--bins name=1,2,4]

Bins are one-based native normal bins. All normal bins and stored metadata
windows are reported by default. Flows are checked but never mixed into regions.
Exit 0 means arithmetic checks passed, 1 means failed/unavailable core evidence,
and 2 means invalid CLI/input/output. None of these grants scientific acceptance.
Only an explicitly requested, previously nonexistent report file is written.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import sys

SCHEMA = "pion_background_chain_audit/v1"
COMPONENTS = ("pi_n", "pi_sidis", "pi_delta")
SOURCES = ("prompt", "rand", "dummy", "dummy_rand")
MOMENTS = ("signed", "absolute", "sumsq", "count")
# Relative to operand/support scales, including cancellation terms; rounding
# only, no absolute yield floor and no statistical/physics acceptance threshold.
ROUNDING_REL = 1e-10
UNKNOWN = "NOT VERIFIED"
LIMITS = [
    "Diagnostic only; scientific agreement and production acceptance are NOT VERIFIED.",
    "Conditional fixed-transfer-weight pull, not a formal fit p-value.",
    "Conditional uncertainty omits pion-fit covariance, SIMC-template uncertainty, "
    "fitted-event correlation, and transfer-weight uncertainty, even if retained.",
    "Producer fingerprints establish internal consistency, not independent run provenance.",
    "Method A/B have no authority here; no correction, fallback, or production decision.",
]


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def digest(value):
    data = (json.dumps(value, sort_keys=True, allow_nan=False, separators=(",", ":")) + "\n").encode()
    return hashlib.sha256(data).hexdigest()


def unknown(reason):
    return {"status": UNKNOWN, "reason": reason}


class EvidenceError(ValueError):
    pass


def require(condition, reason):
    if not condition:
        raise EvidenceError(reason)


def close(a, b, *scales):
    if not all(finite(x) for x in (a, b, *scales)):
        return False
    scale = max(abs(a), abs(b), *(abs(x) for x in scales))
    return abs(a - b) <= ROUNDING_REL * scale


def vector(value, size, name, nonnegative=False):
    require(isinstance(value, list) and len(value) == size, name + ": missing/wrong bin count")
    require(all(finite(x) for x in value), name + ": nonfinite/non-numeric bin")
    require(not nonnegative or all(x >= 0 for x in value), name + ": negative support")
    return value


def histogram(value, name, edges=None):
    require(isinstance(value, dict) and value.get("status") == "available", name + ": unavailable")
    axis = value.get("edges")
    require(isinstance(axis, list) and len(axis) >= 2 and all(finite(x) for x in axis)
            and all(a < b for a, b in zip(axis, axis[1:])), name + ": missing/invalid axis")
    require(edges is None or axis == edges, name + ": axis mismatch")
    vector(value.get("values"), len(axis) + 1, name + ".values")
    errors = vector(value.get("errors"), len(axis) + 1, name + ".errors", True)
    # Exact scalar encoding of the producer's ROOT-free ownership fingerprint.
    scalars = [float(len(axis) - 1), float(axis[0]), float(axis[-1]), *map(float, axis)]
    for v, e in zip(value["values"], errors):
        scalars.extend((float(v), float(e)))
    fp = hashlib.sha256(b"".join(struct.pack("!d", x) for x in scalars)).hexdigest()
    require(value.get("axis_signature") == digest(axis), name + ": axis signature mismatch/missing")
    require(value.get("content_error_fingerprint") == fp, name + ": fingerprint mismatch/missing")
    return value


def compare(left, right, name, scales=None):
    require(len(left) == len(right), name + ": bin count mismatch")
    for i, (a, b) in enumerate(zip(left, right)):
        require(close(a, b, *(scales[i] if scales else ())), name + ": mismatch at bin " + str(i))


def check_flags(tree, name):
    """Honor retained invariant failures, but never treat a saved PASS as proof."""
    if isinstance(tree, dict):
        require(tree.get("status") != "FAIL", name + ": retained invariant FAIL")
        for key, value in tree.items():
            check_flags(value, name + "." + key)
    elif isinstance(tree, list):
        for i, value in enumerate(tree):
            check_flags(value, name + "." + str(i))


def moments(value, edges, name):
    require(isinstance(value, dict) and value.get("status") == "available", name + ": unavailable")
    require(value.get("edges") == edges, name + ": axis mismatch")
    require(not value.get("invalid_records"), name + ": invalid source records")
    size = len(edges) + 1
    combined = value.get("combined", {})
    for key in MOMENTS:
        vector(combined.get(key), size, name + ".combined." + key, key != "signed")
    source = value.get("sources", {})
    require(set(source) == set(SOURCES), name + ": missing signed source classes")
    for s in SOURCES:
        for key in MOMENTS:
            vector(source[s].get(key), size, name + "." + s + "." + key, key != "signed")
    for key in MOMENTS:
        compare(combined[key], [math.fsum(source[s][key][i] for s in SOURCES) for i in range(size)],
                name + ": source sum " + key,
                [(combined["absolute"][i],) if key == "signed" else () for i in range(size)])
    for cell in [combined, *source.values()]:
        for i in range(size):
            a, c, v, s = (cell[k][i] for k in ("absolute", "count", "sumsq", "signed"))
            require(c == int(c), name + ": noninteger event count")
            require(a >= abs(s) or close(a, abs(s)), name + ": absolute support below signed magnitude")
            require(c != 0 or (a == s == v == 0), name + ": nonzero support with no events")
            require(v != 0 or (a == 0 and s == 0), name + ": nonzero weight with zero Sumw2")
            if v > 0:
                require(a * a <= c * v or close(a * a, c * v), name + ": impossible source moments")
    return combined, source


def moment_closure(h, cell, name):
    compare(h["values"], cell["signed"], name + ": signed content",
            [(x,) for x in cell["absolute"]])
    compare([e * e for e in h["errors"]], cell["sumsq"], name + ": error squared")
    if isinstance(h.get("stored_sumw2"), list):
        compare(vector(h["stored_sumw2"], len(h["values"]), name + ".stored_sumw2", True),
                cell["sumsq"], name + ": stored Sumw2")


def metadata_regions(fit, edges):
    """Use retained coordinate windows only; do not assume physical MM ranges.

    Windows select native bins by inclusive center, without fractional rebinning.
    Each window is a named diagnostic region, not a reconstructed fit selection.
    """
    regions = {}
    def add(name, value):
        if (isinstance(value, list) and len(value) == 2 and all(finite(x) for x in value)
                and value[0] < value[1]):
            regions[name] = {"window": list(value), "bins": [i for i in range(1, len(edges))
                if value[0] <= (edges[i-1] + edges[i]) / 2 <= value[1]],
                "selection": "native bin centers inclusive; no fractional bins or fit-mask inference"}
        elif isinstance(value, list):
            for i, item in enumerate(value):
                add(name + "." + str(i), item)
        elif isinstance(value, dict):
            for key, item in value.items():
                add(name + "." + key, item)
    def walk(node, path):
        if not isinstance(node, dict):
            return
        if finite(node.get("fit_min")) and finite(node.get("fit_max")):
            add(path + ".fit_range", [node["fit_min"], node["fit_max"]])
        for key, value in node.items():
            if key in {"anchor_windows", "stage_amplitude_windows", "fit_window", "lambda_gauge_window"}:
                add(path + "." + key, value)
            elif isinstance(value, dict):
                walk(value, path + "." + key)
    walk(fit, "fit")
    return regions


def fit_summary(fit):
    result = {"producer_validation_flags": {k: copy.deepcopy(fit.get(k, unknown("not retained")))
              for k in ("fit_validation_pion", "fit_validation_kaon", "fit_status_pion", "fit_status_kaon")},
              "flag_interpretation": "existing producer flags only; accepted is not a new result"}
    fields = ("success", "status", "message", "chi2", "chi2_data", "ndf", "chi2_ndf", "fit_p_value",
              "uncertainties", "amplitude_uncertainties", "bound_hit_flags", "covariance_matrix",
              "coefficient_covariance_matrix", "fit_metrics", "validation", "fit_bin_indices")
    for side in ("pion", "kaon"):
        diag = fit.get("diagnostics", {}).get(side, {})
        if not isinstance(diag, dict):
            diag = {}
        blocks = {"parent": diag, "joint_refinement": diag.get("joint_refinement", {}),
                  "protected_fit": diag.get("pi_delta_signal_protected_fit", {})}
        result[side] = {}
        for name, block in blocks.items():
            block = block if isinstance(block, dict) else {}
            data = {k: copy.deepcopy(block[k]) for k in fields if k in block}
            flags = block.get("bound_hit_flags")
            data["boundary_hit_count"] = (sum(v is True for v in flags.values())
                if isinstance(flags, dict) and flags else unknown("boundary flags not retained"))
            data["amplitude_uncertainty_status"] = ("retained; not propagated" if
                block.get("uncertainties") or block.get("amplitude_uncertainties") else UNKNOWN)
            result[side][name] = data
    return result


def support(cell, bins):
    sums = {k: math.fsum(cell[k][i] for i in bins) for k in MOMENTS}
    s, a, v = (sums[k] for k in ("signed", "absolute", "sumsq"))
    valid = v > 0 and all(cell["sumsq"][i] > 0 or cell["absolute"][i] == 0 for i in bins)
    sums.update(signed_effective_support=(s / math.sqrt(v)) ** 2 if valid else None,
                absolute_effective_support=(a / math.sqrt(v)) ** 2 if valid else None,
                effective_support_status="available" if valid else UNKNOWN,
                cancellation_fraction=1 - abs(s) / a if a > 0 else None)
    return sums


def region_report(bins, qualified, obs, model, kaon, applied, weights, cell, sources):
    q = [i for i in bins if i in qualified]
    def discrepancy(indices):
        d = [applied[i] - kaon[i] for i in indices]
        return {"D_R": math.fsum(d), "L_R": math.fsum(abs(x) for x in d)}
    direct = discrepancy(bins)
    actual = discrepancy(q)
    predicted = [weights[i] * (obs[i] - model[i]) for i in q]
    variance_ok = bool(q) and all(cell["sumsq"][i] > 0 or cell["absolute"][i] == 0 for i in q)
    sigma = math.hypot(*(abs(weights[i]) * math.sqrt(cell["sumsq"][i]) for i in q)) if variance_ok else None
    norm = math.fsum(abs(kaon[i]) for i in bins)
    return {"bins": bins, "qualified_bins": q, "excluded_bins": [i for i in bins if i not in qualified],
        "sums": {name: math.fsum(values[i] for i in bins) for name, values in
                 (("P_obs", obs), ("P_model", model), ("K_model", kaon), ("C_applied", applied))},
        "source_support": support(cell, bins),
        "by_source": {s: support(c, bins) for s, c in sources.items()},
        "direct_all_normal_bins": direct,
        "qualified": {"status": "available" if q else UNKNOWN, "direct": actual,
            "reconstructed": {"D_R": math.fsum(predicted), "L_R": math.fsum(abs(x) for x in predicted)},
            "conditional_fixed_weight_sigma": sigma,
            "conditional_pull": actual["D_R"] / sigma if sigma is not None and sigma > 0 else None,
            "uncertainty_status": "conditional only" if sigma is not None and sigma > 0 else UNKNOWN},
        "sum_abs_K_model": norm,
        "normalized_signed_discrepancy": direct["D_R"] / norm if norm > 0 else None,
        "normalized_absolute_discrepancy": direct["L_R"] / norm if norm > 0 else None,
        "normalization_status": "available" if norm > 0 else UNKNOWN}


def analyze_parent(parent, cache_coordinate, requested_regions):
    result = {"identity": copy.deepcopy(parent.get("identity")), "arithmetic_status": "FAIL",
              "scientific_agreement": UNKNOWN, "regions": {},
              "uncertainty_budget": {name: unknown("not independently calibrated or propagated") for name in
                  ("pion_fit_covariance", "SIMC_template_uncertainty", "fitted_event_correlation", "transfer_weight_uncertainty")}}
    try:
        require(parent.get("status") == "available", "parent unavailable")
        ident = parent.get("identity", {})
        require(isinstance(ident, dict) and ident.get("pion_parent_id"), "missing parent identity")
        require(parent.get("identity_closure") == "PASS", "parent identity closure missing/failed")
        require(ident.get("t_bin_index") == parent.get("t_bin_index"), "t identity mismatch")
        te = ident.get("t_edges")
        require(isinstance(te, list) and len(te) == 2 and all(finite(x) for x in te) and te[0] < te[1],
                "missing/invalid canonical t edges")
        require(isinstance(cache_coordinate, str) and bool(cache_coordinate)
                and ident.get("coordinate_fingerprint") == cache_coordinate, "coordinate identity mismatch/missing")
        control = histogram(parent.get("control"), "control")
        edges, obs = control["edges"], control["values"]
        require(control["content_error_fingerprint"] == ident.get("pion_control_input_fingerprint"),
                "control input identity mismatch/missing")
        size = len(obs)
        cell, sources = moments(parent.get("source_spectra"), edges, "source_spectra")
        require(parent["source_spectra"].get("selection") == "nommcuts", "control source selection mismatch/missing")
        moment_closure(control, cell, "control/source")
        fit_input_closure = parent.get("fit_input_closure")
        require(isinstance(fit_input_closure, dict) and fit_input_closure.get("status") == "PASS",
                "fit-input closure missing/invalid/failed")
        for key in ("source_closure", "cache_control_closure", "fit_input_closure"):
            check_flags(parent.get(key), key)
        transfer = parent.get("transfer", {})
        hs = transfer.get("histograms", {})
        p, k, wh = [histogram(hs.get(key), key, edges) for key in
            ("H_pion_control_model", "H_kaon_pion_model", "H_pion_weight_vs_MM")]
        model, kaon, weights = p["values"], k["values"], wh["values"]
        compare(vector(transfer.get("weights"), size, "transfer.weights"), weights, "weight array/histogram")
        floor, lo, hi = (transfer.get(key) for key in ("denominator_floor", "clip_min", "clip_max"))
        require(finite(floor) and floor >= 0, "denominator policy missing/invalid")
        def bound_known(key, value):
            if key not in transfer or isinstance(value, dict) and value.get("status") == "unavailable":
                return False
            require(value is None or finite(value), "clipping policy invalid")
            return True
        lo_known, hi_known = bound_known("clip_min", lo), bound_known("clip_max", hi)
        require(not (lo_known and hi_known and lo is not None and hi is not None) or lo <= hi,
                "clipping policy invalid")
        retained_bins = {b["bin"]: b for b in transfer.get("bins", []) if isinstance(b, dict) and "bin" in b}
        fit = parent.get("fit", {})
        require(isinstance(fit, dict), "fit evidence unavailable")
        shapes = [histogram(fit.get("histograms", {}).get("H_simc_shape_" + c), "shape." + c, edges)
                  for c in COMPONENTS]
        for side, values in (("pion", model), ("kaon", kaon)):
            amps = transfer.get("selected_amplitudes", {}).get(side, {})
            require(all(finite(amps.get(c)) for c in COMPONENTS), "selected amplitudes missing/invalid")
            terms = [[amps[c] * h["values"][i] for c, h in zip(COMPONENTS, shapes)] for i in range(size)]
            compare(values, [math.fsum(t) for t in terms], side + " component-model algebra", terms)
        check_flags(transfer.get("closure"), "transfer.closure")
        final = parent.get("final", {})
        applied_h = histogram(final.get("histograms", {}).get("H_pion_subtraction_template_MM_nosub"),
                              "final applied", edges)
        applied = applied_h["values"]
        reference = histogram(final.get("weight_reference"), "final weight reference", edges)
        compare(reference["values"], model, "final/transfer reference")
        compare(vector(final.get("weights"), size, "final.weights"), weights,
                "final/transfer weights (fallback/proposal cannot substitute)")
        final_cell, _ = moments(final.get("nommcuts_weighted_sources"), edges, "final weighted sources")
        require(final["nommcuts_weighted_sources"].get("selection") == "nommcuts", "final selection mismatch")
        source_fingerprint = parent["source_spectra"].get("records_fingerprint")
        final_fingerprint = final["nommcuts_weighted_sources"].get("records_fingerprint")
        require(isinstance(source_fingerprint, str) and bool(source_fingerprint.strip())
                and isinstance(final_fingerprint, str) and bool(final_fingerprint.strip())
                and source_fingerprint == final_fingerprint, "source record provenance mismatch/missing/invalid")
        moment_closure(applied_h, final_cell, "final reproduction")
        reproduction = final.get("nommcuts_event_reproduction", {})
        require(reproduction.get("content", {}).get("status") == "PASS"
                and reproduction.get("error_squared", {}).get("status") == "PASS",
                "stored independent final event reproduction missing/failed")
        check_flags(reproduction, "final event reproduction")
        check_flags(final.get("nommcuts_after_closure"), "final after closure")
        qualified, bins = [], []
        for i in range(size):
            flow = i == 0 or i == size - 1
            supported = model[i] > floor
            raw = kaon[i] / model[i] if supported else None
            finite_ratio = raw is not None and finite(raw)
            policy_raw = raw if finite_ratio else 0.0
            expected = policy_raw
            if not flow:
                if lo_known and lo is not None:
                    expected = max(lo, expected)
                if hi_known and hi is not None:
                    expected = min(hi, expected)
            else:
                expected = 0.0  # producer leaves transfer flow entries zero
            if not flow:
                require(not lo_known or lo is None or weights[i] >= lo or close(weights[i], lo),
                        "known clipping bound violated at bin " + str(i))
                require(not hi_known or hi is None or weights[i] <= hi or close(weights[i], hi),
                        "known clipping bound violated at bin " + str(i))
            policy_known = flow or lo_known and hi_known
            if policy_known:
                require(close(weights[i], expected), "stored numerical policy mismatch at bin " + str(i))
                clipped = not flow and expected != policy_raw
            else:
                # snapshot(None) loses the distinction between disabled and
                # missing bounds. Do not invent an unbounded policy. An unchanged
                # raw ratio proves no clipping effect; changed weights require
                # the producer's per-bin clipping observation and stay qualified
                # out. Full configured-policy verification remains unavailable.
                clipped = not close(weights[i], policy_raw)
                if clipped:
                    saved = retained_bins.get(i, {})
                    require(saved.get("configured_clipping") is True and saved.get("stored_weight_policy_match") is True,
                            "changed ratio without retained clipping evidence at bin " + str(i))
                else:
                    require(close(weights[i], policy_raw),
                            "unsupported weight without numerical policy at bin " + str(i))
            q = not flow and supported and finite_ratio and not clipped
            if q:
                qualified.append(i)
                require(close(weights[i] * model[i], kaon[i]), "model ratio algebra at bin " + str(i))
                direct = applied[i] - kaon[i]
                predicted = weights[i] * (obs[i] - model[i])
                require(close(direct, predicted, applied[i], kaon[i], weights[i] * obs[i], weights[i] * model[i]),
                        "final residual identity at bin " + str(i))
            factor = 0.0 if i == 0 else weights[i]
            for s in SOURCES:
                actual = final["nommcuts_weighted_sources"]["sources"][s]
                base = sources[s]
                for key, f in (("signed", factor), ("absolute", abs(factor)), ("sumsq", factor * factor), ("count", 1)):
                    require(close(actual[key][i], base[key][i] * f, base["absolute"][i] * abs(f)
                                  if key == "signed" else 0), "signed-event transfer " + s + "." + key + " bin " + str(i))
            bins.append({"bin": i, "flow": "underflow" if i == 0 else "overflow" if flow else None,
                "P_model": model[i], "K_model": kaon[i], "weight": weights[i],
                "denominator_supported": supported, "denominator_floor_touch": model[i] <= floor,
                "unsupported_positive_numerator": not supported and kaon[i] > floor,
                "clipped": clipped, "ratio_nonfinite": supported and not finite_ratio,
                "identity_qualified": q})
        regions = {"normal": {"bins": list(range(1, size - 1)), "selection": "all native normal bins"}}
        windows = metadata_regions(fit, edges)
        regions.update(windows)
        for name, indices in requested_regions.items():
            require(name not in regions, "duplicate region name: " + name)
            require(indices and len(set(indices)) == len(indices) and
                    all(isinstance(i, int) and not isinstance(i, bool) and 1 <= i < size - 1 for i in indices),
                    "requested region outside native normal bins: " + name)
            regions[name] = {"bins": sorted(indices), "selection": "explicit native normal-bin indices"}
        region_results = {}
        for name, spec in regions.items():
            region_results[name] = {**spec, **region_report(spec["bins"], qualified, obs, model, kaon, applied, weights, cell, sources),
                "exclusions": {key: [i for i in spec["bins"] if bins[i][key]] for key in
                               ("denominator_floor_touch", "unsupported_positive_numerator", "clipped", "ratio_nonfinite")}}
        result.update(arithmetic_status="PASS", bins=bins,
            numerical_policy={"denominator_floor": floor,
                "clip_min": copy.deepcopy(lo), "clip_max": copy.deepcopy(hi),
                "configured_clipping_status": "available" if lo_known and hi_known else UNKNOWN,
                "limitation": "Unavailable serialized bounds are not assumed disabled; ratio changes require retained clipping evidence."},
            counts={key: sum(bool(b[key]) for b in bins if b["flow"] is None) for key in
                    ("denominator_floor_touch", "unsupported_positive_numerator", "clipped", "ratio_nonfinite", "identity_qualified")},
            flows={label: {"P_obs": obs[i], "P_model": model[i], "K_model": kaon[i], "C_applied": applied[i],
                          "support": support(cell, [i]), "included_in_regions": False}
                   for label, i in (("underflow", 0), ("overflow", size - 1))},
            metadata_window_status="available" if windows else UNKNOWN,
            fit=fit_summary(fit), simc_provenance=copy.deepcopy(parent.get("simc_provenance", unknown("not retained"))),
            simc_raw_occupancy=unknown("normalized shape bins are not raw SIMC event counts"),
            regions=region_results)
    except (EvidenceError, KeyError, TypeError, ValueError, OverflowError, AttributeError) as exc:
        result.update(first_failed_invariant=str(exc), interpretation="STOP; no region interpretation",
                      core_evidence=unknown(str(exc)))
    return result


def analyze(document, regions=None):
    """Return new scalar diagnostics; input dictionaries are never mutated."""
    report = {"schema_version": "pion_transfer_validity/v1", "arithmetic_status": "FAIL",
              "scientific_agreement": UNKNOWN, "limits": list(LIMITS), "parents": [],
              "definitions": {"signed_effective_support": "(sum signed source weights)^2 / sumw2; cancellation-sensitive precision diagnostic, not an independent event count",
                  "absolute_effective_support": "(sum absolute source weights)^2 / sumw2",
                  "cancellation_fraction": "1 - abs(sum signed source weights) / sum absolute source weights",
                  "D_R": "sum(C_applied - K_model)", "L_R": "sum(abs(C_applied - K_model)); no signed cancellation",
                  "rounding_relative_tolerance": ROUNDING_REL}}
    try:
        require(isinstance(document, dict), "JSON top level must be an object")
        require(document.get("schema_version") == SCHEMA, "unsupported audit schema")
        audit = document.get("audit", document)
        require(isinstance(audit, dict) and audit.get("schema_version") == SCHEMA, "unsupported nested audit schema")
        parents = audit.get("parents")
        require(isinstance(parents, list) and parents and all(isinstance(p, dict) for p in parents), "no compatible parents")
        require(audit.get("production_objects_mutated") is False and audit.get("diagnostic_only") is True,
                "producer read-only diagnostic provenance missing/contradictory")
        identities = [p.get("identity", p) for p in parents]
        require(audit.get("collection_identity") == digest(identities), "collection identity mismatch/missing")
        for key in ("pion_parent_id", "t_bin_index"):
            tokens = [digest(p.get("identity", {}).get(key)) for p in parents]
            require(len(set(tokens)) == len(tokens), "duplicate parent " + key)
        head = audit.get("observed_source_HEAD")
        source_known = isinstance(head, str) and re.fullmatch(r"[0-9a-fA-F]{40}", head) is not None
        report["producer_provenance"] = {"observed_source_HEAD": head,
            "source_identity_status": "retained; not independently attested" if source_known else UNKNOWN,
            "producer_status": document.get("status", UNKNOWN),
            "capture_completeness": copy.deepcopy(audit.get("capture_completeness", unknown("not retained"))),
            "scope": copy.deepcopy(audit.get("scope")), "collection_identity": audit["collection_identity"]}
        coordinate = audit.get("cache_identity", {}).get("coordinate_fingerprint")
        report["parents"] = [analyze_parent(p, coordinate, regions or {}) for p in parents]
        # One failed parent stops interpretation of the entire collection.
        failed = [p for p in report["parents"] if p["arithmetic_status"] != "PASS"]
        if failed:
            report["first_failed_invariant"] = failed[0]["first_failed_invariant"]
            for p in report["parents"]:
                p["regions"] = {}
                p["interpretation"] = "STOP; collection contains a failed invariant"
        else:
            report["arithmetic_status"] = "PASS"
    except (EvidenceError, KeyError, TypeError, ValueError, OverflowError, AttributeError) as exc:
        report["first_failed_invariant"] = str(exc)
        report["core_evidence"] = unknown(str(exc))
        report["interpretation"] = "STOP; no collection interpretation"
    return report


def strict_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key: " + key)
        result[key] = value
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit_json", type=Path)
    parser.add_argument("--report", type=Path, help="create a new report; never overwrite existing files")
    parser.add_argument("--bins", action="append", default=[], metavar="NAME=1,2", help="native normal-bin region")
    args = parser.parse_args(argv)
    try:
        regions = {}
        for spec in args.bins:
            name, indices = spec.split("=", 1)
            if not name or name in regions:
                raise ValueError("empty/duplicate region name")
            regions[name] = [int(i) for i in indices.split(",")]
        raw = args.audit_json.read_bytes()
        document = json.loads(raw, object_pairs_hook=strict_object,
                              parse_constant=lambda s: (_ for _ in ()).throw(ValueError("nonfinite JSON: " + s)))
        report = analyze(document, regions)
        report["input_sha256"] = hashlib.sha256(raw).hexdigest()
        encoded = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
        if args.report is not None:
            if args.report.resolve() == args.audit_json.resolve():
                raise ValueError("report must not alias input")
            with args.report.open("x", encoding="utf-8", newline="\n") as handle:
                handle.write(encoded)
        print(encoded, end="")
        return 0 if report["arithmetic_status"] == "PASS" else 1
    except (OSError, ValueError) as exc:
        print("pion-transfer audit: " + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
