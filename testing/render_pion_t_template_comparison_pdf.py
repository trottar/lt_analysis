#!/usr/bin/env python3
"""Replay one frozen six-page A/B JSON as ROOT presentation only; never refit."""
from __future__ import annotations

import argparse
from array import array
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "src/cuts"), str(REPO / "src/utility")]
import pion_t_template_comparison as comparison


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def validate_histogram(value):
    require(isinstance(value, dict) and value.get("status") == "available", "histogram unavailable")
    edges, values, errors = (value.get(k) for k in ("edges", "values", "errors"))
    require(all(isinstance(a, list) for a in (edges, values, errors)), "histogram arrays missing")
    require(len(edges) >= 2 and len(values) == len(errors) == len(edges)+1, "histogram flow geometry")
    require(all(isinstance(x, (float, int)) and not isinstance(x, bool) and math.isfinite(x)
                for a in (edges, values, errors) for x in a), "nonfinite/non-numeric histogram")
    require(all(a < b for a, b in zip(edges, edges[1:])) and min(errors) >= 0, "invalid axis/errors")
    require(value.get("bin_layout") == "underflow, normal bins 1..N, overflow; no width weighting",
            "histogram bin layout differs")
    require(value.get("axis_signature") == comparison.digest(edges), "histogram axis signature differs")
    numbers = [float(len(edges)-1), float(edges[0]), float(edges[-1]), *map(float, edges)]
    numbers += [float(x) for pair in zip(values, errors) for x in pair]
    fingerprint = hashlib.sha256(b"".join(struct.pack("!d", x) for x in numbers)).hexdigest()
    require(value.get("content_error_fingerprint") == fingerprint, "histogram content/error fingerprint differs")
    sumw2 = value.get("stored_sumw2")
    require(isinstance(sumw2, list) and len(sumw2) == len(values), "explicit histogram Sumw2 missing")
    require(all(isinstance(x, (float, int)) and not isinstance(x, bool) and math.isfinite(x) and x >= 0
                for x in sumw2), "invalid stored Sumw2")
    require(all(math.isclose(math.sqrt(v), e, rel_tol=1e-12, abs_tol=1e-15)
                for v, e in zip(sumw2, errors)), "stored errors/Sumw2 differ")
    return value


def validate_document(document, source_commit):
    require(re.fullmatch(r"[0-9a-f]{40}", source_commit or "") is not None, "exact source commit required")
    require(document.get("schema_version") == comparison.SCHEMA and
            document.get("source_HEAD") == source_commit and document.get("status") == "available" and
            document.get("diagnostic_only") is True and document.get("production_applied") is False and
            document.get("scientifically_accepted") is False and not document.get("source_identity_rejection_reasons"),
            "scientific source/schema/diagnostic identity differs")
    require(document.get("phi_setting") in ("Center", "Left") and document.get("setting") == {
        "Q2": "4p4", "W": "2p74", "EPSSET": "low", "ParticleType": "kaon"}, "setting scope differs")
    parents = document.get("parents")
    require(isinstance(parents, list) and [p.get("t_bin_index") for p in parents] == [0, 1, 2],
            "three ordered canonical parents required")
    fits = []
    for parent in parents:
        require(parent.get("status") == "available" and parent.get("authoritative_parent_unchanged") is True,
                "parent unavailable or mutated")
        pair = []
        for key in ("baseline_A", "candidate_B"):
            fit = parent.get(key)
            require(isinstance(fit, dict), "stored A/B fit missing")
            selected = {}
            for side, data_key in (("pion", "H_pion_control_input"), ("kaon", "H_kaon_nosub_input")):
                selected[data_key] = validate_histogram(fit.get(data_key))
                for suffix in ("", "_refined"):
                    total = "H_{}_fit_total{}".format(side, suffix)
                    selected[total] = validate_histogram(fit.get(total))
                    for component in (*comparison.COMPONENTS, "k_lambda", "k_sigma0"):
                        name = "H_{}_fit_{}_scaled{}".format(side, component, suffix)
                        value = fit.get(name)
                        # Serializer represents absent optional ROOT components as unavailable.
                        if value is None or (isinstance(value, dict) and value.get("status") == "unavailable"):
                            continue
                        selected[name] = validate_histogram(value)
                    require(all(value["edges"] == selected[data_key]["edges"] for name, value in selected.items()
                                if name.startswith("H_"+side+"_")), "display histogram axes differ")
            pair.append(selected)
        for data_key in ("H_pion_control_input", "H_kaon_nosub_input"):
            require(pair[0][data_key] == pair[1][data_key], "A/B frozen data differs")
        fits.append((parent, *pair))
    return fits


def restore_histogram(value, name, root):
    validate_histogram(value)
    hist = root.TH1D(name, "", len(value["edges"])-1, array("d", value["edges"]))
    hist.SetDirectory(0)
    hist.Sumw2()
    for i, (content, error, variance) in enumerate(zip(value["values"], value["errors"], value["stored_sumw2"])):
        hist.SetBinContent(i, content)
        hist.SetBinError(i, error)
        hist.GetSumw2().SetAt(variance, i)
    retained = comparison.histogram(hist)
    require(retained.get("content_error_fingerprint") == value["content_error_fingerprint"] and
            retained.get("stored_sumw2") == value["stored_sumw2"], "ROOT histogram restoration differs")
    return hist


def replay(input_json, output_pdf, report_json, source_commit):
    input_json, output_pdf, report_json = map(Path, (input_json, output_pdf, report_json))
    require(input_json.is_file(), "input JSON missing")
    require(output_pdf.suffix == ".pdf" and report_json.suffix == ".json", "output extensions differ")
    outputs = [p.resolve() for p in (output_pdf, report_json)]
    require(len(set(outputs)) == 2 and input_json.resolve() not in outputs, "input/output collision")
    require(all(not p.is_relative_to(REPO.resolve()) for p in outputs), "outputs must stay outside source checkout")
    require(all(p.parent.is_dir() and not os.path.lexists(p) for p in (output_pdf, report_json)),
            "outputs must be new files in existing directories")
    raw = input_json.read_bytes()
    document = json.loads(raw)
    selected = validate_document(document, source_commit)
    import ROOT  # ROOT is used only after validating every stored scientific input.
    ROOT.gROOT.SetBatch(True)
    rows = [(parent, *[{key: restore_histogram(value, "stored_{}_{}_{}".format(pi, fi, key), ROOT)
                       for key, value in fit.items()} for fi, fit in enumerate((a, b))])
            for pi, (parent, a, b) in enumerate(selected)]
    comparison.render_comparison_pdf(output_pdf, rows, "{} Q4p4W2p74 lowe | scientific source {} | stored JSON replay | diagnostic only".format(
        document["phi_setting"], source_commit[:12]))
    require(input_json.read_bytes() == raw, "scientific input changed during rendering")
    require(output_pdf.is_file() and output_pdf.read_bytes().startswith(b"%PDF-"), "ROOT PDF missing/invalid header")
    renderer = Path(comparison.__file__).resolve()
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()
    report = {"schema_version": "pion_t_template_comparison_presentation_replay/v1",
        "scientific_input": {"path": str(input_json.resolve()), "sha256": hashlib.sha256(raw).hexdigest(),
                             "source_HEAD": source_commit},
        "renderer": {"observed_git_HEAD": head,
                     "files": {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in (renderer, Path(__file__).resolve())}},
        "output_pdf": {"path": str(output_pdf.resolve()), "sha256": hashlib.sha256(output_pdf.read_bytes()).hexdigest()},
        "expected_pages": 6, "diagnostic_only": True, "refitted": False,
        "production_applied": False, "pdf_legibility_verified": False}
    with report_json.open("x", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-json", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path, required=True)
    parser.add_argument("--report-json", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    args = parser.parse_args(argv)
    report = replay(args.input_json, args.output_pdf, args.report_json, args.source_commit)
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
