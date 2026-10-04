"""Produce aggregate current-lineage diagnostic JSON and an eight-page review PDF."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import textwrap

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))
import pion_hgcer_method_a_current_lineage_diagnostic as diagnostic
import pion_hgcer_method_a_parallel_full_procedure as lineage


def output_paths(outdir, kinematic, phi="Left", epsilon="lowe"):
    if kinematic != diagnostic.KINEMATIC or phi != "Left" or epsilon != "lowe":
        raise diagnostic.CurrentLineageDiagnosticError("diagnostic_scope_unsupported")
    root = Path(outdir).resolve()
    if not root.is_dir():
        raise diagnostic.CurrentLineageDiagnosticError("existing_outdir_required")
    return root / (diagnostic.STEM + ".json"), root / (diagnostic.STEM + ".pdf")


def _git_value(args):
    return subprocess.check_output(["git", *args], cwd=REPO_ROOT, text=True).strip()


def load_diagnostic(outdir, kinematic):
    paths = lineage.accepted_f6_3_artifact_paths(outdir, kinematic)
    f1, f3, f4, hashes = lineage.load_accepted_f6_3_authority(paths)
    result = diagnostic.build_diagnostic(
        f1, f3, f4, f1_input_file_hashes={k: v for k, v in hashes.items() if k not in ("f3", "f4")},
        f3_input_file_sha256=hashes["f3"], f4_input_file_sha256=hashes["f4"], kinematic=kinematic)
    return diagnostic.build_artifact(result, provenance={
        "git_head": _git_value(["rev-parse", "HEAD"]),
        "git_status_short": _git_value(["status", "--short"]),
        "input_paths": paths,
        "runtime_claim": "detached diagnostic only; no production promotion",
    })


def write_review_pdf(path, artifact):
    """Render persisted aggregates only; no transient factors or ROOT required."""
    diagnostic.guard_aggregate_persistence(artifact)
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    result = artifact["diagnostic"]
    a = result["measurement_a"]["parents"]
    b = result["measurement_b"]["parents"]
    if result["available"] is not True or [p["canonical_t_index"] for p in a] != [0, 1, 2] or [p["canonical_t_index"] for p in b] != [0, 1, 2]:
        raise diagnostic.CurrentLineageDiagnosticError("pdf_parent_inventory_invalid")

    def text_page(title, lines):
        figure = plt.figure(figsize=(11.7, 8.3))
        figure.suptitle(title, fontsize=14)
        wrapped = "\n".join(textwrap.fill(str(line), 110, break_long_words=False) for line in lines)
        figure.text(.04, .91, wrapped, va="top", fontsize=8, family="monospace", linespacing=1.5)
        return figure

    def save(figure):
        pdf.savefig(figure)
        plt.close(figure)

    def plot_bins(axis, table, fields, labels):
        rows = table["bins"]
        for field, label in zip(fields, labels):
            axis.step(range(len(rows)), [row[field] for row in rows], where="mid", label=label)
        axis.set_xlabel("fixed display bin (first/last = under/overflow)")
        axis.legend(fontsize=7)
        axis.grid(alpha=.2)

    metadata = {"Title": "Left/lowe current-lineage Method-A diagnostic", "Author": "KaonLT",
                "CreationDate": None, "ModDate": None}
    with PdfPages(path, metadata=metadata) as pdf:
        authority = result["input_authority"]
        lines = ["Q4p4W2p74 / Left / lowe; t0,t1,t2; primary target = canonical index t1",
                 "Detached, aggregate-only, non-authoritative. F.4 payload reproduced exactly.",
                 "NPE=0: operational kaon PID category; NPE>0: pion tree; NPE>2: physical pion control.",
                 "No HGC-free true-pion tag in consumed artifacts; no direct absolute mis-ID calibration.",
                 "Measurement A: observed weak-positive relative topology evidence; proxy validity NOT VERIFIED.",
                 "Measurement B: F.4 signed-normalization sensitivity; N_abs descriptive only.",
                 "Current candidate validation source: " + result["candidate_validation_source_head"],
                 "F.3 file SHA-256: " + result["input_file_hashes"]["f3"],
                 "F.3 map fingerprint: " + result["f3_map_fingerprint"],
                 "F.4 correction fingerprint: " + result["f4_exact_reproduction"]["fingerprint"]]
        lines.extend(f"observed {key}: {value}" for key, value in sorted(authority["observed"].items()) if not isinstance(value, (dict, list)))
        lines.append("F.1 observed file hashes:")
        lines.extend(f"  {item['setting_id']}: {item['source_file_sha256']}" for item in authority_f1(artifact))
        save(text_page("1 — Authority, scope and limitations", lines))

        figure, axes = plt.subplots(3, 2, figsize=(11.7, 8.3))
        figure.suptitle("2 — Measurement A: weak-positive and control response")
        quantiles = ("min", "p01", "p10", "p25", "p50", "p75", "p90", "p99", "max")
        for index, parent in enumerate(a):
            bands = parent["positive_npe_bands"]
            axes[index, 0].bar(range(4), [p["count"] for p in bands])
            axes[index, 0].set_xticks(range(4), [p["definition"] for p in bands], fontsize=7)
            axes[index, 0].set_ylabel(f"t{index} count")
            for name, p in sorted(parent["classes"].items()):
                s = p["raw_relative_response"]
                if s["count"]:
                    axes[index, 1].plot(range(9), [s[k] for k in quantiles], marker=".", label=f"{name} n={p['count']}")
            for p in bands:
                s = p["raw_relative_response"]
                if s["count"]:
                    axes[index, 1].plot(range(9), [s[k] for k in quantiles], linestyle=":", linewidth=.8,
                                        label=p["definition"])
            axes[index, 1].set_xticks(range(9), quantiles, rotation=30, fontsize=6)
            axes[index, 1].set_ylabel("unweighted r quantiles")
            axes[index, 1].legend(fontsize=6)
        figure.tight_layout(rect=(0, 0, 1, .94)); save(figure)

        figure, axes = plt.subplots(3, 3, figsize=(11.7, 8.3))
        figure.suptitle("3 — Measurement A: hgcer3 topology (counts, not mis-ID calibration)")
        for index, parent in enumerate(a):
            for col, feature in enumerate(diagnostic.FEATURES):
                for label, p in sorted(parent["classes"].items()):
                    rows = p["coordinate_bins"][feature]["bins"]
                    axes[index, col].step(range(len(rows)), [r["count"] for r in rows], where="mid", label=label)
                edges = parent["classes"]["weak_positive"]["coordinate_bins"][feature]["edges"]
                axes[index, col].set_title(f"t{index} {feature} [{edges[0]:g}, {edges[-1]:g}]", fontsize=8)
                axes[index, col].set_xlabel("fixed bin incl. under/overflow", fontsize=7)
                axes[index, col].legend(fontsize=6)
        support_lines = [f"t{i}: weak in/OOD={p['classes']['weak_positive']['in_support_count']}/{p['classes']['weak_positive']['ood_count']}; matched/train-only/physical-only="
                         f"{p['training_physical_population_shift']['matched_count']}/{p['training_physical_population_shift']['training_control_only_count']}/{p['training_physical_population_shift']['physical_control_only_count']}"
                         for i, p in enumerate(a)]
        figure.text(.03, .02, "; ".join(support_lines), fontsize=6)
        figure.tight_layout(rect=(0, .07, 1, .94)); save(figure)

        lines = ["Measurement B only: N_abs is a descriptive comparator, never an applied normalization.",
                 "parent       B_signed       A_abs       |B|/A       U_signed       V_abs        N_signed       N_abs",
                 "----------------------------------------------------------------------------------------------------"]
        for p in b:
            n = p["normalization"]
            lines.append("t{}  ".format(p["canonical_t_index"]) + "  ".join(f"{n[k]:.7g}" for k in ("B_signed", "A_abs", "cancellation_ratio", "U_signed", "V_abs", "N_signed", "N_abs")))
            lines.append(f"  N_signed/N_abs={n['normalization_ratio']:.8g}; relative difference={n['relative_normalization_difference']:.8g}; accepted closure residual={p['accepted_closure_residual']:.8g}")
        save(text_page("4 — Measurement B: signed/absolute parent normalization", lines))

        target = b[1]
        lines = ["t1 authoritative statistical/source classes (not true-species labels)",
                 "source count B_signed A_abs U_signed V_abs signed_fraction abs_fraction cancellation"]
        for row in target["source_decomposition"]:
            lines.append(row["source_label"] + " " + str(row["count"]) + " " + " ".join(
                "unavailable" if row[k] is None else f"{row[k]:.8g}" for k in ("B_signed", "A_abs", "U_signed", "V_abs", "signed_fraction", "absolute_support_fraction", "cancellation_ratio")))
        for name in ("r", "C_signed", "w0", "b", "abs_b"):
            lines.append(name + " tails: " + ", ".join(f"{k}={target['distributions'][name][k]:.7g}" for k in quantiles))
        lines.append("Application support: " + json.dumps(target["support"], sort_keys=True))
        save(text_page("5 — Measurement B: t1 source decomposition and tails", lines))

        figure, axes = plt.subplots(2, 1, figsize=(11.7, 8.3))
        figure.suptitle("6 — Measurement B: t1 child/MM signed sensitivity")
        for field, label in (("baseline_signed_sum", "baseline"), ("adjusted_signed_sum", "accepted F.4 adjusted")):
            axes[0].step(range(len(target["phi"]["bins"])), [r[field] for r in target["phi"]["bins"]], where="mid", label=label)
        axes[0].set_xticks(range(len(target["phi"]["bins"])), [f"{lo:g} to {hi:g}" for lo, hi in zip(target["phi"]["edges"], target["phi"]["edges"][1:])], fontsize=7)
        axes[0].set_xlabel("frozen canonical phi child [degrees]"); axes[0].legend()
        plot_bins(axes[1], target["coordinate_bins"]["analysis_MM"], ("baseline_signed_sum", "adjusted_signed_sum"), ("baseline", "accepted F.4 adjusted"))
        axes[1].set_title("analysis_MM; fixed 0.8–1.4 bins plus under/overflow")
        figure.tight_layout(rect=(0, 0, 1, .94)); save(figure)

        figure, axes = plt.subplots(3, 1, figsize=(11.7, 8.3))
        figure.suptitle("7 — Measurement B: t1 coordinate sensitivity")
        for axis, feature in zip(axes, diagnostic.FEATURES):
            plot_bins(axis, target["coordinate_bins"][feature], ("baseline_signed_sum", "adjusted_signed_sum"), ("baseline", "accepted F.4 adjusted"))
            edges = target["coordinate_bins"][feature]["edges"]
            axis.set_title(f"{feature}; fixed range {edges[0]:g} to {edges[-1]:g}", fontsize=9)
        figure.tight_layout(rect=(0, 0, 1, .94)); save(figure)

        save(text_page("8 — Interpretation boundary", [
            "Measured: positive HGC response topology, support and population matching; accepted signed/absolute support and normalization contrasts.",
            "INFERENCE: these measurements may inform the weak-positive proxy question and signed-normalization sensitivity hypothesis.",
            "NOT VERIFIED: direct pion-to-kaon mis-ID calibration, proxy validity, detector hardware cause or production correctness.",
            "NPE=0 is kaon-selected. No HGC-independent pion truth tag is invented.",
            "N_abs is descriptive only. No alternative factor, child normalization, yield or cross section is constructed.",
            "No overall validity score, winner, tuned threshold or promotion recommendation.",
            "Method A remains detached. Baseline w0 and production remain authoritative.",
            "No absolute-SIMC amplitude interpretation. Await scientific review of returned measurements.",
        ]))


def authority_f1(artifact):
    """Render only aggregate input identities, never event identities."""
    record = artifact["diagnostic"]["input_file_hashes"]["f1"]
    return [{"setting_id": key, "source_file_sha256": value} for key, value in sorted(record.items())]


def build_argument_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--kinematic", choices=(diagnostic.KINEMATIC,), required=True)
    parser.add_argument("--phi", choices=("Left",), default="Left")
    parser.add_argument("--epsilon", choices=("lowe",), default="lowe")
    return parser


def main(argv=None):
    args = build_argument_parser().parse_args(argv)
    created = []
    try:
        json_path, pdf_path = output_paths(args.outdir, args.kinematic, args.phi, args.epsilon)
        if json_path.exists() or pdf_path.exists():
            raise diagnostic.CurrentLineageDiagnosticError("output_already_exists")
        artifact = load_diagnostic(args.outdir, args.kinematic)
        diagnostic.guard_aggregate_persistence(artifact)
        with tempfile.TemporaryDirectory(dir=args.outdir, prefix=".current-lineage-") as temp:
            rendered = Path(temp) / pdf_path.name
            write_review_pdf(rendered, artifact)
            with json_path.open("xb") as handle:
                created.append(json_path)
                handle.write((json.dumps(artifact, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8"))
            with pdf_path.open("xb") as handle:
                created.append(pdf_path)
                with rendered.open("rb") as source:
                    shutil.copyfileobj(source, handle)
    except (ValueError, OSError, lineage.MethodAParallelFullProcedureError, subprocess.SubprocessError) as exc:
        for path in created:
            path.unlink()
        print("error: " + str(exc), file=sys.stderr)
        return 1
    print("wrote " + str(json_path)); print("wrote " + str(pdf_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
