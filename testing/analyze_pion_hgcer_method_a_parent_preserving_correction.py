"""Run detached Phase F.4 parent-preserving Method-A correction construction."""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Mapping, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))
import pion_hgcer_method_a_parent_preserving_correction as correction  # noqa: E402

CANONICAL_SETTINGS = (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"), ("Center", "highe"), ("Right", "highe"))


def _f1_filename(phi: str, kinematic: str, epsilon: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-contract_{}_{}.json".format(phi, kinematic, epsilon)


def _f3_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-map.json".format(kinematic)


def _pdf_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.pdf".format(kinematic)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_value(arguments: Sequence[str]) -> str | None:
    try:
        completed = subprocess.run(["git", *arguments], cwd=REPO_ROOT, text=True, capture_output=True, check=False)
    except OSError:
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


def load_inputs(outdir: Path, kinematic: str) -> tuple[list[dict[str, object]], dict[str, str], dict[str, str], dict[str, object], str, str]:
    if not outdir.is_dir():
        raise ValueError("f4_outdir_invalid")
    artifacts: list[dict[str, object]] = []; hashes: dict[str, str] = {}; paths: dict[str, str] = {}; seen: set[Path] = set()
    for phi, epsilon in CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon); path = outdir / _f1_filename(phi, kinematic, epsilon)
        if path in seen: raise ValueError("f4_f1_input_duplicate")
        seen.add(path)
        if not path.is_file(): raise ValueError("f4_f1_input_missing:{}".format(setting_id))
        try: payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError) as exc: raise ValueError("f4_f1_input_json_invalid:{}".format(setting_id)) from exc
        if not isinstance(payload, dict) or not isinstance(payload.get("setting"), Mapping): raise ValueError("f4_f1_input_json_invalid:{}".format(setting_id))
        setting = payload["setting"]
        if setting.get("phi_setting") != phi or setting.get("epsilon_filename_token") != epsilon or setting.get("kinematic_token") != kinematic or setting.get("particle_type") != "kaon": raise ValueError("f4_f1_input_identity_invalid:{}".format(setting_id))
        artifacts.append(payload); hashes[setting_id] = _sha256(path); paths[setting_id] = str(path.resolve())
    f3_path = outdir / _f3_filename(kinematic)
    if f3_path in seen: raise ValueError("f4_f3_input_path_collides")
    if not f3_path.is_file(): raise ValueError("f4_f3_input_missing")
    try: f3 = json.loads(f3_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc: raise ValueError("f4_f3_input_json_invalid") from exc
    if not isinstance(f3, dict): raise ValueError("f4_f3_input_json_invalid")
    return artifacts, hashes, paths, f3, _sha256(f3_path), str(f3_path.resolve())


def _text_page(lines: Sequence[str], title: str) -> object:
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    figure = plt.figure(figsize=(11, 8.5)); figure.suptitle(title, fontsize=15, fontweight="bold")
    figure.text(0.03, 0.94, "\n".join(lines), va="top", family="monospace", fontsize=7.5)
    return figure


def write_review_pdf(path: Path, artifact: Mapping[str, object], review_data: Sequence[Mapping[str, object]]) -> None:
    """Render aggregate-only seven-page F.4 review PDF."""
    import matplotlib
    matplotlib.use("Agg", force=True)
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    result = artifact.get("correction")
    if not isinstance(result, Mapping) or result.get("available") is not True: raise ValueError("f4_pdf_correction_unavailable")
    parents = result.get("parents")
    if not isinstance(parents, list) or len(parents) != 15: raise ValueError("f4_pdf_parent_inventory_invalid")
    review_by_parent: dict[tuple[str, int], Mapping[str, object]] = {}
    for review in review_data:
        if not isinstance(review, Mapping): raise ValueError("f4_pdf_review_data_invalid")
        key = (str(review.get("setting_id")), int(review.get("canonical_t_index")))
        if key in review_by_parent: raise ValueError("f4_pdf_review_data_duplicate")
        factors = np.asarray(review.get("correction_factors"), dtype=float)
        if factors.ndim != 1 or factors.size == 0 or not np.all(np.isfinite(factors)) or np.any(factors <= 0.0): raise ValueError("f4_pdf_review_factors_invalid")
        review_by_parent[key] = review
    if len(review_by_parent) != 15: raise ValueError("f4_pdf_review_data_inventory_invalid")
    with PdfPages(path) as pdf:
        lines = ["Detached parent-preserving Method-A-only correction", "basis = hgcer3", "A = exp(beta dot z) in F3 support; A = 1 outside support", "N = sum(b A) / sum(b); C = A / N; sum(b C) = sum(b)", "F3 source SHA-256 = {}".format(result["f3_source_file_sha256"]), "F3 map fingerprint = {}".format(result["f3_map_fingerprint"]), "15/15 parents required; no source or child normalization.", "No Method B, template, yield, probability, or production application."]
        pdf.savefig(_text_page(lines, "F4.1 — authority and frozen correction formula")); plt.close()
        for phi, epsilon in CANONICAL_SETTINGS:
            setting = [row for row in parents if isinstance(row, Mapping) and isinstance(row.get("setting"), Mapping) and row["setting"].get("phi_setting") == phi and row["setting"].get("epsilon_filename_token") == epsilon]
            if len(setting) != 3: raise ValueError("f4_pdf_setting_inventory_invalid")
            figure, axes = plt.subplots(3, 2, figsize=(12, 10)); figure.suptitle("F4 — {} {} detached parent corrections".format(phi, epsilon), fontsize=14, fontweight="bold")
            for index, row in enumerate(sorted(setting, key=lambda value: int(value["canonical_t_index"]))):
                summary = row["correction_factor_summary"]; shape = row["raw_shape_factor_summary"]; left, right = axes[index]
                review = review_by_parent.get((str(row["setting_id"]), int(row["canonical_t_index"])))
                if review is None: raise ValueError("f4_pdf_review_data_missing")
                correction_values = np.asarray(review["correction_factors"], dtype=float)
                child_text = "; ".join("phi {}: {:+.3g}".format(child["phi_index"], child["signed_delta"]) for child in row["canonical_phi_diagnostics"])
                left.axis("off")
                left.text(0, 1, "t{t}  B={B:.8g}\nU={U:.8g}\nN={N:.8g}\nin/OOD={inside}/{ood} ({frac:.4f})\nclosure={res:.3g}  PASS={passed}\nA p01/p50/p99={a1:.4g}/{a5:.4g}/{a9:.4g}\nC p01/p50/p99={c1:.4g}/{c5:.4g}/{c9:.4g}\nchild signed deltas: {child}".format(t=row["canonical_t_index"], B=row["baseline_parent_sum"], U=row["raw_shape_parent_sum"], N=row["parent_normalization"], inside=row["in_support_count"], ood=row["ood_count"], frac=row["ood_fraction"], res=row["closure_residual"], passed=row["closure_passed"], a1=shape["p01"], a5=shape["p50"], a9=shape["p99"], c1=summary["p01"], c5=summary["p50"], c9=summary["p99"], child=child_text), va="top", family="monospace", fontsize=8)
                log_correction_values = np.sort(np.log10(correction_values))
                right.step(log_correction_values, np.arange(1, log_correction_values.size + 1, dtype=float) / log_correction_values.size, where="post", color="tab:blue")
                for name, color in (("p01", "tab:green"), ("p50", "tab:orange"), ("p99", "tab:red"), ("max", "tab:purple")):
                    right.axvline(np.log10(float(summary[name])), color=color, linewidth=0.8, linestyle="--", label=name)
                right.set_xlabel("log10(C)"); right.set_ylabel("ECDF"); right.set_title("unclipped transient correction-factor ECDF")
                right.legend(fontsize=6, loc="lower right")
                child = row["canonical_phi_diagnostics"]
                right.text(0.01, 0.02, "child aggregates only; no child normalization\n{} child bins".format(len(child)), transform=right.transAxes, fontsize=7)
            figure.tight_layout(rect=(0, 0.03, 1, 0.95)); pdf.savefig(figure); plt.close(figure)
        closure = ["parent                   B             U             N        adjusted      residual     OOD     C p01/p50/p99/max", "---------------------------------------------------------------------------------------------------------------"]
        for row in parents:
            setting = row["setting"]; summary = row["correction_factor_summary"]
            closure.append("{:<24} {:>10.4g} {:>10.4g} {:>10.4g} {:>10.4g} {:>10.2g} {:>7.4f} {:>7.3g}/{:>7.3g}/{:>7.3g}/{:>7.3g}".format("{}-{} t{}".format(setting["phi_setting"], setting["epsilon_filename_token"], row["canonical_t_index"]), row["baseline_parent_sum"], row["raw_shape_parent_sum"], row["parent_normalization"], row["adjusted_parent_sum"], row["closure_residual"], row["ood_fraction"], summary["p01"], summary["p50"], summary["p99"], summary["max"]))
        closure.extend(("", "All rows are detached aggregate diagnostics. No production correction or template exists."))
        pdf.savefig(_text_page(closure, "F4.7 — fifteen-parent closure")); plt.close()


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path); parser.add_argument("--kinematic", required=True)
    parser.add_argument("--output-json", type=Path); parser.add_argument("--output-pdf", type=Path)
    return parser


def main(argv: Sequence[str] | None = None, *, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    expected_json = (arguments.outdir / correction.pion_hgcer_method_a_parent_preserving_correction_filename(arguments.kinematic)).resolve(); expected_pdf = (arguments.outdir / _pdf_filename(arguments.kinematic)).resolve()
    output_json = expected_json if arguments.output_json is None else arguments.output_json.resolve(); output_pdf = expected_pdf if arguments.output_pdf is None else arguments.output_pdf.resolve()
    if output_json != expected_json or output_pdf != expected_pdf: print("error: f4_output_path_not_deterministic", file=sys.stderr); return 2
    if output_json == output_pdf: print("error: f4_output_paths_collide", file=sys.stderr); return 2
    if output_json.exists() or output_pdf.exists(): print("error: f4_output_path_already_exists", file=sys.stderr); return 2
    try:
        f1, hashes, f1_paths, f3, f3_sha, f3_path = load_inputs(arguments.outdir, arguments.kinematic)
        artifact, review_data = correction.build_pion_hgcer_method_a_parent_preserving_correction_artifact_with_review_data(f1, f3, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, input_paths={"f1": f1_paths, "f3": f3_path}, generated_at_utc=_datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z"), git_head=_git_value(("rev-parse", "HEAD")), git_status_short=_git_value(("status", "--short")), accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
        if artifact["correction"]["available"] is not True: raise ValueError("f4_correction_unavailable")  # type: ignore[index]
        correction.write_pion_hgcer_method_a_parent_preserving_correction_json(output_json, artifact); write_review_pdf(output_pdf, artifact, review_data)
    except (correction.MethodAParentPreservingCorrectionError, OSError, ValueError) as exc:
        print("error: {}".format(exc), file=sys.stderr); return 1
    print("wrote {}".format(output_json)); print("wrote {}".format(output_pdf)); return 0


if __name__ == "__main__": raise SystemExit(main())
