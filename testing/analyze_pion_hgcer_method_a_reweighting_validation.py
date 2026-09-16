"""Run detached Phase F.6.1 Method-A reweighting validation from frozen inputs."""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Mapping, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))
import pion_hgcer_method_a_reweighting_validation as validation  # noqa: E402


CANONICAL_SETTINGS = (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"), ("Center", "highe"), ("Right", "highe"))


def _f1_filename(phi: str, kinematic: str, epsilon: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-contract_{}_{}.json".format(phi, kinematic, epsilon)


def _f3_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-map.json".format(kinematic)


def _f4_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.json".format(kinematic)


def _f5_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-tphi-propagation.json".format(kinematic)


def _pdf_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-reweighting-validation.pdf".format(kinematic)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path, label: str) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError("{}_json_invalid".format(label)) from exc
    if not isinstance(payload, dict):
        raise ValueError("{}_json_invalid".format(label))
    return payload


def _git_value(arguments: Sequence[str]) -> str | None:
    try:
        completed = subprocess.run(["git", *arguments], cwd=REPO_ROOT, text=True, capture_output=True, check=False)
    except OSError:
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


def load_inputs(outdir: Path, kinematic: str) -> tuple[list[dict[str, object]], dict[str, str], dict[str, str], dict[str, object], str, str, dict[str, object], str, str, dict[str, object], str, str]:
    if not outdir.is_dir():
        raise ValueError("f6_1_outdir_invalid")
    artifacts: list[dict[str, object]] = []; hashes: dict[str, str] = {}; paths: dict[str, str] = {}; seen: set[Path] = set()
    for phi, epsilon in CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon); path = outdir / _f1_filename(phi, kinematic, epsilon)
        if path in seen:
            raise ValueError("f6_1_f1_input_duplicate")
        seen.add(path)
        if not path.is_file():
            raise ValueError("f6_1_f1_input_missing:{}".format(setting_id))
        payload = _read_json(path, "f6_1_f1_input"); setting = payload.get("setting")
        if not isinstance(setting, Mapping) or setting.get("phi_setting") != phi or setting.get("epsilon_filename_token") != epsilon or setting.get("kinematic_token") != kinematic or setting.get("particle_type") != "kaon":
            raise ValueError("f6_1_f1_input_identity_invalid:{}".format(setting_id))
        artifacts.append(payload); hashes[setting_id] = _sha256(path); paths[setting_id] = str(path.resolve())
    f3_path = outdir / _f3_filename(kinematic); f4_path = outdir / _f4_filename(kinematic); f5_path = outdir / _f5_filename(kinematic)
    if f3_path in seen or f4_path in seen or f5_path in seen or len({f3_path, f4_path, f5_path}) != 3:
        raise ValueError("f6_1_upstream_input_path_collides")
    if not f3_path.is_file() or not f4_path.is_file() or not f5_path.is_file():
        raise ValueError("f6_1_upstream_input_missing")
    return artifacts, hashes, paths, _read_json(f3_path, "f6_1_f3_input"), _sha256(f3_path), str(f3_path.resolve()), _read_json(f4_path, "f6_1_f4_input"), _sha256(f4_path), str(f4_path.resolve()), _read_json(f5_path, "f6_1_f5_input"), _sha256(f5_path), str(f5_path.resolve())


def _parents(result: Mapping[str, object], setting_id: str) -> list[Mapping[str, object]]:
    parents = result.get("parents")
    values = [row for row in parents if isinstance(row, Mapping) and row.get("setting_id") == setting_id] if isinstance(parents, list) else []
    if len(values) != 3:
        raise ValueError("f6_1_pdf_parent_inventory_invalid")
    return sorted(values, key=lambda row: int(row["canonical_t_index"]))


def _stairs(axis: object, values: Sequence[float], edges: Sequence[float], label: str, color: str) -> None:
    axis.stairs(np.asarray(values, dtype=float), np.asarray(edges, dtype=float), label=label, color=color)  # type: ignore[attr-defined]


def _overlay_shape(axis: object, comparison: Mapping[str, object], title: str) -> None:
    edges = comparison["edges"]
    _stairs(axis, comparison["low_response_unit_area"], edges, "Low response", "tab:green")
    _stairs(axis, comparison["baseline_prompt_control_unit_area"], edges, "Baseline control", "tab:blue")
    _stairs(axis, comparison["method_a_prompt_control_unit_area"], edges, "Method-A control", "tab:orange")
    axis.set_title("{}; H0={:.3g}, HA={:.3g}".format(title, float(comparison["hellinger_baseline"]), float(comparison["hellinger_method_a"])))  # type: ignore[attr-defined]
    axis.set_ylabel("Unit-area shape")  # type: ignore[attr-defined]


def _text_page(pdf: object, title: str, lines: Sequence[str]) -> None:
    import matplotlib.pyplot as plt
    figure = plt.figure(figsize=(11, 8.5)); figure.suptitle(title, fontsize=15, fontweight="bold")
    figure.text(0.03, 0.94, "\n".join(lines), va="top", family="monospace", fontsize=7.5)
    pdf.savefig(figure); plt.close(figure)


def write_review_pdf(path: Path, artifact: Mapping[str, object]) -> None:
    import matplotlib
    matplotlib.use("Agg", force=True)
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    result = artifact.get("validation")
    if not isinstance(result, Mapping) or result.get("available") is not True:
        raise ValueError("f6_1_pdf_validation_unavailable")
    with PdfPages(path) as pdf:
        authority = result["runtime_authority"]
        _text_page(pdf, "F6.1 — detached Method-A reweighting validation", [
            "Purpose: compare the accepted F.4 reweighting with the observed prompt low-HGCer response.",
            "Low-response, baseline-control, and Method-A-control curves are independently unit-area normalized.",
            "Signed full-background curves are not normalized; negative weighted contents are physical subtraction values.",
            "F4 and F5 are recomputed from the frozen inputs before this artifact is available.",
            "No production weights, ROOT objects, yield, cross section, child normalization, smoothing, or Method B input exists.",
            "", "F5 scientific fingerprint = {}".format(authority["accepted"]["f5_propagation_fingerprint"]),
            "F5.2 accepted source SHA-256 = {}".format(authority["accepted"]["f5_source_file_sha256"]),
            "Shape improvement is descriptive only; no automatic gate is applied.",
        ])
        for variables, title in ((_MODEL_VARIABLES(), "model variables"), (("SHMS_xptar", "SHMS_yptar", "phi"), "independent variables"), (("analysis_MM", "Q2", "W"), "missing mass and kinematics")):
            for phi, epsilon in CANONICAL_SETTINGS:
                setting_id = "{}-{}".format(phi, epsilon); rows = _parents(result, setting_id)
                figure, axes = plt.subplots(3, 3, figsize=(14, 10)); figure.suptitle("F6.1 {} — {}".format(setting_id, title), fontsize=14, fontweight="bold")
                for row_index, parent in enumerate(rows):
                    comparisons = parent["prompt_shape_comparisons"]
                    for column, variable in enumerate(variables):
                        _overlay_shape(axes[row_index, column], comparisons[variable], variable)
                        axes[row_index, column].set_xlabel("phi [degrees]" if variable == "phi" else variable)
                        if row_index == 0 and column == 0: axes[row_index, column].legend(fontsize=6)
                figure.tight_layout(rect=(0, 0, 1, 0.95)); pdf.savefig(figure); plt.close(figure)
        for phi, epsilon in CANONICAL_SETTINGS:
            setting_id = "{}-{}".format(phi, epsilon); rows = _parents(result, setting_id)
            figure, axes = plt.subplots(3, 3, figsize=(14, 10)); figure.suptitle("F6.1 {} — HGCer detector position".format(setting_id), fontsize=14, fontweight="bold")
            for row_index, parent in enumerate(rows):
                xy = parent["hgcer_xy"]
                for column, name in enumerate(("low_response_unit_area", "baseline_prompt_control_unit_area", "method_a_prompt_control_unit_area")):
                    mesh = axes[row_index, column].pcolormesh(np.asarray(xy["x_edges"], dtype=float), np.asarray(xy["y_edges"], dtype=float), np.asarray(xy[name], dtype=float).T, shading="flat")
                    figure.colorbar(mesh, ax=axes[row_index, column]); axes[row_index, column].set_title(name.replace("_", " "))
                    axes[row_index, column].set_xlabel("P_hgcer_xAtCer"); axes[row_index, column].set_ylabel("P_hgcer_yAtCer")
            figure.tight_layout(rect=(0, 0, 1, 0.95)); pdf.savefig(figure); plt.close(figure)
        for phi, epsilon in CANONICAL_SETTINGS:
            setting_id = "{}-{}".format(phi, epsilon); rows = _parents(result, setting_id)
            figure, axes = plt.subplots(3, 2, figsize=(12, 10)); figure.suptitle("F6.1 {} — signed physical pion background".format(setting_id), fontsize=14, fontweight="bold")
            for row_index, parent in enumerate(rows):
                for column, variable in enumerate(("analysis_MM", "phi")):
                    signed = parent["signed_background"][variable]
                    _stairs(axes[row_index, column], signed["baseline_signed_contents"], signed["edges"], "Baseline", "tab:blue")
                    _stairs(axes[row_index, column], signed["method_a_signed_contents"], signed["edges"], "Method A", "tab:orange")
                    _stairs(axes[row_index, column], signed["signed_delta_contents"], signed["edges"], "Difference", "tab:purple")
                    axes[row_index, column].axhline(0.0, color="black", linewidth=0.6); axes[row_index, column].set_title(variable); axes[row_index, column].set_ylabel("Signed weighted events")
                    if row_index == 0 and column == 0: axes[row_index, column].legend(fontsize=6)
            figure.tight_layout(rect=(0, 0, 1, 0.95)); pdf.savefig(figure); plt.close(figure)
        summary = result["f5_continuity"]
        _text_page(pdf, "F6.1 — fifteen-parent authority and closure", [
            "All 15 canonical parents completed prompt identity/parity, F4, and F5 gates.",
            "F5 continuity: {} settings, {} parents, {} explicit canonical cells.".format(summary["setting_count"], summary["parent_count"], summary["cell_count"]),
            "F4 exact reproduction = {}".format(result["f4_reproduction"]["exact_payload_match"]),
            "F5 exact reproduction = {}".format(result["f5_reproduction"]["exact_payload_match"]),
            "", "setting / t       prompt controls   low reference   in support   OOD",
            *["{:<18} {:>8} {:>14} {:>12} {:>5}".format("{} / {}".format(parent["setting_id"], int(parent["canonical_t_index"]) + 1), parent["population_counts"]["prompt_physical_control"], parent["population_counts"]["low_response_prompt_training"], parent["support"]["in_support_count"], parent["support"]["ood_count"]) for parent in result["parents"]],
            "", "Detached shadow-only artifact. Human review is required; no shape metric is an automatic acceptance decision.",
        ])


def _MODEL_VARIABLES() -> tuple[str, str, str]:
    return ("SHMS_delta", "P_hgcer_xAtCer", "P_hgcer_yAtCer")


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path); parser.add_argument("--kinematic", required=True)
    parser.add_argument("--output-json", required=True, type=Path); parser.add_argument("--output-pdf", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None, *, accepted_runtime_authority_by_kinematic: Mapping[str, object] | None = None, accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None = None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    expected_json = (arguments.outdir / validation.pion_hgcer_method_a_reweighting_validation_filename(arguments.kinematic)).resolve(); expected_pdf = (arguments.outdir / _pdf_filename(arguments.kinematic)).resolve()
    output_json, output_pdf = arguments.output_json.resolve(), arguments.output_pdf.resolve()
    if output_json != expected_json or output_pdf != expected_pdf:
        print("error: f6_1_output_path_not_deterministic", file=sys.stderr); return 2
    if output_json == output_pdf:
        print("error: f6_1_output_paths_collide", file=sys.stderr); return 2
    if output_json.exists() or output_pdf.exists():
        print("error: f6_1_output_path_already_exists", file=sys.stderr); return 2
    temporary_json: Path | None = None; temporary_pdf: Path | None = None
    promoted_json = False; promoted_pdf = False
    try:
        f1, hashes, f1_paths, f3, f3_sha, f3_path, f4, f4_sha, f4_path, f5, f5_sha, f5_path = load_inputs(arguments.outdir, arguments.kinematic)
        artifact = validation.build_pion_hgcer_method_a_reweighting_validation_artifact(f1, f3, f4, f5, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, f4_input_file_sha256=f4_sha, f5_input_file_sha256=f5_sha, input_paths={"f1": f1_paths, "f3": f3_path, "f4": f4_path, "f5": f5_path}, generated_at_utc=_datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z"), git_head=_git_value(("rev-parse", "HEAD")), git_status_short=_git_value(("status", "--short")), accepted_runtime_authority_by_kinematic=accepted_runtime_authority_by_kinematic, accepted_f4_runtime_authority_by_kinematic=accepted_f4_runtime_authority_by_kinematic, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
        json_fd, json_name = tempfile.mkstemp(prefix=".f6_1_", suffix=".json", dir=arguments.outdir); os.close(json_fd); temporary_json = Path(json_name)
        pdf_fd, pdf_name = tempfile.mkstemp(prefix=".f6_1_", suffix=".pdf", dir=arguments.outdir); os.close(pdf_fd); temporary_pdf = Path(pdf_name)
        validation.write_pion_hgcer_method_a_reweighting_validation_json(temporary_json, artifact); write_review_pdf(temporary_pdf, artifact)
        # Promote the PDF first: a late filesystem failure can never leave a
        # final-path JSON claiming success without its required review PDF.
        os.replace(temporary_pdf, output_pdf); temporary_pdf = None; promoted_pdf = True
        os.replace(temporary_json, output_json); temporary_json = None; promoted_json = True
    except (validation.MethodAReweightingValidationError, OSError, ValueError) as exc:
        for path in (temporary_json, temporary_pdf):
            if path is not None and path.exists(): path.unlink()
        if promoted_json and output_json.exists(): output_json.unlink()
        if promoted_pdf and output_pdf.exists(): output_pdf.unlink()
        print("error: {}".format(exc), file=sys.stderr); return 1
    print("wrote {}".format(output_json)); print("wrote {}".format(output_pdf)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
