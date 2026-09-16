"""Run detached Phase F.6.2 refinement validation from frozen artifacts."""

from __future__ import annotations

import argparse
import datetime as _datetime
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from typing import Mapping, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))
import pion_hgcer_method_a_acceptance_refinement_validation as validation  # noqa: E402
import pion_hgcer_method_a_reweighting_validation as f6_1  # noqa: E402


KINEMATIC = "Q4p4W2p74"
CANONICAL_SETTINGS = (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"), ("Center", "highe"), ("Right", "highe"))


def _f1_filename(phi: str, kinematic: str, epsilon: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-contract_{}_{}.json".format(phi, kinematic, epsilon)


def _f3_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-map.json".format(kinematic)


def _f4_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.json".format(kinematic)


def _f5_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-tphi-propagation.json".format(kinematic)


def _f6_1_filename(kinematic: str) -> str:
    return f6_1.pion_hgcer_method_a_reweighting_validation_filename(kinematic)


def _pdf_filename(kinematic: str) -> str:
    if kinematic != KINEMATIC:
        raise ValueError("f6_2_kinematic_unsupported")
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.pdf".format(kinematic)


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


def load_inputs(outdir: Path, kinematic: str) -> tuple[list[dict[str, object]], dict[str, str], dict[str, str], dict[str, object], str, str, dict[str, object], str, str, dict[str, object], str, str, dict[str, object], str, str]:
    if kinematic != KINEMATIC:
        raise ValueError("f6_2_kinematic_unsupported")
    if not outdir.is_dir():
        raise ValueError("f6_2_outdir_invalid")
    f1_artifacts: list[dict[str, object]] = []; f1_hashes: dict[str, str] = {}; f1_paths: dict[str, str] = {}
    for phi, epsilon in CANONICAL_SETTINGS:
        path = outdir / _f1_filename(phi, kinematic, epsilon)
        if not path.is_file():
            raise ValueError("f6_2_f1_input_missing:{}".format(phi))
        artifact = _read_json(path, "f6_2_f1_{}".format(phi)); setting = artifact.get("setting")
        if not isinstance(setting, Mapping) or setting.get("phi_setting") != phi or setting.get("epsilon_filename_token") != epsilon:
            raise ValueError("f6_2_f1_identity_mismatch:{}".format(phi))
        setting_id = "{}-{}".format(phi, epsilon)
        f1_artifacts.append(artifact); f1_hashes[setting_id] = _sha256(path); f1_paths[setting_id] = str(path.resolve())
    f3_path, f4_path, f5_path, f6_path = (outdir / _f3_filename(kinematic), outdir / _f4_filename(kinematic), outdir / _f5_filename(kinematic), outdir / _f6_1_filename(kinematic))
    for path, label in ((f3_path, "f3"), (f4_path, "f4"), (f5_path, "f5"), (f6_path, "f6_1")):
        if not path.is_file():
            raise ValueError("f6_2_{}_input_missing".format(label))
    return (f1_artifacts, f1_hashes, f1_paths, _read_json(f3_path, "f6_2_f3"), _sha256(f3_path), str(f3_path.resolve()), _read_json(f4_path, "f6_2_f4"), _sha256(f4_path), str(f4_path.resolve()), _read_json(f5_path, "f6_2_f5"), _sha256(f5_path), str(f5_path.resolve()), _read_json(f6_path, "f6_2_f6_1"), _sha256(f6_path), str(f6_path.resolve()))


def _text_page(pdf: object, title: str, lines: Sequence[str]) -> None:
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(11, 8.5)); figure.suptitle(title, fontsize=15, fontweight="bold")
    figure.text(0.06, 0.92, "\n".join(lines), va="top", ha="left", family="monospace", fontsize=8.7, wrap=True)
    pdf.savefig(figure); plt.close(figure)


def _metric_text(value: object) -> str:
    if not isinstance(value, Mapping):
        return "unavailable: malformed"
    if value.get("available") is True:
        if "value" in value:
            return "{:.4g}".format(float(value["value"]))
        if "ci_low" in value and "ci_high" in value:
            return "CI [{:.4g}, {:.4g}]".format(float(value["ci_low"]), float(value["ci_high"]))
        return "available"
    return "unavailable: {}".format(value.get("reason"))


def _overlay(axis: object, payload: Mapping[str, object], title: str) -> None:
    edges = np.asarray(payload.get("edges"), dtype=float)
    if edges.ndim != 1 or edges.size < 2:
        axis.set_title(title + " (invalid edges)"); return
    for name, color in (("L", "#1f77b4"), ("B", "#333333"), ("A", "#7f3c8d")):
        population = payload.get(name)
        if not isinstance(population, Mapping):
            continue
        contents = np.asarray(population.get("unit_area"), dtype=float)
        if contents.shape == (edges.size - 1,):
            axis.stairs(contents, edges, label=name, color=color, linewidth=1.5)
    axis.set_title(title, fontsize=9); axis.legend(fontsize=7); axis.set_ylabel("unit area")


def _matrix(axis: object, payload: Mapping[str, object], population_name: str, title: str, vmax: float) -> None:
    x_edges, y_edges = np.asarray(payload.get("x_edges"), dtype=float), np.asarray(payload.get("y_edges"), dtype=float)
    population = payload.get(population_name)
    contents = np.asarray(population.get("unit_area"), dtype=float) if isinstance(population, Mapping) else np.asarray([])
    if x_edges.ndim == y_edges.ndim == 1 and contents.shape == (x_edges.size - 1, y_edges.size - 1):
        axis.pcolormesh(x_edges, y_edges, contents.T, shading="auto", cmap="viridis", vmin=0.0, vmax=vmax)
    axis.set_title(title, fontsize=8); axis.set_xlabel("{}".format(payload.get("x_variable", "x"))); axis.set_ylabel("{}".format(payload.get("y_variable", "y")))


def _parent_overview_page(pdf: object, parent: Mapping[str, object]) -> None:
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(14, 8.5)); axis.axis("off")
    title = "F.6.2 {} t{} — canonical-phi overview".format(parent.get("setting_id"), parent.get("canonical_t_index")); figure.suptitle(title, fontsize=14, fontweight="bold")
    headings = ("phi", "N low/control/full", "N_eff B/A", "OOD ctrl/full", "DeltaH MM", "R MM", "kappa MM", "MM×xptar", "MM×yptar", "DeltaP K", "RV")
    rows: list[list[str]] = []
    for child in parent.get("children", []):
        counts = child.get("population_counts", {}); support = child.get("support", {}); neff = child.get("effective_sample_size", {}); shapes = child.get("one_dimensional", {}); mm = shapes.get("analysis_MM", {}) if isinstance(shapes, Mapping) else {}; metrics = mm.get("metrics", {}) if isinstance(mm, Mapping) else {}; joint = child.get("joint_distributions", {}); signed = child.get("signed_background", {}); window = signed.get("kaon_window", {}) if isinstance(signed, Mapping) else {}; variance = signed.get("variance_proxy", {}) if isinstance(signed, Mapping) else {}
        x_kappa = ((joint.get("analysis_MM__SHMS_xptar", {}) if isinstance(joint, Mapping) else {}).get("metrics", {}) if isinstance(joint, Mapping) else {}).get("kappa", {})
        y_kappa = ((joint.get("analysis_MM__SHMS_yptar", {}) if isinstance(joint, Mapping) else {}).get("metrics", {}) if isinstance(joint, Mapping) else {}).get("kappa", {})
        rows.append(["{} [{:.0f},{:.0f})".format(child.get("phi_index"), float(child.get("phi_low", 0.0)), float(child.get("phi_high", 0.0))), "{}/{}/{}".format(counts.get("N_low"), counts.get("N_control"), counts.get("N_full")), "{}/{}".format(_metric_text(neff.get("baseline_w0")), _metric_text(neff.get("method_a_w0_times_C"))), "{}/{}".format((support.get("prompt_control", {}) or {}).get("ood_count"), (support.get("full_physical_application", {}) or {}).get("ood_count")), _metric_text(metrics.get("DeltaH")), _metric_text(metrics.get("R")), _metric_text(metrics.get("kappa")), _metric_text(x_kappa), _metric_text(y_kappa), _metric_text(window.get("DeltaP_K")), _metric_text(variance.get("R_V"))])
    table = axis.table(cellText=rows, colLabels=headings, loc="center", cellLoc="center")
    table.auto_set_font_size(False); table.set_fontsize(6.5); table.scale(1.0, 1.55)
    axis.text(0.01, 0.03, "Descriptive diagnostics only. No automatic case or acceptance classification is assigned.", transform=axis.transAxes, fontsize=8)
    pdf.savefig(figure); plt.close(figure)


def _detail_page(pdf: object, child: Mapping[str, object], window: Mapping[str, object]) -> None:
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(3, 3, figsize=(14, 10)); figure.suptitle("F.6.2 {} t{} phi{} [{:.0f}, {:.0f})".format(child.get("setting_id"), child.get("canonical_t_index"), child.get("phi_index"), float(child.get("phi_low", 0.0)), float(child.get("phi_high", 0.0))), fontsize=13, fontweight="bold")
    shapes = child.get("one_dimensional", {}); variables = ("analysis_MM", "SHMS_xptar", "SHMS_yptar")
    for axis, variable in zip(axes[0], variables):
        payload = shapes.get(variable, {}) if isinstance(shapes, Mapping) else {}; _overlay(axis, payload if isinstance(payload, Mapping) else {}, variable)
        if variable == "analysis_MM":
            axis.axvspan(float(window["mm_min"]), float(window["mm_max"]), color="#d9d9d9", alpha=0.6, label="frozen K window")
    joint = child.get("joint_distributions", {}); map_keys = ("SHMS_delta__SHMS_xptar", "SHMS_delta__SHMS_yptar")
    for row_axes, map_key in zip(axes[1:], map_keys):
        payload = joint.get(map_key, {}) if isinstance(joint, Mapping) else {}; payload = payload if isinstance(payload, Mapping) else {}
        populations = [payload.get(name, {}) for name in ("L", "B", "A")]
        values = [np.asarray(item.get("unit_area"), dtype=float) for item in populations if isinstance(item, Mapping)]
        vmax = max([float(np.max(value)) for value in values if value.size and np.all(np.isfinite(value))] or [1.0])
        for axis, name in zip(row_axes, ("L", "B", "A")):
            _matrix(axis, payload, name, "{} — {}".format(map_key, name), vmax)
    mm = shapes.get("analysis_MM", {}) if isinstance(shapes, Mapping) else {}
    metrics = mm.get("metrics", {}) if isinstance(mm, Mapping) else {}
    signed = child.get("signed_background", {})
    bootstrap = child.get("bootstrap", {})
    bootstrap_one_dimensional = bootstrap.get("one_dimensional", {}) if isinstance(bootstrap, Mapping) else {}
    bootstrap_mm = bootstrap_one_dimensional.get("analysis_MM", {}) if isinstance(bootstrap_one_dimensional, Mapping) else {}
    diagnostics = "MM: DeltaH={} R={} kappa={} rho={}; K window: {}; bootstrap kappa: {}".format(
        _metric_text(metrics.get("DeltaH")), _metric_text(metrics.get("R")), _metric_text(metrics.get("kappa")), _metric_text(metrics.get("rho")),
        _metric_text((signed.get("kaon_window", {}) if isinstance(signed, Mapping) else {}).get("DeltaP_K")),
        _metric_text(bootstrap_mm.get("kappa", {}) if isinstance(bootstrap_mm, Mapping) else {}),
    )
    figure.text(0.02, 0.015, diagnostics, fontsize=7.5)
    figure.tight_layout(rect=(0, 0.04, 1, 0.94)); pdf.savefig(figure); plt.close(figure)


def write_review_pdf(path: Path, artifact: Mapping[str, object]) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    from matplotlib.backends.backend_pdf import PdfPages

    validation_payload = artifact.get("validation")
    if not isinstance(validation_payload, Mapping):
        raise ValueError("f6_2_pdf_validation_missing")
    window = validation_payload.get("kaon_window")
    if not isinstance(window, Mapping):
        raise ValueError("f6_2_pdf_window_missing")
    parents = validation_payload.get("parents")
    if not isinstance(parents, Sequence):
        raise ValueError("f6_2_pdf_parents_missing")
    with PdfPages(path) as pdf:
        authority = validation_payload.get("runtime_authority", {}); f6_authority = validation_payload.get("f6_1_authority", {}); policy = validation_payload.get("bootstrap_policy", {})
        _text_page(pdf, "F.6.2 — acceptance-correlated Method-A refinement validation", [
            "Detached aggregate-only review. Baseline disagreement is not automatic failure.",
            "Low-HGCer is an observed reference shape, not an absolute leakage probability.",
            "Pion leakage in the kaon region is permitted and reviewed against MM/acceptance evidence.",
            "F.6.1 source SHA: {}".format(((f6_authority.get("accepted", {}) if isinstance(f6_authority, Mapping) else {}).get("source_file_sha256"))),
            "F.6.1 validation/artifact fingerprints: {} / {}".format(((f6_authority.get("accepted", {}) if isinstance(f6_authority, Mapping) else {}).get("validation_fingerprint")), ((f6_authority.get("accepted", {}) if isinstance(f6_authority, Mapping) else {}).get("artifact_fingerprint"))),
            "Upstream authority: {}".format(authority),
            "Bootstrap policy: {}".format(policy),
            "Fixed kaon window: [{}, {})".format(window.get("mm_min"), window.get("mm_max")),
            "Method B is not a numerical input; no production object, yield, or threshold gate is constructed.",
        ])
        for parent in parents:
            if isinstance(parent, Mapping):
                _parent_overview_page(pdf, parent)
        detail_count = 0
        for parent in parents:
            if not isinstance(parent, Mapping):
                continue
            for child in parent.get("children", []):
                if isinstance(child, Mapping):
                    counts = child.get("population_counts", {})
                    if isinstance(counts, Mapping) and any(int(counts.get(name, 0)) > 0 for name in ("N_low", "N_control", "N_full")):
                        _detail_page(pdf, child, window); detail_count += 1
        empty = sum(1 for parent in parents if isinstance(parent, Mapping) for child in parent.get("children", []) if isinstance(child, Mapping) and not any(int((child.get("population_counts", {}) or {}).get(name, 0)) > 0 for name in ("N_low", "N_control", "N_full")))
        _text_page(pdf, "F.6.2 — policy and availability summary", [
            "Canonical children: 135; detail pages: {}; empty children: {}.".format(detail_count, empty),
            "No composite score, automatic case assignment, numerical threshold, or post-hoc tuning is present.",
            "Unavailable metrics retain their literal reason in JSON and detail-page footer diagnostics.",
            "F.6.2 is detached evidence only; it neither mutates production nor claims final-yield uncertainty reduction.",
        ])


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path); parser.add_argument("--kinematic", required=True)
    parser.add_argument("--output-json", required=True, type=Path); parser.add_argument("--output-pdf", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None, *, accepted_runtime_authority_by_kinematic: Mapping[str, object] | None = None, accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None = None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None, accepted_f6_1_artifact_authority_by_kinematic: Mapping[str, object] | None = None, bootstrap_test_config: Mapping[str, object] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    if arguments.kinematic != KINEMATIC:
        print("error: f6_2_kinematic_unsupported", file=sys.stderr); return 1
    try:
        expected_json = (arguments.outdir / validation.pion_hgcer_method_a_acceptance_refinement_validation_filename(arguments.kinematic)).resolve(); expected_pdf = (arguments.outdir / _pdf_filename(arguments.kinematic)).resolve()
        output_json, output_pdf = arguments.output_json.resolve(), arguments.output_pdf.resolve()
        if output_json != expected_json or output_pdf != expected_pdf:
            print("error: f6_2_output_path_not_deterministic", file=sys.stderr); return 2
        if output_json == output_pdf:
            print("error: f6_2_output_paths_collide", file=sys.stderr); return 2
        if output_json.exists() or output_pdf.exists():
            print("error: f6_2_output_path_already_exists", file=sys.stderr); return 2
    except (OSError, ValueError, validation.MethodAAcceptanceRefinementValidationError) as exc:
        print("error: {}".format(exc), file=sys.stderr); return 1
    temporary_json: Path | None = None; temporary_pdf: Path | None = None; promoted_json = False; promoted_pdf = False
    try:
        f1, hashes, f1_paths, f3, f3_sha, f3_path, f4, f4_sha, f4_path, f5, f5_sha, f5_path, f6, f6_sha, f6_path = load_inputs(arguments.outdir, arguments.kinematic)
        artifact = validation.build_pion_hgcer_method_a_acceptance_refinement_validation_artifact(f1, f3, f4, f5, f6, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, f4_input_file_sha256=f4_sha, f5_input_file_sha256=f5_sha, f6_1_input_file_sha256=f6_sha, input_paths={"f1": f1_paths, "f3": f3_path, "f4": f4_path, "f5": f5_path, "f6_1": f6_path}, generated_at_utc=_datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z"), accepted_runtime_authority_by_kinematic=accepted_runtime_authority_by_kinematic, accepted_f4_runtime_authority_by_kinematic=accepted_f4_runtime_authority_by_kinematic, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic, accepted_f6_1_artifact_authority_by_kinematic=accepted_f6_1_artifact_authority_by_kinematic, bootstrap_test_config=bootstrap_test_config)
        json_fd, json_name = tempfile.mkstemp(prefix=".f6_2_", suffix=".json", dir=arguments.outdir); os.close(json_fd); temporary_json = Path(json_name)
        pdf_fd, pdf_name = tempfile.mkstemp(prefix=".f6_2_", suffix=".pdf", dir=arguments.outdir); os.close(pdf_fd); temporary_pdf = Path(pdf_name)
        validation.write_pion_hgcer_method_a_acceptance_refinement_validation_json(temporary_json, artifact); write_review_pdf(temporary_pdf, artifact)
        os.replace(temporary_pdf, output_pdf); temporary_pdf = None; promoted_pdf = True
        os.replace(temporary_json, output_json); temporary_json = None; promoted_json = True
    except (validation.MethodAAcceptanceRefinementValidationError, OSError, ValueError) as exc:
        for path in (temporary_json, temporary_pdf):
            if path is not None and path.exists():
                path.unlink()
        if promoted_json and output_json.exists(): output_json.unlink()
        if promoted_pdf and output_pdf.exists(): output_pdf.unlink()
        print("error: {}".format(exc), file=sys.stderr); return 1
    print("wrote {}".format(output_json)); print("wrote {}".format(output_pdf)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
