"""Run the detached Phase F.3 Method-A ``hgcer3`` relative response map.

The analyzer reads the exact five F.1 v2 JSON artifacts and one accepted F.2
representation from one directory.  It does not import ROOT, rerun KaonLT, or
modify analysis inputs.
"""

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

import pion_hgcer_method_a_acceptance_map as acceptance_map  # noqa: E402


CANONICAL_SETTINGS = (
    ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
    ("Center", "highe"), ("Right", "highe"),
)


def _f1_filename(phi: str, kinematic: str, epsilon: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-contract_{}_{}.json".format(phi, kinematic, epsilon)


def _f2_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-representation.json".format(kinematic)


def _pdf_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-map.pdf".format(kinematic)


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


def _validate_f1_filename_identity(payload: Mapping[str, object], *, phi: str, epsilon: str, kinematic: str, setting_id: str) -> None:
    setting = payload.get("setting")
    if not isinstance(setting, Mapping):
        raise ValueError("f3_f1_input_identity_setting_invalid:{}".format(setting_id))
    for name, expected in (("phi_setting", phi), ("epsilon_filename_token", epsilon), ("kinematic_token", kinematic), ("particle_type", "kaon")):
        if setting.get(name) != expected:
            raise ValueError("f3_f1_input_identity_{}:{}".format(name, setting_id))


def load_inputs(outdir: Path, kinematic: str) -> tuple[list[dict[str, object]], dict[str, str], dict[str, str], dict[str, object], str, str]:
    """Load the exact six-file F.3 authority closure without discovery/fallbacks."""
    if not outdir.is_dir():
        raise ValueError("f3_outdir_invalid")
    artifacts: list[dict[str, object]] = []
    hashes: dict[str, str] = {}
    paths: dict[str, str] = {}
    seen: set[Path] = set()
    for phi, epsilon in CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon)
        path = outdir / _f1_filename(phi, kinematic, epsilon)
        if path in seen:
            raise ValueError("f3_f1_input_duplicate")
        seen.add(path)
        if not path.is_file():
            raise ValueError("f3_f1_input_missing:{}".format(setting_id))
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("f3_f1_input_json_invalid:{}".format(setting_id)) from exc
        if not isinstance(payload, dict):
            raise ValueError("f3_f1_input_json_invalid:{}".format(setting_id))
        _validate_f1_filename_identity(payload, phi=phi, epsilon=epsilon, kinematic=kinematic, setting_id=setting_id)
        artifacts.append(payload); hashes[setting_id] = _sha256(path); paths[setting_id] = str(path.resolve())
    f2_path = outdir / _f2_filename(kinematic)
    if f2_path in seen:
        raise ValueError("f3_f2_input_path_collides")
    if not f2_path.is_file():
        raise ValueError("f3_f2_input_missing")
    try:
        f2_payload = json.loads(f2_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError("f3_f2_input_json_invalid") from exc
    if not isinstance(f2_payload, dict):
        raise ValueError("f3_f2_input_json_invalid")
    return artifacts, hashes, paths, f2_payload, _sha256(f2_path), str(f2_path.resolve())


def _page_lines(lines: Sequence[str], title: str) -> object:
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(11, 8.5))
    figure.suptitle(title, fontsize=15, fontweight="bold")
    figure.text(0.03, 0.94, "\n".join(lines), va="top", family="monospace", fontsize=7.5)
    return figure


def _training_by_parent(f1_artifacts: Sequence[Mapping[str, object]]) -> dict[tuple[str, int], list[Mapping[str, object]]]:
    parsed = acceptance_map._validate_f1_artifacts(f1_artifacts)
    result: dict[tuple[str, int], list[Mapping[str, object]]] = {}
    for artifact in parsed:
        setting_id = str(artifact["setting_id"])
        for row in acceptance_map._sequence(artifact["training"], "training"):
            item = acceptance_map._mapping(row, "training")
            result.setdefault((setting_id, int(item["t_index"])), []).append(item)
    return result


def _masked_slice(model: Mapping[str, object], training: Sequence[Mapping[str, object]], x_index: int, y_index: int, fixed_index: int, side: int) -> tuple[np.ndarray, np.ndarray, np.ma.MaskedArray]:
    minimum = np.asarray(model["training_feature_minimum"], dtype=float)
    maximum = np.asarray(model["training_feature_maximum"], dtype=float)
    median = np.asarray(model["scaler"]["median"], dtype=float)  # type: ignore[index]
    x_axis = np.linspace(minimum[x_index], maximum[x_index], side)
    y_axis = np.linspace(minimum[y_index], maximum[y_index], side)
    x_grid, y_grid = np.meshgrid(x_axis, y_axis)
    raw = np.tile(median, (side * side, 1))
    raw[:, x_index] = x_grid.ravel(); raw[:, y_index] = y_grid.ravel(); raw[:, fixed_index] = median[fixed_index]
    training_values = np.asarray([[float(row[name]) for name in acceptance_map.ACCEPTED_FEATURES] for row in training], dtype=float)
    response, in_support = acceptance_map.evaluate_relative_response_grid(model, raw, training_values)
    return x_grid, y_grid, np.ma.masked_where(~in_support.reshape(side, side), response.reshape(side, side))


def _setting_page(pdf: object, setting_models: Sequence[Mapping[str, object]], parent_training: Mapping[tuple[str, int], Sequence[Mapping[str, object]]], side: int) -> None:
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(3, 4, figsize=(15, 11))
    setting = setting_models[0]["setting"]
    figure.suptitle("F3 support-aware relative response — {} {}".format(setting["phi_setting"], setting["epsilon_filename_token"]), fontsize=14, fontweight="bold")  # type: ignore[index]
    slices = ((0, 1, 2, "SHMS_delta", "P_hgcer_xAtCer", "P_hgcer_yAtCer"), (0, 2, 1, "SHMS_delta", "P_hgcer_yAtCer", "P_hgcer_xAtCer"), (1, 2, 0, "P_hgcer_xAtCer", "P_hgcer_yAtCer", "SHMS_delta"))
    for row, model in enumerate(setting_models):
        axis = axes[row, 0]
        support = model["application_support"]  # type: ignore[index]
        fit = model["logistic_fit"]  # type: ignore[index]
        axis.axis("off")
        axis.text(0.0, 1.0, "t{t}\n[{low:.6g}, {high:.6g}]\nlow/control = {low_count}/{control_count}\nobjective = {objective:.6g}\nthreshold = {threshold:.6g}\nnonprompt = {nonprompt}\nOOD = {ood} ({fraction:.4f})\nF2 continuity = {continuity}".format(t=model["canonical_t_index"], low=model["canonical_t_low"], high=model["canonical_t_high"], low_count=model["low_count"], control_count=model["control_count"], objective=fit["objective"], threshold=support["support_distance_threshold"], nonprompt=support["nonprompt_application_count"], ood=support["application_ood_count"], fraction=support["application_ood_fraction"], continuity=model["f2_support_continuity"]["passed"]), va="top", family="monospace", fontsize=8)  # type: ignore[index]
        key = (str(model["setting_id"]), int(model["canonical_t_index"]))
        for column, (x_index, y_index, fixed_index, x_label, y_label, fixed_label) in enumerate(slices, start=1):
            plot = axes[row, column]
            x_grid, y_grid, response = _masked_slice(model, parent_training[key], x_index, y_index, fixed_index, side)
            contour = plot.contourf(x_grid, y_grid, response, levels=20)
            figure.colorbar(contour, ax=plot, fraction=0.046, pad=0.04)
            plot.set_title("t{}; {} at median {}".format(model["canonical_t_index"], fixed_label, fixed_label), fontsize=8)
            plot.set_xlabel(x_label, fontsize=7); plot.set_ylabel(y_label, fontsize=7)
    figure.text(0.5, 0.012, "Color is exp(beta dot z), relative only. White regions are outside the cKDTree p99 support mask. No absolute probability or correction is shown.", ha="center", fontsize=8)
    pdf.savefig(figure); plt.close(figure)


def write_review_pdf(path: Path, artifact: Mapping[str, object], f1_artifacts: Sequence[Mapping[str, object]]) -> None:
    """Write the seven-page F.3 review PDF without serializing plot grids."""
    import matplotlib
    matplotlib.use("Agg", force=True)
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    result = artifact.get("acceptance_map")
    if not isinstance(result, Mapping) or result.get("available") is not True:
        raise ValueError("f3_pdf_map_unavailable")
    models = result.get("models")
    if not isinstance(models, list) or len(models) != 15:
        raise ValueError("f3_pdf_models_invalid")
    parent_training = _training_by_parent(f1_artifacts)
    side = int(result["algorithm_config"]["review_grid_side"])  # type: ignore[index]
    with PdfPages(path) as pdf:
        summary = ["Detached F.3 support-aware relative Method-A response map", "accepted basis = hgcer3 = (SHMS_delta, P_hgcer_xAtCer, P_hgcer_yAtCer)", "parents = 15 (five settings x three canonical-t bins)", "F2 representation fingerprint = {}".format(result["f2_representation_fingerprint"]), "F2 source SHA-256 = {}".format(result["f2_source_file_sha256"]), "F3 fingerprint = {}".format(result["fingerprint"]), "", "Map definition: exp(beta dot z); fitted intercept retained only as provenance.", "No absolute probability, parent normalization, correction, event output, production application, or Method-B numerical input exists.", "All plotted grids are transient and cKDTree p99-masked."]
        pdf.savefig(_page_lines(summary, "F3.1 — global map authority and boundaries")); plt.close()
        for phi, epsilon in CANONICAL_SETTINGS:
            setting_models = [model for model in models if isinstance(model, Mapping) and isinstance(model.get("setting"), Mapping) and model["setting"].get("phi_setting") == phi and model["setting"].get("epsilon_filename_token") == epsilon]
            if len(setting_models) != 3:
                raise ValueError("f3_pdf_setting_model_inventory_invalid")
            _setting_page(pdf, sorted(setting_models, key=lambda item: int(item["canonical_t_index"])), parent_training, side)
        closure = ["setting/t               low/control  nonprompt  OOD fraction  support threshold  continuity", "--------------------------------------------------------------------------------"]
        for model in models:
            if not isinstance(model, Mapping):
                continue
            setting = model["setting"]
            support = model["application_support"]
            closure.append("{:<24} {:>4}/{:<4} {:>9} {:>13.6f} {:>18.9g} {:>11}".format("{}-{} t{}".format(setting["phi_setting"], setting["epsilon_filename_token"], model["canonical_t_index"]), model["low_count"], model["control_count"], support["nonprompt_application_count"], support["application_ood_fraction"], support["support_distance_threshold"], str(model["f2_support_continuity"]["passed"])))
        closure.extend(("", "All fifteen parents must be valid before this artifact is available.", "Detached review artifact only; F.4 remains blocked pending farm review."))
        pdf.savefig(_page_lines(closure, "F3.7 — fifteen-parent closure table")); plt.close()


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path, help="directory containing the F.1/F.2 inputs")
    parser.add_argument("--kinematic", required=True, help="deterministic F.1/F.2 kinematic filename token")
    parser.add_argument("--output-json", type=Path, help="optional exact deterministic F.3 JSON output path")
    parser.add_argument("--output-pdf", type=Path, help="optional exact deterministic seven-page F.3 PDF path")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    expected_json = (arguments.outdir / acceptance_map.pion_hgcer_method_a_acceptance_map_filename(arguments.kinematic)).resolve()
    expected_pdf = (arguments.outdir / _pdf_filename(arguments.kinematic)).resolve()
    output_json = expected_json if arguments.output_json is None else arguments.output_json.resolve()
    output_pdf = expected_pdf if arguments.output_pdf is None else arguments.output_pdf.resolve()
    if output_json != expected_json or output_pdf != expected_pdf:
        print("error: f3_output_path_not_deterministic", file=sys.stderr); return 2
    if output_json == output_pdf:
        print("error: f3_output_paths_collide", file=sys.stderr); return 2
    if output_json.exists() or output_pdf.exists():
        print("error: f3_output_path_already_exists", file=sys.stderr); return 2
    try:
        f1_artifacts, f1_hashes, f1_paths, f2_artifact, f2_sha256, f2_path = load_inputs(arguments.outdir, arguments.kinematic)
        input_paths: dict[str, object] = {"f1": f1_paths, "f2": f2_path}
        artifact = acceptance_map.build_pion_hgcer_method_a_acceptance_map_artifact(f1_artifacts, f2_artifact, f1_input_file_hashes=f1_hashes, f2_input_file_sha256=f2_sha256, input_paths=input_paths, generated_at_utc=_datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z"), git_head=_git_value(("rev-parse", "HEAD")), git_status_short=_git_value(("status", "--short")))
        if artifact["acceptance_map"]["available"] is not True:  # type: ignore[index]
            raise ValueError("f3_map_unavailable")
        acceptance_map.write_pion_hgcer_method_a_acceptance_map_json(output_json, artifact)
        write_review_pdf(output_pdf, artifact, f1_artifacts)
    except (acceptance_map.MethodAAcceptanceMapError, OSError, ValueError) as exc:
        print("error: {}".format(exc), file=sys.stderr); return 1
    print("wrote {}".format(output_json)); print("wrote {}".format(output_pdf))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
