"""Run detached Phase F.5 signed ``(t, phi)`` propagation from accepted inputs."""

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
import pion_hgcer_method_a_tphi_propagation as propagation  # noqa: E402


CANONICAL_SETTINGS = (("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"), ("Center", "highe"), ("Right", "highe"))
BASELINE_PION_BACKGROUND_LABEL = "Baseline pion background"
METHOD_A_PION_BACKGROUND_LABEL = "Method-A pion background"
PION_BACKGROUND_YIELD_LABEL = "Pion-background yield [weighted events]"
METHOD_A_CHANGE_LABEL = "Method-A change / t-bin baseline"
FRACTIONAL_T_BIN_CHANGE_LABEL = "Fraction of baseline t-bin pion background"


def _f1_filename(phi: str, kinematic: str, epsilon: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-contract_{}_{}.json".format(phi, kinematic, epsilon)


def _f3_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-map.json".format(kinematic)


def _f4_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-parent-preserving-correction.json".format(kinematic)


def _pdf_filename(kinematic: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-tphi-propagation.pdf".format(kinematic)


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


def _read_json(path: Path, label: str) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError("{}_json_invalid".format(label)) from exc
    if not isinstance(payload, dict):
        raise ValueError("{}_json_invalid".format(label))
    return payload


def load_inputs(outdir: Path, kinematic: str) -> tuple[list[dict[str, object]], dict[str, str], dict[str, str], dict[str, object], str, str, dict[str, object], str, str]:
    if not outdir.is_dir():
        raise ValueError("f5_outdir_invalid")
    artifacts: list[dict[str, object]] = []; hashes: dict[str, str] = {}; paths: dict[str, str] = {}; seen: set[Path] = set()
    for phi, epsilon in CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon); path = outdir / _f1_filename(phi, kinematic, epsilon)
        if path in seen:
            raise ValueError("f5_f1_input_duplicate")
        seen.add(path)
        if not path.is_file():
            raise ValueError("f5_f1_input_missing:{}".format(setting_id))
        payload = _read_json(path, "f5_f1_input")
        setting = payload.get("setting")
        if not isinstance(setting, Mapping) or setting.get("phi_setting") != phi or setting.get("epsilon_filename_token") != epsilon or setting.get("kinematic_token") != kinematic or setting.get("particle_type") != "kaon":
            raise ValueError("f5_f1_input_identity_invalid:{}".format(setting_id))
        artifacts.append(payload); hashes[setting_id] = _sha256(path); paths[setting_id] = str(path.resolve())
    f3_path = outdir / _f3_filename(kinematic); f4_path = outdir / _f4_filename(kinematic)
    if f3_path in seen or f4_path in seen or f3_path == f4_path:
        raise ValueError("f5_upstream_input_path_collides")
    if not f3_path.is_file() or not f4_path.is_file():
        raise ValueError("f5_upstream_input_missing")
    f3 = _read_json(f3_path, "f5_f3_input"); f4 = _read_json(f4_path, "f5_f4_input")
    return artifacts, hashes, paths, f3, _sha256(f3_path), str(f3_path.resolve()), f4, _sha256(f4_path), str(f4_path.resolve())


def _text_page(lines: Sequence[str], title: str) -> object:
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    figure = plt.figure(figsize=(11, 8.5)); figure.suptitle(title, fontsize=15, fontweight="bold")
    figure.text(0.03, 0.94, "\n".join(lines), va="top", family="monospace", fontsize=7.5)
    return figure


def _setting_rows(result: Mapping[str, object], setting_id: str) -> Mapping[str, object]:
    rows = result.get("setting_templates")
    matches = [row for row in rows if isinstance(row, Mapping) and row.get("setting_id") == setting_id] if isinstance(rows, list) else []
    if len(matches) != 1:
        raise ValueError("f5_pdf_setting_inventory_invalid")
    return matches[0]


def _t_bin_label(edges: np.ndarray, index: int) -> str:
    return "{:.3f} < t < {:.3f}: pion-background yield".format(float(edges[index]), float(edges[index + 1]))


def _phi_bin_interval_label(phi_edges: object, index: object) -> str:
    """Render a physical phi interval from the persisted F.5 geometry."""
    edges = np.asarray(phi_edges, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges)) or np.any(edges[1:] <= edges[:-1]):
        raise ValueError("f5_pdf_phi_edges_invalid")
    if isinstance(index, bool) or not isinstance(index, (int, np.integer)) or not 0 <= int(index) < len(edges) - 1:
        raise ValueError("f5_pdf_phi_index_invalid")
    return "[{:.0f}, {:.0f}] deg".format(float(edges[int(index)]), float(edges[int(index) + 1]))


def write_review_pdf(path: Path, artifact: Mapping[str, object]) -> None:
    """Render the fixed twelve-page, aggregate-only F.5 review package."""
    import matplotlib
    matplotlib.use("Agg", force=True)
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    result = artifact.get("propagation")
    if not isinstance(result, Mapping) or result.get("available") is not True:
        raise ValueError("f5_pdf_propagation_unavailable")
    with PdfPages(path) as pdf:
        lines = ["F.5 Method-A pion-background redistribution in (t,phi)", "", "Baseline pion background: sum of the signed pion-background event weights before Method A in each (t,phi) bin.", "Method-A pion background: the same contribution after the validated F.4 Method-A event correction.", "Method-A change: (Method-A pion background - baseline pion background) / total baseline pion background in the parent t bin.", "The total pion-background yield in each t bin is preserved by construction.", "Weighted yields may be negative because the established subtraction chain is signed.", "", "F4 source SHA-256 = {}".format(result["f4_source_file_sha256"]), "F4 correction fingerprint = {}".format(result["f4_correction_fingerprint"]), "Detached review only: no child renormalization, smoothing, interpolation, Method B, ROOT, yield, or production application."]
        pdf.savefig(_text_page(lines, "F5.1 — authority and physical interpretation")); plt.close()
        for phi, epsilon in CANONICAL_SETTINGS:
            setting_id = "{}-{}".format(phi, epsilon); row = _setting_rows(result, setting_id); edges = np.asarray(row["phi_edges"], dtype=float); t_edges = np.asarray(row["t_edges"], dtype=float)
            baseline = np.asarray(row["baseline_signed_contents"], dtype=float); adjusted = np.asarray(row["adjusted_signed_contents"], dtype=float); redistribution = np.asarray(row["redistribution_fraction_of_parent"], dtype=float)
            figure, axes = plt.subplots(3, 2, figsize=(12, 9)); figure.suptitle("F5 — {} pion-background yield vs phi".format(setting_id), fontsize=14, fontweight="bold")
            for index in range(3):
                axes[index, 0].stairs(baseline[index], edges, label=BASELINE_PION_BACKGROUND_LABEL, color="tab:blue"); axes[index, 0].stairs(adjusted[index], edges, label=METHOD_A_PION_BACKGROUND_LABEL, color="tab:orange")
                axes[index, 0].axhline(0.0, color="black", linewidth=0.6); axes[index, 0].set_title(_t_bin_label(t_edges, index)); axes[index, 0].set_xlabel("phi [degrees]"); axes[index, 0].set_ylabel(PION_BACKGROUND_YIELD_LABEL); axes[index, 0].legend(fontsize=7)
                axes[index, 1].stairs(redistribution[index], edges, color="tab:purple"); axes[index, 1].axhline(0.0, color="black", linewidth=0.6); axes[index, 1].set_title(METHOD_A_CHANGE_LABEL); axes[index, 1].set_xlabel("phi [degrees]"); axes[index, 1].set_ylabel(FRACTIONAL_T_BIN_CHANGE_LABEL)
            figure.tight_layout(rect=(0, 0, 1, 0.95)); pdf.savefig(figure); plt.close(figure)
        for phi, epsilon in CANONICAL_SETTINGS:
            setting_id = "{}-{}".format(phi, epsilon); row = _setting_rows(result, setting_id); t_edges = np.asarray(row["t_edges"], dtype=float); phi_edges = np.asarray(row["phi_edges"], dtype=float)
            baseline = np.asarray(row["baseline_signed_contents"], dtype=float); adjusted = np.asarray(row["adjusted_signed_contents"], dtype=float); redistribution = np.asarray(row["redistribution_fraction_of_parent"], dtype=float)
            scale = max(float(np.max(np.abs(baseline))), float(np.max(np.abs(adjusted))), 1.0e-300); rscale = max(float(np.max(np.abs(redistribution))), 1.0e-300)
            figure, axes = plt.subplots(1, 3, figsize=(15, 5)); figure.suptitle("F5 — {} pion-background yield and Method-A redistribution".format(setting_id), fontsize=14, fontweight="bold")
            for axis, values, title, colorbar_label, cmap, low, high in ((axes[0], baseline, BASELINE_PION_BACKGROUND_LABEL, PION_BACKGROUND_YIELD_LABEL, "coolwarm", -scale, scale), (axes[1], adjusted, METHOD_A_PION_BACKGROUND_LABEL, PION_BACKGROUND_YIELD_LABEL, "coolwarm", -scale, scale), (axes[2], redistribution, METHOD_A_CHANGE_LABEL, FRACTIONAL_T_BIN_CHANGE_LABEL, "PiYG", -rscale, rscale)):
                mesh = axis.pcolormesh(phi_edges, t_edges, values, shading="flat", cmap=cmap, vmin=low, vmax=high); figure.colorbar(mesh, ax=axis).set_label(colorbar_label); axis.set_title(title); axis.set_xlabel("phi [degrees]"); axis.set_ylabel("t")
            figure.tight_layout(rect=(0, 0, 1, 0.93)); pdf.savefig(figure); plt.close(figure)
        phi_edges = result.get("phi_edges")
        closure = ["setting / t             baseline pion yield  Method-A pion yield  Method-A - baseline", "setting / t             max |change / t-bin baseline|  sum |change / t-bin baseline|  phi interval of max change [deg]", "------------------------------------------------------------------------------------------------------------------------"]
        for row in result.get("parent_closure", []):
            if not isinstance(row, Mapping):
                raise ValueError("f5_pdf_parent_inventory_invalid")
            phi_interval = _phi_bin_interval_label(phi_edges, row["max_abs_redistribution_phi_index"])
            closure.append("{:<23} {:>11.5g} {:>11.5g} {:>11.3g} {:>12.3g} {:>12.3g} {:>29}".format("{} t bin {}".format(row["setting_id"], int(row["canonical_t_index"]) + 1), row["baseline_parent_sum"], row["adjusted_parent_sum"], row["closure_residual"], row["max_abs_redistribution"], row["sum_abs_redistribution"], phi_interval))
        closure.extend(("", "Detached review only. No production template, correction application, Method B, ROOT, yield, or cross section exists."))
        pdf.savefig(_text_page(closure, "F5.12 — fifteen-parent closure")); plt.close()


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path); parser.add_argument("--kinematic", required=True)
    parser.add_argument("--output-json", required=True, type=Path); parser.add_argument("--output-pdf", required=True, type=Path)
    return parser


def main(argv: Sequence[str] | None = None, *, accepted_f4_runtime_authority_by_kinematic: Mapping[str, object] | None = None, accepted_f3_runtime_authority_by_kinematic: Mapping[str, object] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    expected_json = (arguments.outdir / propagation.pion_hgcer_method_a_tphi_propagation_filename(arguments.kinematic)).resolve(); expected_pdf = (arguments.outdir / _pdf_filename(arguments.kinematic)).resolve()
    output_json, output_pdf = arguments.output_json.resolve(), arguments.output_pdf.resolve()
    if output_json != expected_json or output_pdf != expected_pdf:
        print("error: f5_output_path_not_deterministic", file=sys.stderr); return 2
    if output_json == output_pdf:
        print("error: f5_output_paths_collide", file=sys.stderr); return 2
    if output_json.exists() or output_pdf.exists():
        print("error: f5_output_path_already_exists", file=sys.stderr); return 2
    try:
        f1, hashes, f1_paths, f3, f3_sha, f3_path, f4, f4_sha, f4_path = load_inputs(arguments.outdir, arguments.kinematic)
        artifact = propagation.build_pion_hgcer_method_a_tphi_propagation_artifact(f1, f3, f4, f1_input_file_hashes=hashes, f3_input_file_sha256=f3_sha, f4_input_file_sha256=f4_sha, input_paths={"f1": f1_paths, "f3": f3_path, "f4": f4_path}, generated_at_utc=_datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z"), git_head=_git_value(("rev-parse", "HEAD")), git_status_short=_git_value(("status", "--short")), accepted_f4_runtime_authority_by_kinematic=accepted_f4_runtime_authority_by_kinematic, accepted_f3_runtime_authority_by_kinematic=accepted_f3_runtime_authority_by_kinematic)
        propagation.write_pion_hgcer_method_a_tphi_propagation_json(output_json, artifact); write_review_pdf(output_pdf, artifact)
    except (propagation.MethodATPhiPropagationError, OSError, ValueError) as exc:
        print("error: {}".format(exc), file=sys.stderr); return 1
    print("wrote {}".format(output_json)); print("wrote {}".format(output_pdf)); return 0


if __name__ == "__main__":
    raise SystemExit(main())
