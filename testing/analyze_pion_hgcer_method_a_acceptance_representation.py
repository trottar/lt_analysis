"""Run the detached global F.2 Method-A acceptance-representation audit.

This script reads exactly five F.1 v2 JSON artifacts from one output directory.
It neither opens ROOT files nor reruns or changes the KaonLT analysis.
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


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "cuts"))

from pion_hgcer_method_a_acceptance_representation import (  # noqa: E402
    CANDIDATE_REPRESENTATIONS,
    MethodAAcceptanceRepresentationError,
    build_pion_hgcer_method_a_acceptance_representation_artifact,
    pion_hgcer_method_a_acceptance_representation_filename,
    write_pion_hgcer_method_a_acceptance_representation_json,
)


CANONICAL_SETTINGS = (
    ("Left", "lowe"),
    ("Left", "highe"),
    ("Center", "lowe"),
    ("Center", "highe"),
    ("Right", "highe"),
)


def _f1_filename(phi: str, kinematic: str, epsilon: str) -> str:
    return "{}_kaon_pion-background_hgcer_method-a-acceptance-contract_{}_{}.json".format(
        phi, kinematic, epsilon
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_value(arguments: Sequence[str]) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *arguments], cwd=REPO_ROOT, text=True, capture_output=True,
            check=False,
        )
    except OSError:
        return None
    return completed.stdout.strip() if completed.returncode == 0 else None


def load_f1_inputs(outdir: Path, kinematic: str) -> tuple[list[dict[str, object]], dict[str, str], dict[str, str]]:
    """Load the exact, non-recursive five-file F.1 input set."""
    if not outdir.is_dir():
        raise ValueError("f2_outdir_invalid")
    artifacts, hashes, paths = [], {}, {}
    seen_paths: set[Path] = set()
    for phi, epsilon in CANONICAL_SETTINGS:
        setting_id = "{}-{}".format(phi, epsilon)
        path = outdir / _f1_filename(phi, kinematic, epsilon)
        if path in seen_paths:
            raise ValueError("f2_f1_input_duplicate")
        seen_paths.add(path)
        if not path.is_file():
            raise ValueError("f2_f1_input_missing:{}".format(setting_id))
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("f2_f1_input_json_invalid:{}".format(setting_id)) from exc
        if not isinstance(payload, dict):
            raise ValueError("f2_f1_input_json_invalid:{}".format(setting_id))
        artifacts.append(payload)
        hashes[setting_id] = _sha256(path)
        paths[setting_id] = str(path.resolve())
    return artifacts, hashes, paths


def _page_lines(lines: Sequence[str], *, title: str) -> object:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(11, 8.5))
    figure.suptitle(title, fontsize=15, fontweight="bold")
    figure.text(0.03, 0.94, "\n".join(lines), va="top", family="monospace", fontsize=7.5)
    return figure


def write_review_pdf(path: Path, artifact: Mapping[str, object]) -> None:
    """Write the four-page detached F.2 review PDF without event-level output."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt

    representation = artifact["representation"]
    if not isinstance(representation, Mapping):
        raise ValueError("f2_representation_invalid")
    groups = representation["groups"]
    summaries = representation["candidate_summaries"]
    recommendation = representation["recommendation"]
    if not isinstance(groups, list) or not isinstance(summaries, list) or not isinstance(recommendation, Mapping):
        raise ValueError("f2_representation_invalid")
    with PdfPages(path) as pdf:
        inventory = [
            "setting/t        positive   low  control  low_frac  nonprompt_app",
            "----------------------------------------------------------------",
        ]
        for group in groups:
            if not isinstance(group, Mapping):
                continue
            positive = int(group["training_positive_count"])
            low = int(group["low_count"])
            setting = group["setting"] if isinstance(group.get("setting"), Mapping) else {}
            name = "{}-{} t{}".format(setting.get("phi_setting"), setting.get("epsilon_filename_token"), group.get("canonical_t_index"))
            inventory.append("{:<20} {:>8} {:>5} {:>8} {:>9.4f} {:>14}".format(
                name, positive, low, int(group["control_count"]), low / positive if positive else 0.0,
                int(group["nonprompt_application_count"]),
            ))
        inventory.extend(("", "Diagnostic representation audit only.", "No map, correction, event probability, or pion-weight adjustment is constructed."))
        pdf.savefig(_page_lines(inventory, title="F2.1 — response-support inventory"))
        plt.close()

        information = ["candidate/group                       ROC AUC   AUC loss vs full5", "---------------------------------------------------------------"]
        for group in groups:
            if not isinstance(group, Mapping):
                continue
            setting = group["setting"] if isinstance(group.get("setting"), Mapping) else {}
            name = "{}-{} t{}".format(setting.get("phi_setting"), setting.get("epsilon_filename_token"), group.get("canonical_t_index"))
            metrics = group.get("candidate_metrics", {})
            if isinstance(metrics, Mapping):
                for candidate in CANDIDATE_REPRESENTATIONS:
                    metric = metrics.get(candidate["candidate_id"], {})
                    if isinstance(metric, Mapping):
                        auc = metric.get("roc_auc")
                        loss = metric.get("auc_loss_relative_to_full5_reference")
                        information.append("{:<26} {:<18} {:>8} {:>16}".format(
                            name, candidate["candidate_id"], "--" if auc is None else "{:.5f}".format(float(auc)),
                            "--" if loss is None else "{:.5f}".format(float(loss)),
                        ))
        information.extend(("", "Information gates: median AUC loss <= 0.02; max AUC loss <= 0.05."))
        pdf.savefig(_page_lines(information, title="F2.2 — response-information comparison"))
        plt.close()

        support = ["candidate/group                BLL penalty   app OOD frac   sparse", "---------------------------------------------------------------"]
        for group in groups:
            if not isinstance(group, Mapping):
                continue
            setting = group["setting"] if isinstance(group.get("setting"), Mapping) else {}
            name = "{}-{} t{}".format(setting.get("phi_setting"), setting.get("epsilon_filename_token"), group.get("canonical_t_index"))
            metrics = group.get("candidate_metrics", {})
            if isinstance(metrics, Mapping):
                for candidate in CANDIDATE_REPRESENTATIONS:
                    metric = metrics.get(candidate["candidate_id"], {})
                    if not isinstance(metric, Mapping):
                        continue
                    audit = metric.get("application_support", {})
                    penalty = metric.get("balanced_log_loss_penalty_relative_to_full5_reference")
                    ood = audit.get("application_ood_fraction") if isinstance(audit, Mapping) else None
                    sparse = audit.get("statistically_sparse") if isinstance(audit, Mapping) else None
                    support.append("{:<26} {:<18} {:>11} {:>14} {:>8}".format(
                        name, candidate["candidate_id"], "--" if penalty is None else "{:.5f}".format(float(penalty)),
                        "--" if ood is None else "{:.5f}".format(float(ood)), str(sparse),
                    ))
        support.extend(("", "Gates: median/max BLL penalty <= 0.01/0.02; OOD fraction <= 0.10.", "Groups with 1-19 non-prompt rows require manual review."))
        pdf.savefig(_page_lines(support, title="F2.3 — calibration probe and application support"))
        plt.close()

        decision = ["candidate                  info gate  support gate  overall  eligible"]
        decision.append("---------------------------------------------------------------------")
        for summary in summaries:
            if isinstance(summary, Mapping):
                decision.append("{:<26} {:>9} {:>13} {:>8} {:>9}".format(
                    str(summary.get("candidate_id")), str(summary.get("information_gate_passed")),
                    str(summary.get("application_support_gate_passed")),
                    str(summary.get("overall_candidate_passed")),
                    str(summary.get("automatic_recommendation_eligible")),
                ))
        decision.extend((
            "", "recommendation_status = {}".format(recommendation.get("recommendation_status")),
            "recommended_basis = {}".format(recommendation.get("recommended_basis")),
            "basis_frozen = false", "manual_review_required = true", "",
            "Diagnostic representation audit only.",
            "No Method-A map, correction, event probability, or pion-weight adjustment is constructed.",
        ))
        pdf.savefig(_page_lines(decision, title="F2.4 — representation recommendation"))
        plt.close()


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path, help="directory containing the five F.1 JSON artifacts")
    parser.add_argument("--kinematic", required=True, help="F.1 kinematic filename token, for example Q4p4W2p74")
    parser.add_argument("--output-json", required=True, type=Path, help="new deterministic F.2 JSON output path")
    parser.add_argument("--output-pdf", required=True, type=Path, help="new F.2 four-page review PDF path")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    expected_json_name = pion_hgcer_method_a_acceptance_representation_filename(arguments.kinematic)
    if arguments.output_json.name != expected_json_name:
        print("error: f2_output_json_basename_invalid", file=sys.stderr)
        return 2
    if arguments.output_json.exists() or arguments.output_pdf.exists():
        print("error: f2_output_path_already_exists", file=sys.stderr)
        return 2
    try:
        f1_artifacts, input_hashes, input_paths = load_f1_inputs(arguments.outdir, arguments.kinematic)
        artifact = build_pion_hgcer_method_a_acceptance_representation_artifact(
            f1_artifacts,
            input_file_hashes=input_hashes,
            input_paths=input_paths,
            generated_at_utc=_datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
            git_head=_git_value(("rev-parse", "HEAD")),
            git_status_short=_git_value(("status", "--short")),
        )
        write_pion_hgcer_method_a_acceptance_representation_json(arguments.output_json, artifact)
        write_review_pdf(arguments.output_pdf, artifact)
    except (MethodAAcceptanceRepresentationError, OSError, ValueError) as exc:
        print("error: {}".format(exc), file=sys.stderr)
        return 1
    print("wrote {}".format(arguments.output_json))
    print("wrote {}".format(arguments.output_pdf))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
