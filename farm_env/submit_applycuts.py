#!/usr/bin/env python3
"""
Submit KaonLT applyCuts jobs to JLab SWIF2/ifarm.

Behavior:
- User supplies only q2 and w, e.g. `5p5 3p02`.
- The script discovers all matching JSON files under `input/kaon`.
- Each manifest variant keeps its own run list and submission context.
- Exactly one SWIF2 job is planned per manifest variant + run.
- A job is only planned when the full replay ROOT file already exists.
- A job is skipped when both kaon and pion skim outputs already exist.
- Each submitted job runs `applyCuts_Prod.sh` once, which in turn runs both
  kaon and pion passes for the requested run.

Dry-run is the default. Use --submit to create/modify the workflow and add jobs.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

try:
    from ltsep_paths import LtsepPaths, resolve_ltsep_paths
except ModuleNotFoundError:
    from farm_env.ltsep_paths import LtsepPaths, resolve_ltsep_paths

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_JSON_DIR = str(REPO_ROOT / "input" / "kaon")
DEFAULT_APPLYCUTS_SCRIPT = str(REPO_ROOT / "applyCuts_Prod.sh")
DEFAULT_WORKFLOW_PREFIX = "kaonlt"
DEFAULT_ACCOUNT = "c-kaonlt"
DEFAULT_PARTITION = "production"
DEFAULT_MAX_CONCURRENT = 200
DEFAULT_CORES = 1
DEFAULT_RAM = "6g"
DEFAULT_DISK = "25g"
DEFAULT_TIME = "4h"
DEFAULT_SWIF2 = os.environ.get("SWIF2_BIN", "swif2")

RUN_LINE_RE = re.compile(r"^\s*(\d+)\s*$")
SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9_.-]+")
VARIANT_RE = re.compile(
    r"^Q(?P<q2>[0-9p]+)W(?P<w>[0-9p]+)"
    r"(?P<phi>center|left|right)_(?P<epsilon>high|low)e(?P<dummy>_dummy)?$"
)


@dataclass(frozen=True)
class JsonVariant:
    path: Path
    family: str
    q2: str
    w: str
    phi: str
    epsilon: str
    target: str
    runs_files: Tuple[Path, ...]


@dataclass(frozen=True)
class RunPlan:
    variant: JsonVariant
    run: int
    job_name: str
    replay_file: Optional[Path]
    skim_kaon: Path
    skim_pion: Path
    status: str
    reason: str

    @property
    def should_add(self) -> bool:
        return self.status == "ADD"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit one KaonLT applyCuts SWIF2 job per manifest variant + run for a Q/W family."
    )
    parser.add_argument("q2", help="Q value without prefix, e.g. 5p5")
    parser.add_argument("w", help="W value without prefix, e.g. 3p02")
    parser.add_argument(
        "--json-dir",
        default=DEFAULT_JSON_DIR,
        help=f"Directory containing KaonLT JSON files (default: {DEFAULT_JSON_DIR})",
    )
    parser.add_argument(
        "--applycuts-script",
        default=DEFAULT_APPLYCUTS_SCRIPT,
        help=f"applyCuts wrapper to execute per variant/run (default: {DEFAULT_APPLYCUTS_SCRIPT})",
    )
    parser.add_argument(
        "--swif2-bin",
        default=DEFAULT_SWIF2,
        help="SWIF2 client executable (default: swif2 or $SWIF2_BIN)",
    )
    parser.add_argument(
        "--workflow-prefix",
        default=DEFAULT_WORKFLOW_PREFIX,
        help=f"Workflow name prefix (default: {DEFAULT_WORKFLOW_PREFIX})",
    )
    parser.add_argument(
        "--workflow-name",
        default=None,
        help="Override the auto-generated workflow name.",
    )
    parser.add_argument(
        "--account",
        default=DEFAULT_ACCOUNT,
        help=f"SWIF2 account / Slurm account (default: {DEFAULT_ACCOUNT})",
    )
    parser.add_argument(
        "--partition",
        default=DEFAULT_PARTITION,
        help="SWIF2 partition / Slurm partition (default: production).",
    )
    parser.add_argument(
        "--max-concurrent",
        type=int,
        default=DEFAULT_MAX_CONCURRENT,
        help=f"Workflow max concurrent jobs when creating a new workflow (default: {DEFAULT_MAX_CONCURRENT})",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=DEFAULT_CORES,
        help=f"CPU cores per job (default: {DEFAULT_CORES})",
    )
    parser.add_argument(
        "--ram",
        default=DEFAULT_RAM,
        help=f"RAM per job (default: {DEFAULT_RAM})",
    )
    parser.add_argument(
        "--disk",
        default=DEFAULT_DISK,
        help=f"Disk per job (default: {DEFAULT_DISK})",
    )
    parser.add_argument(
        "--time",
        default=DEFAULT_TIME,
        help=f"Wall time per job (default: {DEFAULT_TIME})",
    )
    parser.add_argument(
        "--family-regex",
        default=None,
        help="Optional regex applied to matching JSON basenames for additional filtering.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        default=True,
        help="Skip jobs whose names already exist in the target workflow (default: enabled).",
    )
    parser.add_argument(
        "--no-skip-existing",
        dest="skip_existing",
        action="store_false",
        help="Attempt to add jobs even if their names already exist in the workflow.",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Actually create/update the workflow and add jobs. Default is dry-run.",
    )
    parser.add_argument(
        "--no-run",
        action="store_true",
        help="With --submit, do not issue 'swif2 run' after adding jobs.",
    )
    return parser.parse_args()


def expand_path(path_text: str) -> Path:
    return Path(os.path.expandvars(path_text)).expanduser()


def normalize_family(q2: str, w: str) -> str:
    q2_text = str(q2).strip().lower()
    w_text = str(w).strip().lower()
    if q2_text.startswith("q"):
        q2_text = q2_text[1:]
    if w_text.startswith("w"):
        w_text = w_text[1:]
    return f"Q{q2_text}W{w_text}"


def safe_name(text: str) -> str:
    cleaned = SAFE_NAME_RE.sub("_", text).strip("._-")
    return cleaned or "kaonlt"


def load_json(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def extract_runs_files(data: Dict) -> List[Path]:
    jobs = data.get("jobs", [])
    runs_files: List[Path] = []
    if not isinstance(jobs, list):
        return runs_files
    for entry in jobs:
        if not isinstance(entry, dict):
            continue
        runs_file = entry.get("runs_file")
        if isinstance(runs_file, str) and runs_file.strip():
            runs_files.append(expand_path(runs_file))
    return runs_files


def parse_variant_name(path: Path) -> Optional[Dict[str, str]]:
    match = VARIANT_RE.match(path.stem)
    if match is None:
        return None
    info = match.groupdict()
    info["target"] = "dummy" if info.get("dummy") else "lh2"
    return info


def discover_json_variants(json_dir: Path, family_prefix: str, family_regex: Optional[str]) -> List[JsonVariant]:
    if not json_dir.exists():
        raise FileNotFoundError(f"JSON directory does not exist: {json_dir}")
    if not json_dir.is_dir():
        raise NotADirectoryError(f"JSON directory is not a directory: {json_dir}")

    regex = re.compile(family_regex) if family_regex else None
    candidates = sorted(json_dir.glob(f"{family_prefix}*.json"))
    variants: List[JsonVariant] = []

    for path in candidates:
        if regex and not regex.search(path.name):
            continue
        parsed = parse_variant_name(path)
        if parsed is None:
            continue
        data = load_json(path)
        runs_files = extract_runs_files(data)
        if not runs_files:
            continue
        variants.append(
            JsonVariant(
                path=path,
                family=path.stem,
                q2=parsed["q2"],
                w=parsed["w"],
                phi=parsed["phi"],
                epsilon=parsed["epsilon"],
                target=parsed["target"],
                runs_files=tuple(runs_files),
            )
        )

    return variants


def read_runs_file(path: Path) -> List[int]:
    if not path.exists():
        raise FileNotFoundError(f"Runs file does not exist: {path}")
    runs: List[int] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith(">"):
                continue
            match = RUN_LINE_RE.match(line)
            if match:
                runs.append(int(match.group(1)))
    if not runs:
        raise ValueError(f"No run numbers found in runs file: {path}")
    return runs


def runs_for_variant(variant: JsonVariant) -> List[int]:
    runs: Set[int] = set()
    for runs_file in variant.runs_files:
        runs.update(read_runs_file(runs_file))
    return sorted(runs)


def find_replay_file(paths: LtsepPaths, run: int) -> Optional[Path]:
    exact = paths.replay_source_dir / f"{paths.anatype}_coin_replay_production_{run}_-1.root"
    if exact.exists():
        return exact

    matches = sorted(paths.replay_source_dir.glob(f"*coin_replay_production_{run}_-1.root"))
    if matches:
        return matches[0]
    return None


def skim_files_for_run(paths: LtsepPaths, run: int) -> Tuple[Path, Path]:
    base = paths.skim_source_dir
    return (
        base / f"kaon_{run}_-1_Raw_Data.root",
        base / f"pion_{run}_-1_Raw_Data.root",
    )


def build_run_plans(paths: LtsepPaths, variants: Sequence[JsonVariant]) -> List[RunPlan]:
    plans: List[RunPlan] = []
    for variant in variants:
        for run in runs_for_variant(variant):
            replay_file = find_replay_file(paths, run)
            skim_kaon, skim_pion = skim_files_for_run(paths, run)

            if replay_file is None:
                plans.append(
                    RunPlan(
                        variant=variant,
                        run=run,
                        job_name=safe_name(f"{variant.family}_run{run}_applycuts"),
                        replay_file=None,
                        skim_kaon=skim_kaon,
                        skim_pion=skim_pion,
                        status="SKIP",
                        reason="missing_replay_root",
                    )
                )
                continue

            if skim_kaon.exists() and skim_pion.exists():
                plans.append(
                    RunPlan(
                        variant=variant,
                        run=run,
                        job_name=safe_name(f"{variant.family}_run{run}_applycuts"),
                        replay_file=replay_file,
                        skim_kaon=skim_kaon,
                        skim_pion=skim_pion,
                        status="SKIP",
                        reason="skim_outputs_already_exist",
                    )
                )
                continue

            plans.append(
                RunPlan(
                    variant=variant,
                    run=run,
                    job_name=safe_name(f"{variant.family}_run{run}_applycuts"),
                    replay_file=replay_file,
                    skim_kaon=skim_kaon,
                    skim_pion=skim_pion,
                    status="ADD",
                    reason="replay_ready",
                )
            )

    return plans


def workflow_name(args: argparse.Namespace, family_prefix: str) -> str:
    if args.workflow_name:
        return safe_name(args.workflow_name)
    user = os.environ.get("USER", "user")
    return safe_name(f"{args.workflow_prefix}_{family_prefix}_applycuts_{user}")


def run_command(cmd: Sequence[str], capture: bool = False) -> subprocess.CompletedProcess:
    return subprocess.run(
        list(cmd),
        text=True,
        capture_output=capture,
        check=False,
    )


def swif_status_exists(swif2_bin: str, workflow: str) -> bool:
    result = run_command([swif2_bin, "status", workflow, "-summary", "-display", "json"], capture=True)
    return result.returncode == 0


def ensure_workflow(swif2_bin: str, workflow: str, max_concurrent: int, submit: bool) -> Tuple[bool, List[str]]:
    if swif_status_exists(swif2_bin, workflow):
        return False, []

    cmd = [swif2_bin, "create", workflow, "-max-concurrent", str(max_concurrent)]
    if not submit:
        return True, cmd

    result = run_command(cmd, capture=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to create workflow {workflow}\n"
            f"STDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    return True, cmd


def get_existing_job_names(swif2_bin: str, workflow: str) -> Set[str]:
    result = run_command([swif2_bin, "status", workflow, "-jobs", "-display", "json"], capture=True)
    if result.returncode != 0:
        return set()
    try:
        payload = json.loads(result.stdout)
    except json.JSONDecodeError:
        return set()
    jobs = payload.get("jobs", [])
    names: Set[str] = set()
    if isinstance(jobs, list):
        for job in jobs:
            if isinstance(job, dict):
                name = job.get("job_name")
                if isinstance(name, str):
                    names.add(name)
    return names


def build_add_job_command(
    swif2_bin: str,
    workflow: str,
    account: str,
    partition: str,
    applycuts_script: str,
    plan: RunPlan,
    family_prefix: str,
    args: argparse.Namespace,
) -> List[str]:
    return [
        swif2_bin,
        "add-job",
        workflow,
        "-account",
        account,
        "-partition",
        partition,
        "-name",
        plan.job_name,
        "-cores",
        str(args.cores),
        "-ram",
        args.ram,
        "-disk",
        args.disk,
        "-time",
        args.time,
        "-tag",
        "run",
        str(plan.run),
        "-tag",
        "family",
        family_prefix,
        "-tag",
        "variant",
        plan.variant.family,
        "-tag",
        "stage",
        "applycuts",
        applycuts_script,
        plan.variant.epsilon,
        plan.variant.phi,
        plan.variant.q2,
        plan.variant.w,
        plan.variant.target,
        str(plan.run),
    ]


def print_summary(
    family_prefix: str,
    variants: Sequence[JsonVariant],
    plans: Sequence[RunPlan],
    workflow: str,
    args: argparse.Namespace,
    paths: LtsepPaths,
    existing_job_names: Set[str],
) -> None:
    print(f"Family prefix    : {family_prefix}")
    print(f"JSON directory   : {expand_path(args.json_dir)}")
    print(f"Workflow         : {workflow}")
    print(f"applyCuts script : {expand_path(args.applycuts_script)}")
    print(f"Replay source    : {paths.replay_source_dir}")
    print(f"Skim source      : {paths.skim_source_dir}")
    print(f"Account          : {args.account}")
    print(f"Partition        : {args.partition}")
    print(f"Variants matched : {len(variants)}")
    print(f"Runs inspected   : {len(plans)}")
    print()

    print("Matched JSON variants")
    print("-" * 80)
    for variant in variants:
        print(f"* {variant.path.name}")
        print(f"    phi={variant.phi} epsilon={variant.epsilon} target={variant.target}")
        for runs_file in variant.runs_files:
            print(f"    runs_file: {runs_file}")
    print()

    print("Planned SWIF2 jobs")
    print("-" * 80)
    applycuts_script = str(expand_path(args.applycuts_script))
    for plan in plans:
        if plan.should_add:
            if args.skip_existing and plan.job_name in existing_job_names:
                status = "SKIP existing"
            else:
                status = "ADD"
        else:
            status = f"SKIP {plan.reason}"
        print(f"[{status}] run={plan.run} variant={plan.variant.family} job={plan.job_name}")
        if plan.replay_file is None:
            print(f"         replay  : missing under {paths.replay_source_dir}")
        else:
            print(f"         replay  : {plan.replay_file}")
        print(f"         skim    : {plan.skim_kaon.name}, {plan.skim_pion.name}")
        if plan.should_add:
            cmd = build_add_job_command(
                args.swif2_bin,
                workflow,
                args.account,
                args.partition,
                applycuts_script,
                plan,
                family_prefix,
                args,
            )
            print("         cmd     : " + " ".join(shlex.quote(x) for x in cmd))
    print()


def submit_jobs(
    workflow: str,
    family_prefix: str,
    plans: Sequence[RunPlan],
    existing_job_names: Set[str],
    args: argparse.Namespace,
) -> int:
    created, create_cmd = ensure_workflow(args.swif2_bin, workflow, args.max_concurrent, submit=True)
    if created:
        print("Created workflow:")
        print("  " + " ".join(shlex.quote(x) for x in create_cmd))
        print()

    applycuts_script = str(expand_path(args.applycuts_script))
    rc = 0
    for plan in plans:
        if not plan.should_add:
            continue
        if args.skip_existing and plan.job_name in existing_job_names:
            print(f"Skipping existing job: {plan.job_name}")
            continue
        cmd = build_add_job_command(
            args.swif2_bin,
            workflow,
            args.account,
            args.partition,
            applycuts_script,
            plan,
            family_prefix,
            args,
        )
        print("Adding job:")
        print("  " + " ".join(shlex.quote(x) for x in cmd))
        result = run_command(cmd, capture=True)
        if result.stdout:
            print(result.stdout.rstrip())
        if result.stderr:
            print(result.stderr.rstrip(), file=sys.stderr)
        if result.returncode != 0:
            rc = result.returncode
            print(f"ERROR: add-job failed for {plan.job_name}", file=sys.stderr)
            break
        print()

    if rc != 0:
        return rc

    if not args.no_run:
        cmd = [args.swif2_bin, "run", workflow]
        print("Starting workflow:")
        print("  " + " ".join(shlex.quote(x) for x in cmd))
        result = run_command(cmd, capture=True)
        if result.stdout:
            print(result.stdout.rstrip())
        if result.stderr:
            print(result.stderr.rstrip(), file=sys.stderr)
        rc = result.returncode

    return rc


def validate_applycuts_script(path_text: str) -> None:
    path = expand_path(path_text)
    if not path.exists():
        raise FileNotFoundError(f"applyCuts script does not exist: {path}")
    if not os.access(path, os.X_OK):
        raise PermissionError(f"applyCuts script is not executable: {path}")


def main() -> int:
    args = parse_args()
    family_prefix = normalize_family(args.q2, args.w)
    json_dir = expand_path(args.json_dir)
    validate_applycuts_script(args.applycuts_script)
    paths = resolve_ltsep_paths(__file__)

    variants = discover_json_variants(json_dir, family_prefix, args.family_regex)
    if not variants:
        raise SystemExit(f"No JSON variants found for family prefix {family_prefix} in {json_dir}")

    plans = build_run_plans(paths, variants)
    if not plans:
        raise SystemExit("No runs were discovered from the matched JSON files")

    workflow = workflow_name(args, family_prefix)
    existing_job_names = get_existing_job_names(args.swif2_bin, workflow) if swif_status_exists(args.swif2_bin, workflow) else set()

    print_summary(
        family_prefix,
        variants,
        plans,
        workflow,
        args,
        paths,
        existing_job_names,
    )

    print("Monitoring commands")
    print("-" * 80)
    print(f"{args.swif2_bin} status {workflow} -jobs -display json")
    print(f"{args.swif2_bin} status {workflow} -summary -display json")
    print(f"{args.swif2_bin} show-job {workflow} <job-name>")
    print()

    if not args.submit:
        print("Dry run only. Use --submit to create/update the workflow and add jobs.")
        return 0

    return submit_jobs(workflow, family_prefix, plans, existing_job_names, args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
