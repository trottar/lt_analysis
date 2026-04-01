#!/usr/bin/env python3
"""
Submit KaonLT replay jobs to JLab SWIF2/ifarm.

Behavior:
- User supplies only q2 and w, e.g. `5p5 3p02`.
- The script discovers all matching JSON files under the KaonLT input directory
  whose names begin with `Q{q2}W{w}`.
- Each JSON may contain one or more jobs[*].runs_file entries.
- All runs from all matched JSON files are merged into one unique run set.
- Missing raw cache files are requested for staging before replay submission.
- Exactly one SWIF2 job is submitted per unique run number.
- Each job runs:
    replay_env/run_replay.sh <run>
- Resource requests are chosen per run. If the user supplies one or more
  --size-glob-template values, the script estimates the input size for a run by
  summing matching file sizes and applies a tiered heuristic. Otherwise it falls
  back to conservative Hall-C replay defaults.

Dry-run is the default. Use --submit to create/modify the workflow and add jobs.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_MANIFEST_DIR = str(REPO_ROOT / "input" / "kaon")
DEFAULT_REPLAY_SCRIPT = str(REPO_ROOT / "replay_env" / "run_replay.sh")
DEFAULT_JASMINE_SCRIPT = str(REPO_ROOT / "farm_env" / "jasmine_put_from_manifest.py")
DEFAULT_WORKFLOW_PREFIX = "kaonlt"
DEFAULT_ACCOUNT = "hallc"
DEFAULT_PARTITION = "production"
DEFAULT_MAX_CONCURRENT = 200
DEFAULT_CORES = 1
DEFAULT_FALLBACK_RAM = "8g"
DEFAULT_FALLBACK_DISK = "40g"
DEFAULT_FALLBACK_TIME = "8h"
DEFAULT_SWIF2 = os.environ.get("SWIF2_BIN", "swif2")
DEFAULT_PYTHON_BIN = os.environ.get("PYTHON_BIN", "python3")
DEFAULT_JASMINE_RAM = "4g"
DEFAULT_JASMINE_DISK = "20g"
DEFAULT_JASMINE_TIME = "12h"
DEFAULT_JASMINE_STAGE_ROOT = "/scratch/$USER/jasmine_stage"
DEFAULT_RAW_CACHE_GLOB_TEMPLATE = "/cache/hallc/raw/coin_all_{run5}.dat"
DEFAULT_CACHE_REQUEST_TEMPLATE = "jcache get {cache_file}"

RUN_LINE_RE = re.compile(r"^\s*(\d+)\s*$")
SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9_.-]+")


@dataclass(frozen=True)
class JsonVariant:
    path: Path
    family: str
    runs_files: Tuple[Path, ...]


@dataclass(frozen=True)
class ResourceRequest:
    ram: str
    disk: str
    time: str
    cores: int
    reason: str
    estimated_input_bytes: Optional[int]


@dataclass(frozen=True)
class RunPlan:
    run: int
    job_name: str
    variants: Tuple[str, ...]
    resources: ResourceRequest
    representative_manifest: Path
    cache_file: Optional[Path]
    cache_ready: bool
    cache_reason: str
    cache_request_command: Tuple[str, ...] = ()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Submit one KaonLT replay SWIF2 job per unique run for a Q/W family."
    )
    parser.add_argument("q2", help="Q value without prefix, e.g. 5p5")
    parser.add_argument("w", help="W value without prefix, e.g. 3p02")
    parser.add_argument(
        "--manifest-dir",
        "--json-dir",
        dest="manifest_dir",
        default=DEFAULT_MANIFEST_DIR,
        help=f"Directory containing KaonLT manifest JSON files (default: {DEFAULT_MANIFEST_DIR})",
    )
    parser.add_argument(
        "--replay-script",
        default=DEFAULT_REPLAY_SCRIPT,
        help=f"Replay wrapper to execute per run (default: {DEFAULT_REPLAY_SCRIPT})",
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
        "--account",
        default=DEFAULT_ACCOUNT,
        help=f"SWIF2 account / Slurm account (default: {DEFAULT_ACCOUNT})",
    )
    parser.add_argument(
        "--partition",
        default=DEFAULT_PARTITION,
        help=(
            "SWIF2 partition / Slurm partition (default: production). "
            "Use priority for small debug tests."
        ),
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
        "--fallback-ram",
        default=DEFAULT_FALLBACK_RAM,
        help=f"RAM to request when input size is unknown (default: {DEFAULT_FALLBACK_RAM})",
    )
    parser.add_argument(
        "--fallback-disk",
        default=DEFAULT_FALLBACK_DISK,
        help=f"Disk to request when input size is unknown (default: {DEFAULT_FALLBACK_DISK})",
    )
    parser.add_argument(
        "--fallback-time",
        default=DEFAULT_FALLBACK_TIME,
        help=f"Wall time to request when input size is unknown (default: {DEFAULT_FALLBACK_TIME})",
    )
    parser.add_argument(
        "--size-glob-template",
        action="append",
        default=[],
        help=(
            "Optional glob template used to estimate per-run input size. Repeatable. "
            "Use {run} for the run number and shell wildcards if needed. Example: "
            "'/cache/hallc/raw/*_{run}_*.dat'"
        ),
    )
    parser.add_argument(
        "--raw-cache-glob-template",
        action="append",
        default=[],
        help=(
            "Optional glob template used to locate the replay raw data in cache. "
            "Repeatable. Supports {run} and {run5}. "
            f"Default: {DEFAULT_RAW_CACHE_GLOB_TEMPLATE}"
        ),
    )
    parser.add_argument(
        "--cache-request-template",
        default=DEFAULT_CACHE_REQUEST_TEMPLATE,
        help=(
            "Command template used to request a missing raw cache file. "
            "Supports {run}, {run5}, and {cache_file}. "
            f"Default: {DEFAULT_CACHE_REQUEST_TEMPLATE}"
        ),
    )
    parser.add_argument(
        "--no-cache-request",
        dest="request_missing_cache",
        action="store_false",
        default=True,
        help="Do not request raw-file staging when the cache file is missing.",
    )
    parser.add_argument(
        "--workflow-name",
        default=None,
        help="Override the auto-generated workflow name.",
    )
    parser.add_argument(
        "--family-regex",
        default=None,
        help="Optional regex applied to matching JSON basenames for additional filtering.",
    )
    parser.add_argument(
        "--allow-empty-size-match",
        action="store_true",
        help="Do not warn when a size-glob template matches nothing for a run.",
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
    parser.add_argument(
        "--no-auto-jasmine",
        dest="auto_jasmine",
        action="store_false",
        default=True,
        help="Do not add dependent Jasmine tape-upload jobs.",
    )
    parser.add_argument(
        "--jasmine-script",
        default=DEFAULT_JASMINE_SCRIPT,
        help=f"Jasmine helper script to execute per run (default: {DEFAULT_JASMINE_SCRIPT})",
    )
    parser.add_argument(
        "--python-bin",
        default=DEFAULT_PYTHON_BIN,
        help=f"Python executable used for Jasmine jobs (default: {DEFAULT_PYTHON_BIN})",
    )
    parser.add_argument(
        "--jasmine-ram",
        default=DEFAULT_JASMINE_RAM,
        help=f"RAM per Jasmine upload job (default: {DEFAULT_JASMINE_RAM})",
    )
    parser.add_argument(
        "--jasmine-disk",
        default=DEFAULT_JASMINE_DISK,
        help=f"Disk per Jasmine upload job (default: {DEFAULT_JASMINE_DISK})",
    )
    parser.add_argument(
        "--jasmine-time",
        default=DEFAULT_JASMINE_TIME,
        help=f"Wall time per Jasmine upload job (default: {DEFAULT_JASMINE_TIME})",
    )
    parser.add_argument(
        "--jasmine-stage-root",
        default=DEFAULT_JASMINE_STAGE_ROOT,
        help=f"Stage root passed to Jasmine jobs (default: {DEFAULT_JASMINE_STAGE_ROOT})",
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


def render_template(template: str, run: int, extra: Optional[Dict[str, str]] = None) -> str:
    values = {
        "run": str(run),
        "run5": f"{int(run):05d}",
    }
    if extra:
        values.update(extra)
    return os.path.expandvars(template).format(**values)


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
        data = load_json(path)
        runs_files = extract_runs_files(data)
        if not runs_files:
            continue
        variants.append(
            JsonVariant(
                path=path,
                family=path.stem,
                runs_files=tuple(runs_files),
            )
        )

    return variants


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



def collect_runs(variants: Sequence[JsonVariant]) -> Tuple[Dict[int, Set[str]], Dict[int, Path]]:
    runs_to_variants: Dict[int, Set[str]] = {}
    runs_to_manifest: Dict[int, Path] = {}
    for variant in variants:
        for runs_file in variant.runs_files:
            runs = read_runs_file(runs_file)
            for run in runs:
                runs_to_variants.setdefault(run, set()).add(variant.family)
                if run not in runs_to_manifest:
                    runs_to_manifest[run] = variant.path
    return runs_to_variants, runs_to_manifest



def estimate_run_input_bytes(
    run: int,
    templates: Sequence[str],
    allow_empty: bool,
) -> Tuple[Optional[int], List[str]]:
    if not templates:
        return None, []

    total_bytes = 0
    matched_paths: Set[Path] = set()
    warnings: List[str] = []

    for template in templates:
        pattern = render_template(template, run)
        for match in glob.glob(os.path.expanduser(pattern)):
            path = Path(match)
            if path.exists() and path.is_file():
                matched_paths.add(path.resolve())

        if not allow_empty and not glob.glob(os.path.expanduser(pattern)):
            warnings.append(f"Run {run}: size-glob template matched nothing: {pattern}")

    for path in matched_paths:
        try:
            total_bytes += path.stat().st_size
        except FileNotFoundError:
            continue

    if not matched_paths:
        return None, warnings
    return total_bytes, warnings


def raw_cache_templates(args: argparse.Namespace) -> List[str]:
    templates = list(args.raw_cache_glob_template)
    if templates:
        return templates
    return [DEFAULT_RAW_CACHE_GLOB_TEMPLATE]


def find_cached_raw_file(run: int, templates: Sequence[str]) -> Optional[Path]:
    matches: List[Path] = []
    for template in templates:
        pattern = render_template(template, run)
        for match in glob.glob(os.path.expanduser(pattern)):
            path = Path(match)
            if path.exists() and path.is_file():
                matches.append(path.resolve())
    if not matches:
        return None
    matches = sorted(set(matches))
    return matches[0]


def build_cache_request_command(args: argparse.Namespace, run: int, cache_file: Path) -> List[str]:
    command_text = render_template(
        args.cache_request_template,
        run,
        extra={"cache_file": str(cache_file)},
    )
    return shlex.split(command_text)



def choose_resources(
    run: int,
    estimated_input_bytes: Optional[int],
    cores: int,
    fallback_ram: str,
    fallback_disk: str,
    fallback_time: str,
) -> ResourceRequest:
    if estimated_input_bytes is None:
        return ResourceRequest(
            ram=fallback_ram,
            disk=fallback_disk,
            time=fallback_time,
            cores=cores,
            reason="fallback_no_size_estimate",
            estimated_input_bytes=None,
        )

    size_gb = estimated_input_bytes / float(1024 ** 3)

    # Hall-C replay heuristic:
    # - small runs get lighter requests
    # - larger runs scale wall time and disk upward
    if size_gb < 5.0:
        ram, disk, time_label, tier = "4g", "15g", "3h", "tier_lt5GB"
    elif size_gb < 15.0:
        ram, disk, time_label, tier = "6g", "25g", "5h", "tier_5to15GB"
    elif size_gb < 30.0:
        ram, disk, time_label, tier = "8g", "40g", "8h", "tier_15to30GB"
    else:
        ram, disk, time_label, tier = "12g", "60g", "12h", "tier_gt30GB"

    return ResourceRequest(
        ram=ram,
        disk=disk,
        time=time_label,
        cores=cores,
        reason=tier,
        estimated_input_bytes=estimated_input_bytes,
    )



def build_run_plans(
    family_prefix: str,
    runs_to_variants: Dict[int, Set[str]],
    runs_to_manifest: Dict[int, Path],
    args: argparse.Namespace,
) -> Tuple[List[RunPlan], List[str]]:
    warnings: List[str] = []
    plans: List[RunPlan] = []

    for run in sorted(runs_to_variants):
        cache_file = find_cached_raw_file(run, raw_cache_templates(args))
        cache_ready = cache_file is not None
        cache_reason = "cache_ready" if cache_ready else "missing_cache"
        cache_request_command: Tuple[str, ...] = ()
        if not cache_ready:
            nominal_cache_file = Path(
                render_template(raw_cache_templates(args)[0], run)
            ).expanduser()
            cache_file = nominal_cache_file
            if args.request_missing_cache:
                cache_request_command = tuple(build_cache_request_command(args, run, cache_file))

        est_bytes, size_warnings = estimate_run_input_bytes(
            run,
            args.size_glob_template,
            args.allow_empty_size_match,
        )
        warnings.extend(size_warnings)
        resources = choose_resources(
            run,
            est_bytes,
            args.cores,
            args.fallback_ram,
            args.fallback_disk,
            args.fallback_time,
        )
        job_name = safe_name(f"{family_prefix}_run{run}")
        plans.append(
            RunPlan(
                run=run,
                job_name=job_name,
                variants=tuple(sorted(runs_to_variants[run])),
                resources=resources,
                representative_manifest=runs_to_manifest[run],
                cache_file=cache_file,
                cache_ready=cache_ready,
                cache_reason=cache_reason,
                cache_request_command=cache_request_command,
            )
        )

    return plans, warnings



def workflow_name(args: argparse.Namespace, family_prefix: str) -> str:
    if args.workflow_name:
        return safe_name(args.workflow_name)
    user = os.environ.get("USER", "user")
    return safe_name(f"{args.workflow_prefix}_{family_prefix}_{user}")



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
    replay_script: str,
    plan: RunPlan,
    family_prefix: str,
) -> List[str]:
    cmd = [
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
        str(plan.resources.cores),
        "-ram",
        plan.resources.ram,
        "-disk",
        plan.resources.disk,
        "-time",
        plan.resources.time,
        "-tag",
        "run",
        str(plan.run),
        "-tag",
        "family",
        family_prefix,
        "-tag",
        "variant_count",
        str(len(plan.variants)),
        replay_script,
        str(plan.run),
    ]
    return cmd


def jasmine_job_name(plan: RunPlan) -> str:
    return safe_name(f"{plan.job_name}_jasmine")


def jasmine_stage_dir(args: argparse.Namespace, plan: RunPlan) -> str:
    return str(Path(os.path.expandvars(args.jasmine_stage_root)).expanduser() / jasmine_job_name(plan))


def build_jasmine_add_job_command(
    swif2_bin: str,
    workflow: str,
    account: str,
    partition: str,
    plan: RunPlan,
    family_prefix: str,
    args: argparse.Namespace,
) -> List[str]:
    jasmine_script = str(expand_path(args.jasmine_script))
    return [
        swif2_bin,
        "add-job",
        workflow,
        "-account",
        account,
        "-partition",
        partition,
        "-name",
        jasmine_job_name(plan),
        "-cores",
        "1",
        "-ram",
        args.jasmine_ram,
        "-disk",
        args.jasmine_disk,
        "-time",
        args.jasmine_time,
        "-tag",
        "run",
        str(plan.run),
        "-tag",
        "family",
        family_prefix,
        "-tag",
        "stage",
        "jasmine_replay",
        "-antecedent",
        plan.job_name,
        args.python_bin,
        jasmine_script,
        "--manifest-path",
        str(plan.representative_manifest),
        "--run",
        str(plan.run),
        "--product-kind",
        "replay",
        "--stage-dir",
        jasmine_stage_dir(args, plan),
        "--submit",
    ]


def format_bytes(num_bytes: Optional[int]) -> str:
    if num_bytes is None:
        return "unknown"
    value = float(num_bytes)
    for unit in ["B", "KiB", "MiB", "GiB", "TiB"]:
        if value < 1024.0 or unit == "TiB":
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return str(num_bytes)



def print_summary(
    family_prefix: str,
    variants: Sequence[JsonVariant],
    runs_to_variants: Dict[int, Set[str]],
    plans: Sequence[RunPlan],
    workflow: str,
    args: argparse.Namespace,
    size_warnings: Sequence[str],
    existing_job_names: Set[str],
) -> None:
    print(f"Family prefix   : {family_prefix}")
    print(f"Manifest directory: {expand_path(args.manifest_dir)}")
    print(f"Workflow         : {workflow}")
    print(f"Replay script    : {expand_path(args.replay_script)}")
    print(f"Auto Jasmine     : {'yes' if args.auto_jasmine else 'no'}")
    if args.auto_jasmine:
        print(f"Jasmine script   : {expand_path(args.jasmine_script)}")
    print(f"Account          : {args.account}")
    print(f"Partition        : {args.partition}")
    print(f"Runs discovered  : {len(runs_to_variants)} unique")
    print(f"Variants matched : {len(variants)}")
    print()

    print("Matched JSON variants")
    print("-" * 80)
    for variant in variants:
        print(f"* {variant.path.name}")
        for runs_file in variant.runs_files:
            print(f"    runs_file: {runs_file}")
    print()

    duplicate_runs = {run: fams for run, fams in runs_to_variants.items() if len(fams) > 1}
    if duplicate_runs:
        print("Runs present in multiple JSON variants")
        print("-" * 80)
        for run in sorted(duplicate_runs):
            print(f"* run {run}: {', '.join(sorted(duplicate_runs[run]))}")
        print()

    if size_warnings:
        print("Size-estimation warnings")
        print("-" * 80)
        for msg in size_warnings:
            print(f"* {msg}")
        print()

    print("Planned SWIF2 jobs")
    print("-" * 80)
    for plan in plans:
        if not plan.cache_ready:
            print(
                f"[WAIT {plan.cache_reason}] run={plan.run} job={plan.job_name} "
                f"cache={plan.cache_file}"
            )
            if plan.cache_request_command:
                print("         request : " + " ".join(shlex.quote(x) for x in plan.cache_request_command))
            else:
                print("         request : disabled")
            continue
        status = "SKIP existing" if plan.job_name in existing_job_names and args.skip_existing else "ADD"
        print(
            f"[{status}] run={plan.run} job={plan.job_name} "
            f"ram={plan.resources.ram} disk={plan.resources.disk} time={plan.resources.time} "
            f"cores={plan.resources.cores} size={format_bytes(plan.resources.estimated_input_bytes)} "
            f"heuristic={plan.resources.reason}"
        )
        print(f"         variants: {', '.join(plan.variants)}")
        cmd = build_add_job_command(
            args.swif2_bin,
            workflow,
            args.account,
            args.partition,
            str(expand_path(args.replay_script)),
            plan,
            family_prefix,
        )
        print("         cmd     : " + " ".join(shlex.quote(x) for x in cmd))
        if args.auto_jasmine:
            upload_job = jasmine_job_name(plan)
            upload_status = "SKIP existing" if upload_job in existing_job_names and args.skip_existing else "ADD"
            upload_cmd = build_jasmine_add_job_command(
                args.swif2_bin,
                workflow,
                args.account,
                args.partition,
                plan,
                family_prefix,
                args,
            )
            print(f"[{upload_status}] upload job={upload_job} antecedent={plan.job_name}")
            print(f"         manifest: {plan.representative_manifest}")
            print("         cmd     : " + " ".join(shlex.quote(x) for x in upload_cmd))
    print()



def submit_jobs(
    workflow: str,
    family_prefix: str,
    plans: Sequence[RunPlan],
    existing_job_names: Set[str],
    args: argparse.Namespace,
) -> int:
    replay_script = str(expand_path(args.replay_script))
    existing_names = set(existing_job_names)
    rc = 0
    eligible_plans = [plan for plan in plans if plan.cache_ready]

    missing_cache_plans = [plan for plan in plans if not plan.cache_ready]
    for plan in missing_cache_plans:
        if not plan.cache_request_command:
            print(f"Cache missing for run {plan.run}; no request command configured.")
            continue
        print(f"Requesting cache staging for run {plan.run}:")
        print("  " + " ".join(shlex.quote(x) for x in plan.cache_request_command))
        result = run_command(plan.cache_request_command, capture=True)
        if result.stdout:
            print(result.stdout.rstrip())
        if result.stderr:
            print(result.stderr.rstrip(), file=sys.stderr)
        if result.returncode != 0:
            rc = result.returncode
            print(f"ERROR: cache request failed for run {plan.run}", file=sys.stderr)
            return rc
        print()

    if not eligible_plans:
        print("No cache-ready replay jobs to submit yet.")
        print("Missing raw files were requested for staging where configured.")
        return rc

    created, create_cmd = ensure_workflow(args.swif2_bin, workflow, args.max_concurrent, submit=True)
    if created:
        print("Created workflow:")
        print("  " + " ".join(shlex.quote(x) for x in create_cmd))
        print()

    for plan in plans:
        if not plan.cache_ready:
            continue
        replay_exists = plan.job_name in existing_names
        if args.skip_existing and replay_exists:
            print(f"Skipping existing job: {plan.job_name}")
        else:
            cmd = build_add_job_command(
                args.swif2_bin,
                workflow,
                args.account,
                args.partition,
                replay_script,
                plan,
                family_prefix,
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
            existing_names.add(plan.job_name)
            print()

        if not args.auto_jasmine:
            continue

        upload_job = jasmine_job_name(plan)
        if args.skip_existing and upload_job in existing_names:
            print(f"Skipping existing job: {upload_job}")
            continue

        upload_cmd = build_jasmine_add_job_command(
            args.swif2_bin,
            workflow,
            args.account,
            args.partition,
            plan,
            family_prefix,
            args,
        )
        print("Adding job:")
        print("  " + " ".join(shlex.quote(x) for x in upload_cmd))
        result = run_command(upload_cmd, capture=True)
        if result.stdout:
            print(result.stdout.rstrip())
        if result.stderr:
            print(result.stderr.rstrip(), file=sys.stderr)
        if result.returncode != 0:
            rc = result.returncode
            print(f"ERROR: add-job failed for {upload_job}", file=sys.stderr)
            break
        existing_names.add(upload_job)
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



def validate_replay_script(path_text: str) -> None:
    path = expand_path(path_text)
    if not path.exists():
        raise FileNotFoundError(f"Replay script does not exist: {path}")
    if not os.access(path, os.X_OK):
        raise PermissionError(f"Replay script is not executable: {path}")


def validate_jasmine_script(path_text: str) -> None:
    path = expand_path(path_text)
    if not path.exists():
        raise FileNotFoundError(f"Jasmine script does not exist: {path}")



def main() -> int:
    args = parse_args()
    family_prefix = normalize_family(args.q2, args.w)
    json_dir = expand_path(args.manifest_dir)
    validate_replay_script(args.replay_script)
    if args.auto_jasmine:
        validate_jasmine_script(args.jasmine_script)

    variants = discover_json_variants(json_dir, family_prefix, args.family_regex)
    if not variants:
        raise SystemExit(f"No JSON variants found for family prefix {family_prefix} in {json_dir}")

    runs_to_variants, runs_to_manifest = collect_runs(variants)
    if not runs_to_variants:
        raise SystemExit("No runs were discovered from the matched JSON files")

    plans, size_warnings = build_run_plans(family_prefix, runs_to_variants, runs_to_manifest, args)
    workflow = workflow_name(args, family_prefix)

    existing_job_names = get_existing_job_names(args.swif2_bin, workflow) if swif_status_exists(args.swif2_bin, workflow) else set()

    print_summary(
        family_prefix,
        variants,
        runs_to_variants,
        plans,
        workflow,
        args,
        size_warnings,
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
