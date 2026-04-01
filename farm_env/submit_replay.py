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
- Replay output is registered with SWIF `-output` and transferred to MSS after
  the job finishes.

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
DEFAULT_SOFTENV_WRAPPER = str(REPO_ROOT / "farm_env" / "run_with_softenv.sh")
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
DEFAULT_RAW_CACHE_GLOB_TEMPLATE = "/cache/hallc/spring17/raw/coin_all_{run5}.dat"
DEFAULT_RAW_MSS_TEMPLATE = "/mss/hallc/spring17/raw/coin_all_{run5}.dat"
DEFAULT_CACHE_REQUEST_TEMPLATE = "jcache get {mss_file}"

RUN_LINE_RE = re.compile(r"^\s*(\d+)\s*$")
SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9_.-]+")


@dataclass(frozen=True)
class JsonVariant:
    path: Path
    family: str
    partition: str
    runs_files: Tuple[Path, ...]
    replay_destinations: Tuple[Path, ...]
    report_destinations: Tuple[Path, ...]
    archive_prefix_templates: Tuple[str, ...]


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
    replay_destination: Path
    report_destination: Optional[Path]
    report_tarball_basename: Optional[str]
    cache_file: Optional[Path]
    cache_ready: bool
    cache_reason: str
    replay_output_exists: bool = False
    report_output_exists: bool = False
    tarball_output_exists: bool = False
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
        default=None,
        help=(
            "SWIF2 partition / Slurm partition. Defaults to defaults.partition "
            f"from the manifest, or {DEFAULT_PARTITION} when unspecified."
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
            "'/cache/hallc/spring17/raw/*_{run}_*.dat'"
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
            "Supports {run}, {run5}, {cache_file}, and {mss_file}. "
            f"Default: {DEFAULT_CACHE_REQUEST_TEMPLATE}"
        ),
    )
    parser.add_argument(
        "--raw-mss-template",
        default=DEFAULT_RAW_MSS_TEMPLATE,
        help=(
            "Template for the MSS-backed raw file used when requesting cache staging. "
            "Supports {run} and {run5}. "
            f"Default: {DEFAULT_RAW_MSS_TEMPLATE}"
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
        default=False,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--jasmine-script",
        default=DEFAULT_JASMINE_SCRIPT,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--python-bin",
        default=DEFAULT_PYTHON_BIN,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--jasmine-ram",
        default=DEFAULT_JASMINE_RAM,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--jasmine-disk",
        default=DEFAULT_JASMINE_DISK,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--jasmine-time",
        default=DEFAULT_JASMINE_TIME,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--jasmine-stage-root",
        default=DEFAULT_JASMINE_STAGE_ROOT,
        help=argparse.SUPPRESS,
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


def replay_output_basename(run: int) -> str:
    return f"Kaon_coin_replay_production_{run}_-1.root"


def replay_report_basename(run: int) -> str:
    return f"Kaon_output_coin_production_Summary_{run}_-1.report"


def default_report_tarball_basename(run: int) -> str:
    return f"FullReplay_run_{run}_KaonLT.tar"


def replay_mss_output_file(plan: RunPlan) -> Path:
    return plan.replay_destination / replay_output_basename(plan.run)


def replay_report_mss_output_file(plan: RunPlan) -> Optional[Path]:
    if plan.report_destination is None:
        return None
    return plan.report_destination / replay_report_basename(plan.run)


def replay_report_tarball_mss_output_file(plan: RunPlan) -> Optional[Path]:
    if plan.report_destination is None or plan.report_tarball_basename is None:
        return None
    return plan.report_destination / plan.report_tarball_basename


def plan_has_missing_mss_output(plan: RunPlan) -> bool:
    report_missing = plan.report_destination is not None and not plan.report_output_exists
    tarball_missing = plan.report_tarball_basename is not None and not plan.tarball_output_exists
    root_missing = not plan.replay_output_exists
    return root_missing or report_missing or tarball_missing


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
        replay_destinations = extract_replay_destinations(data)
        report_destinations = extract_report_destinations(data)
        archive_prefix_templates = extract_archive_prefix_templates(data)
        partition = extract_partition(data)
        if not runs_files or not replay_destinations:
            continue
        variants.append(
            JsonVariant(
                path=path.resolve(),
                family=path.stem,
                partition=partition,
                runs_files=tuple(runs_files),
                replay_destinations=tuple(replay_destinations),
                report_destinations=tuple(report_destinations),
                archive_prefix_templates=tuple(archive_prefix_templates),
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


def extract_replay_destinations(data: Dict) -> List[Path]:
    jobs = data.get("jobs", [])
    destinations: List[Path] = []
    if not isinstance(jobs, list):
        return destinations
    for entry in jobs:
        if not isinstance(entry, dict):
            continue
        destination = entry.get("destination")
        if isinstance(destination, str) and destination.strip():
            destinations.append(expand_path(destination))
    return destinations


def extract_report_destinations(data: Dict) -> List[Path]:
    jobs = data.get("jobs", [])
    destinations: List[Path] = []
    if not isinstance(jobs, list):
        return destinations
    for entry in jobs:
        if not isinstance(entry, dict):
            continue
        destination = entry.get("tar_destination")
        if isinstance(destination, str) and destination.strip():
            destinations.append(expand_path(destination))
    return destinations


def extract_archive_prefix_templates(data: Dict) -> List[str]:
    jobs = data.get("jobs", [])
    templates: List[str] = []
    if not isinstance(jobs, list):
        return templates
    for entry in jobs:
        if not isinstance(entry, dict):
            continue
        template = entry.get("archive_prefix_template", entry.get("archive_prefix"))
        if isinstance(template, str) and template.strip():
            templates.append(template.strip())
    return templates


def extract_partition(data: Dict) -> str:
    defaults = data.get("defaults", {})
    if not isinstance(defaults, dict):
        return DEFAULT_PARTITION
    partition = str(defaults.get("partition", DEFAULT_PARTITION)).strip()
    return partition or DEFAULT_PARTITION



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



def collect_runs(
    variants: Sequence[JsonVariant],
) -> Tuple[
    Dict[int, Set[str]],
    Dict[int, Path],
    Dict[int, Set[Path]],
    Dict[int, Set[Path]],
    Dict[int, Set[str]],
]:
    runs_to_variants: Dict[int, Set[str]] = {}
    runs_to_manifest: Dict[int, Path] = {}
    runs_to_destinations: Dict[int, Set[Path]] = {}
    runs_to_report_destinations: Dict[int, Set[Path]] = {}
    runs_to_archive_prefix_templates: Dict[int, Set[str]] = {}
    for variant in variants:
        for runs_file in variant.runs_files:
            runs = read_runs_file(runs_file)
            for run in runs:
                runs_to_variants.setdefault(run, set()).add(variant.family)
                if run not in runs_to_manifest:
                    runs_to_manifest[run] = variant.path
                runs_to_destinations.setdefault(run, set()).update(variant.replay_destinations)
                if variant.report_destinations:
                    runs_to_report_destinations.setdefault(run, set()).update(variant.report_destinations)
                if variant.archive_prefix_templates:
                    runs_to_archive_prefix_templates.setdefault(run, set()).update(variant.archive_prefix_templates)
    return (
        runs_to_variants,
        runs_to_manifest,
        runs_to_destinations,
        runs_to_report_destinations,
        runs_to_archive_prefix_templates,
    )


def resolve_run_destinations(runs_to_destinations: Dict[int, Set[Path]]) -> Dict[int, Path]:
    resolved: Dict[int, Path] = {}
    conflicts: List[str] = []
    for run, destinations in sorted(runs_to_destinations.items()):
        unique_destinations = sorted(set(destinations))
        if len(unique_destinations) != 1:
            conflicts.append(
                f"run {run}: " + ", ".join(str(destination) for destination in unique_destinations)
            )
            continue
        resolved[run] = unique_destinations[0]
    if conflicts:
        joined = "\n".join(conflicts)
        raise ValueError(
            "Replay manifests disagree on the MSS destination for one or more runs.\n"
            f"{joined}"
        )
    return resolved


def resolve_optional_run_destinations(runs_to_destinations: Dict[int, Set[Path]]) -> Dict[int, Optional[Path]]:
    resolved: Dict[int, Optional[Path]] = {}
    conflicts: List[str] = []
    for run, destinations in sorted(runs_to_destinations.items()):
        unique_destinations = sorted(set(destinations))
        if not unique_destinations:
            resolved[run] = None
            continue
        if len(unique_destinations) != 1:
            conflicts.append(
                f"run {run}: " + ", ".join(str(destination) for destination in unique_destinations)
            )
            continue
        resolved[run] = unique_destinations[0]
    if conflicts:
        joined = "\n".join(conflicts)
        raise ValueError(
            "Replay manifests disagree on the REPORT_OUTPUT destination for one or more runs.\n"
            f"{joined}"
        )
    return resolved


def resolve_optional_run_templates(runs_to_templates: Dict[int, Set[str]]) -> Dict[int, Optional[str]]:
    resolved: Dict[int, Optional[str]] = {}
    conflicts: List[str] = []
    for run, templates in sorted(runs_to_templates.items()):
        unique_templates = sorted(set(templates))
        if not unique_templates:
            resolved[run] = None
            continue
        if len(unique_templates) != 1:
            conflicts.append(f"run {run}: " + ", ".join(unique_templates))
            continue
        resolved[run] = unique_templates[0]
    if conflicts:
        joined = "\n".join(conflicts)
        raise ValueError(
            "Replay manifests disagree on the archive_prefix_template for one or more runs.\n"
            f"{joined}"
        )
    return resolved


def resolve_family_partition(variants: Sequence[JsonVariant]) -> str:
    partitions = sorted({variant.partition for variant in variants})
    if not partitions:
        return DEFAULT_PARTITION
    if len(partitions) != 1:
        raise ValueError(
            "Matched replay manifests disagree on defaults.partition: "
            + ", ".join(partitions)
        )
    return partitions[0]



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
    mss_file = render_template(args.raw_mss_template, run)
    command_text = render_template(
        args.cache_request_template,
        run,
        extra={"cache_file": str(cache_file), "mss_file": mss_file},
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
    runs_to_destinations: Dict[int, Path],
    runs_to_report_destinations: Dict[int, Optional[Path]],
    runs_to_archive_prefix_templates: Dict[int, Optional[str]],
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
        replay_destination = runs_to_destinations[run]
        report_destination = runs_to_report_destinations.get(run)
        archive_prefix_template = runs_to_archive_prefix_templates.get(run)
        report_tarball_basename = (
            render_template(archive_prefix_template, run) + ".tar"
            if archive_prefix_template
            else default_report_tarball_basename(run)
        )
        replay_target = replay_destination / replay_output_basename(run)
        replay_output_exists = replay_target.exists()
        report_output_exists = False
        tarball_output_exists = False
        if report_destination is not None:
            report_target = report_destination / replay_report_basename(run)
            report_output_exists = report_target.exists()
            tarball_target = report_destination / report_tarball_basename
            tarball_output_exists = tarball_target.exists()
        plans.append(
            RunPlan(
                run=run,
                job_name=job_name,
                variants=tuple(sorted(runs_to_variants[run])),
                resources=resources,
                representative_manifest=runs_to_manifest[run],
                replay_destination=replay_destination,
                report_destination=report_destination,
                report_tarball_basename=report_tarball_basename,
                cache_file=cache_file,
                cache_ready=cache_ready,
                cache_reason=cache_reason,
                replay_output_exists=replay_output_exists,
                report_output_exists=report_output_exists,
                tarball_output_exists=tarball_output_exists,
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
    replay_output_name = replay_output_basename(plan.run)
    replay_report_name = replay_report_basename(plan.run)
    replay_report_tarball_name = plan.report_tarball_basename
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
    ]
    if not plan.replay_output_exists:
        cmd.extend(
            [
                "-output",
                replay_output_name,
                to_mss_output_file_uri(plan.replay_destination / replay_output_name),
            ]
        )
    if plan.report_destination is not None and not plan.report_output_exists:
        cmd.extend(
            [
                "-output",
                replay_report_name,
                to_mss_output_file_uri(plan.report_destination / replay_report_name),
            ]
        )
    if plan.report_destination is not None and replay_report_tarball_name is not None and not plan.tarball_output_exists:
        cmd.extend(
            [
                "-output",
                replay_report_tarball_name,
                to_mss_output_file_uri(plan.report_destination / replay_report_tarball_name),
            ]
        )
    if replay_report_tarball_name is not None:
        cmd.extend(
            [
                "/usr/bin/env",
                f"SWIF_REPORT_TARBALL_BASENAME={replay_report_tarball_name}",
                replay_script,
                str(plan.run),
            ]
        )
    else:
        cmd.extend([replay_script, str(plan.run)])
    return cmd


def to_mss_output_dir_uri(destination: Path) -> str:
    text = str(destination).rstrip("/")
    if text.startswith("mss:"):
        return text if text.endswith("/") else text + "/"
    return f"mss:{text}/"


def to_mss_output_file_uri(destination: Path) -> str:
    text = str(destination)
    if text.startswith("mss:"):
        return text
    return f"mss:{text}"


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
    softenv_wrapper = str(expand_path(DEFAULT_SOFTENV_WRAPPER))
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
        softenv_wrapper,
        args.python_bin,
        str(expand_path(args.jasmine_script)),
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
    partition: str,
    args: argparse.Namespace,
    size_warnings: Sequence[str],
    existing_job_names: Set[str],
) -> None:
    print(f"Family prefix   : {family_prefix}")
    print(f"Manifest directory: {expand_path(args.manifest_dir)}")
    print(f"Workflow         : {workflow}")
    print(f"Replay script    : {expand_path(args.replay_script)}")
    print("Replay MSS copy  : SWIF -output")
    print(f"Account          : {args.account}")
    if args.partition:
        print(f"Partition        : {partition} (CLI override)")
    else:
        print(f"Partition        : {partition} (manifest defaults.partition)")
    print(f"Runs discovered  : {len(runs_to_variants)} unique")
    print(f"Variants matched : {len(variants)}")
    print()

    print("Matched JSON variants")
    print("-" * 80)
    for variant in variants:
        print(f"* {variant.path.name}")
        print(f"    partition: {variant.partition}")
        for runs_file in variant.runs_files:
            print(f"    runs_file: {runs_file}")
        for destination in variant.replay_destinations:
            print(f"    replay_mss_destination: {destination}")
        for destination in variant.report_destinations:
            print(f"    report_mss_destination: {destination}")
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
        if not plan_has_missing_mss_output(plan):
            print(f"[SKIP existing_mss_output] run={plan.run} job={plan.job_name}")
            print(f"         existing_mss_output: {replay_mss_output_file(plan)}")
            if plan.report_destination is not None:
                print(f"         existing_mss_output: {replay_report_mss_output_file(plan)}")
            if plan.report_tarball_basename is not None and replay_report_tarball_mss_output_file(plan) is not None:
                print(f"         existing_mss_output: {replay_report_tarball_mss_output_file(plan)}")
            continue
        status_parts: List[str] = []
        if plan.replay_output_exists:
            status_parts.append("root_exists")
        if plan.report_output_exists:
            status_parts.append("report_exists")
        if plan.tarball_output_exists:
            status_parts.append("tarball_exists")
        status = "SKIP existing" if plan.job_name in existing_job_names and args.skip_existing else "ADD"
        if status_parts:
            status += " (" + ", ".join(status_parts) + ")"
        print(
            f"[{status}] run={plan.run} job={plan.job_name} "
            f"ram={plan.resources.ram} disk={plan.resources.disk} time={plan.resources.time} "
            f"cores={plan.resources.cores} size={format_bytes(plan.resources.estimated_input_bytes)} "
            f"heuristic={plan.resources.reason}"
        )
        print(f"         variants: {', '.join(plan.variants)}")
        print(f"         replay_mss_destination: {plan.replay_destination}")
        if plan.replay_output_exists:
            print(f"         existing_root_mss_output: {replay_mss_output_file(plan)}")
        if plan.report_destination is not None:
            print(f"         report_mss_destination: {plan.report_destination}")
            if plan.report_output_exists:
                print(f"         existing_report_mss_output: {replay_report_mss_output_file(plan)}")
            if plan.report_tarball_basename is not None:
                print(f"         report_tarball: {plan.report_tarball_basename}")
                if plan.tarball_output_exists:
                    print(f"         existing_tarball_mss_output: {replay_report_tarball_mss_output_file(plan)}")
        cmd = build_add_job_command(
            args.swif2_bin,
            workflow,
            args.account,
            partition,
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
                partition,
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
    partition: str,
    args: argparse.Namespace,
) -> int:
    replay_script = str(expand_path(args.replay_script))
    existing_names = set(existing_job_names)
    rc = 0
    eligible_plans = [plan for plan in plans if plan.cache_ready and plan_has_missing_mss_output(plan)]

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
        if not plan_has_missing_mss_output(plan):
            print(f"Skipping run {plan.run}; all requested MSS outputs already exist:")
            print(f"  {replay_mss_output_file(plan)}")
            if plan.report_destination is not None:
                print(f"  {replay_report_mss_output_file(plan)}")
            if replay_report_tarball_mss_output_file(plan) is not None:
                print(f"  {replay_report_tarball_mss_output_file(plan)}")
            print()
            continue
        replay_exists = plan.job_name in existing_names
        if args.skip_existing and replay_exists:
            print(f"Skipping existing job: {plan.job_name}")
        else:
            cmd = build_add_job_command(
                args.swif2_bin,
                workflow,
                args.account,
                partition,
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
                stderr_text = result.stderr or ""
                stdout_text = result.stdout or ""
                if "File already exists on tape" in stderr_text or "File already exists on tape" in stdout_text:
                    print(f"Skipping run {plan.run}; SWIF reports an existing MSS output conflict.")
                    print()
                    continue
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
            partition,
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


def validate_softenv_wrapper(path_text: str) -> None:
    path = expand_path(path_text)
    if not path.exists():
        raise FileNotFoundError(f"Softenv wrapper does not exist: {path}")
    if not os.access(path, os.X_OK):
        raise PermissionError(f"Softenv wrapper is not executable: {path}")



def main() -> int:
    args = parse_args()
    # Jasmine uploads are intentionally run manually from an interactive ifarm
    # session, not from SWIF worker jobs.
    args.auto_jasmine = False
    family_prefix = normalize_family(args.q2, args.w)
    json_dir = expand_path(args.manifest_dir)
    validate_replay_script(args.replay_script)

    variants = discover_json_variants(json_dir, family_prefix, args.family_regex)
    if not variants:
        raise SystemExit(f"No JSON variants found for family prefix {family_prefix} in {json_dir}")
    partition = args.partition or resolve_family_partition(variants)

    (
        runs_to_variants,
        runs_to_manifest,
        runs_to_destination_sets,
        runs_to_report_destination_sets,
        runs_to_archive_prefix_template_sets,
    ) = collect_runs(variants)
    if not runs_to_variants:
        raise SystemExit("No runs were discovered from the matched JSON files")
    runs_to_destinations = resolve_run_destinations(runs_to_destination_sets)
    runs_to_report_destinations = resolve_optional_run_destinations(runs_to_report_destination_sets)
    runs_to_archive_prefix_templates = resolve_optional_run_templates(runs_to_archive_prefix_template_sets)

    plans, size_warnings = build_run_plans(
        family_prefix,
        runs_to_variants,
        runs_to_manifest,
        runs_to_destinations,
        runs_to_report_destinations,
        runs_to_archive_prefix_templates,
        args,
    )
    workflow = workflow_name(args, family_prefix)

    existing_job_names = get_existing_job_names(args.swif2_bin, workflow) if swif_status_exists(args.swif2_bin, workflow) else set()

    print_summary(
        family_prefix,
        variants,
        runs_to_variants,
        plans,
        workflow,
        partition,
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

    return submit_jobs(workflow, family_prefix, plans, existing_job_names, partition, args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
