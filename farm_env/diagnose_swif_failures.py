#!/usr/bin/env python3
"""
Diagnose failed SWIF2 jobs for a workflow.

This helper:
- filters the workflow to only current problem jobs
- calls `swif2 show-job` for each failed job
- groups failures into resource-related vs other categories
- prints suggested next actions for each group
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

DEFAULT_SWIF2 = os.environ.get("SWIF2_BIN", "swif2")

MEMORY_PAT = re.compile(r"(out[_ -]?of[_ -]?memory|oom|memory limit|maxrss|exceeded memory)", re.I)
TIME_PAT = re.compile(r"(time[_ -]?limit|timed out|due to time limit|wall ?time)", re.I)
DISK_PAT = re.compile(r"(no space left on device|disk full|disk quota|scratch|no space)", re.I)
SUBMIT_TIMEOUT_PAT = re.compile(r"Command timed out:\s*sbatch\b", re.I)
OUTPUT_EXISTS_PAT = re.compile(r"File already exists on tape:", re.I)
SITE_PAT = re.compile(r"(site[_ -]?launch[_ -]?fail|invalid account|partition)", re.I)
CACHE_PAT = re.compile(r"(jcache|cache|coin_all_|raw file|missing raw|staging)", re.I)
PATH_PAT = re.compile(r"(not found|no such file|does not exist|permission denied|macro .* not found)", re.I)

CACHE_ROOT_CANDIDATES = tuple(
    dict.fromkeys(
        [
            os.environ.get("CACHEPATH", "").strip(),
            "/lustre/expphy/cache/hallc/kaonlt",
            "/cache/hallc/kaonlt",
        ]
    )
)


@dataclass(frozen=True)
class FailedJob:
    job_name: str
    problem: str
    details: str
    category: str
    suggestion: str
    slurm_exitcode: str = ""
    slurm_state: str = ""
    slurm_id: str = ""


@dataclass(frozen=True)
class RootOutputCheck:
    job_name: str
    local_name: str
    remote_path: str
    cache_path: str
    status: str
    detail: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Diagnose failed SWIF2 jobs for a workflow.")
    parser.add_argument("workflow", help="SWIF2 workflow name")
    parser.add_argument("--swif2-bin", default=DEFAULT_SWIF2, help="SWIF2 client executable")
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Limit how many failed jobs to inspect with show-job (default: inspect all)",
    )
    return parser.parse_args()


def run_command(cmd: Sequence[str], capture: bool = False) -> subprocess.CompletedProcess:
    return subprocess.run(list(cmd), text=True, capture_output=capture, check=False)


def load_problem_job_names(swif2_bin: str, workflow: str) -> List[str]:
    cmd = [swif2_bin, "status", workflow, "-jobs", "-display", "json"]
    result = run_command(cmd, capture=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to query workflow {workflow}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    payload = json.loads(result.stdout)
    jobs = payload.get("jobs", [])
    names: List[str] = []
    for job in jobs:
        if not isinstance(job, dict):
            continue
        if str(job.get("job_attempt_status", "")).lower() != "problem":
            continue
        job_name = str(job.get("job_name", "")).strip()
        if job_name:
            names.append(job_name)
    return names


def load_done_job_names(swif2_bin: str, workflow: str) -> List[str]:
    cmd = [swif2_bin, "status", workflow, "-jobs", "-display", "json"]
    result = run_command(cmd, capture=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to query workflow {workflow}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    payload = json.loads(result.stdout)
    jobs = payload.get("jobs", [])
    names: List[str] = []
    for job in jobs:
        if not isinstance(job, dict):
            continue
        if str(job.get("job_attempt_status", "")).lower() != "done":
            continue
        job_name = str(job.get("job_name", "")).strip()
        if job_name:
            names.append(job_name)
    return names


def parse_show_job_output(text: str) -> Dict[str, str]:
    data: Dict[str, str] = {}
    for raw in text.splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        data[key.strip()] = value.strip()
    return data


def parse_show_job_outputs(text: str) -> List[Tuple[str, str]]:
    outputs: List[Tuple[str, str]] = []
    pending_local: Optional[str] = None
    for raw in text.splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key == "local":
            pending_local = value
        elif key == "remote":
            outputs.append((pending_local or "", value))
            pending_local = None
    return outputs


def classify_failure(problem: str, details: str) -> tuple[str, str]:
    haystack = f"{problem} {details}"
    if problem.strip().upper() == "SWIF_OUTPUT_FAIL" and OUTPUT_EXISTS_PAT.search(details):
        return (
            "output_exists_on_tape",
            "These jobs completed, but SWIF could not export because the target MSS file already exists; do not retry unchanged jobs.",
        )
    if SUBMIT_TIMEOUT_PAT.search(details) or (
        problem.strip().upper() == "SWIF_SYSTEM_ERROR" and "sbatch" in details and "timed out" in details.lower()
    ):
        return (
            "scheduler_submit_timeout",
            "This is a SWIF/Slurm submit timeout before the job starts; retry these jobs first rather than increasing wall time.",
        )
    if MEMORY_PAT.search(haystack):
        return ("resource_memory", "Increase RAM for these jobs, then rerun them.")
    if TIME_PAT.search(haystack):
        return ("resource_time", "Increase wall time for these jobs, then rerun them.")
    if DISK_PAT.search(haystack):
        return ("resource_disk", "Increase disk for these jobs, then rerun them.")
    if SITE_PAT.search(haystack):
        return ("site_or_scheduler", "Check account/partition/site config, then retry the failed jobs.")
    if CACHE_PAT.search(haystack):
        return ("cache_or_staging", "Check cache/tape staging first, then retry once inputs are ready.")
    if PATH_PAT.search(haystack):
        return ("path_or_missing_file", "Fix the missing-path/file issue in the script or environment, then retry.")
    if problem.strip().upper() == "SLURM_FAILED" and "Exited with code 1" in details:
        return ("script_exit_1", "Inspect the job .out/.err logs for the real script error before retrying.")
    return ("other", "Inspect the per-job .out/.err logs, then decide whether to retry or change resources.")


def inspect_failed_job(swif2_bin: str, workflow: str, job_name: str) -> FailedJob:
    cmd = [swif2_bin, "show-job", "-workflow", workflow, "-name", job_name]
    result = run_command(cmd, capture=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to inspect job {job_name}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    parsed = parse_show_job_output(result.stdout)
    problem = parsed.get("job_attempt_problem", "")
    details = parsed.get("job_attempt_problem_details", "")
    category, suggestion = classify_failure(problem, details)
    return FailedJob(
        job_name=job_name,
        problem=problem,
        details=details,
        category=category,
        suggestion=suggestion,
        slurm_exitcode=parsed.get("slurm_exitcode", ""),
        slurm_state=parsed.get("slurm_state", ""),
        slurm_id=parsed.get("slurm_id", ""),
    )


def validate_root_file(path: Path) -> Tuple[str, str]:
    if not path.exists():
        return ("missing_cache_file", "file does not exist in cache")
    try:
        size_bytes = path.stat().st_size
    except OSError as exc:
        return ("cache_stat_error", str(exc))
    if size_bytes <= 0:
        return ("zero_size_file", "file exists but has zero size")

    root_exc: Optional[Exception] = None
    try:
        import ROOT  # type: ignore

        ROOT.gROOT.SetBatch(True)
        previous_error_level = int(getattr(ROOT, "gErrorIgnoreLevel", 0))
        try:
            ROOT.gErrorIgnoreLevel = ROOT.kError
            root_file = ROOT.TFile.Open(str(path), "READ")
        finally:
            ROOT.gErrorIgnoreLevel = previous_error_level
        if not root_file:
            return ("unreadable_root", "ROOT.TFile.Open returned null")
        try:
            if not root_file.IsOpen():
                return ("unreadable_root", "ROOT opened a handle but file is not open")
            if root_file.IsZombie():
                return ("zombie_root_file", "ROOT marked file as zombie")
            if root_file.TestBit(ROOT.TFile.kRecovered):
                return ("recovered_root_file", "ROOT marked file as recovered")
            keys = root_file.GetListOfKeys()
            key_count = int(keys.GetSize()) if keys is not None else 0
            if key_count <= 0:
                return ("root_file_has_no_keys", f"ROOT opened file but found no keys; size={size_bytes}")
            return ("healthy_root_file", f"ROOT opened successfully; keys={key_count}; size={size_bytes}")
        finally:
            root_file.Close()
    except Exception as exc:  # pragma: no cover - depends on runtime env
        root_exc = exc

    try:
        import uproot  # type: ignore

        with uproot.open(path) as handle:
            key_count = len(handle.keys())
        if key_count <= 0:
            return ("root_file_has_no_keys", f"uproot opened file but found no keys; size={size_bytes}")
        return ("healthy_root_file", f"uproot opened successfully; keys={key_count}; size={size_bytes}")
    except Exception as exc:  # pragma: no cover - depends on runtime env
        if root_exc is not None:
            return ("unreadable_root", f"ROOT failed: {root_exc}; uproot failed: {exc}")
        return ("unreadable_root", f"uproot failed: {exc}")


def derive_cache_path_from_mss(mss_path: Path) -> Path:
    mss_text = str(mss_path)
    candidates: List[Path] = []
    if mss_text.startswith("/mss/hallc/kaonlt"):
        suffix = mss_text[len("/mss/hallc/kaonlt"):].lstrip("/")
        for root in CACHE_ROOT_CANDIDATES:
            if root:
                candidates.append(Path(root).expanduser() / suffix)
    elif mss_text.startswith("/mss/"):
        candidates.append(Path("/cache" + mss_text[4:]).expanduser())
        candidates.append(Path("/lustre/expphy/cache" + mss_text[4:]).expanduser())
    else:
        candidates.append(mss_path)

    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0] if candidates else mss_path


def inspect_success_root_outputs(swif2_bin: str, workflow: str, job_name: str) -> List[RootOutputCheck]:
    cmd = [swif2_bin, "show-job", "-workflow", workflow, "-name", job_name]
    result = run_command(cmd, capture=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Failed to inspect job {job_name}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )

    outputs = parse_show_job_outputs(result.stdout)
    if not outputs:
        return []

    checks: List[RootOutputCheck] = []
    for local_name, remote_name in outputs:
        if not remote_name.lower().endswith(".root"):
            continue
        remote_path = Path(remote_name)
        cache_path = derive_cache_path_from_mss(remote_path)
        status, detail = validate_root_file(cache_path)
        checks.append(
            RootOutputCheck(
                job_name=job_name,
                local_name=local_name,
                remote_path=str(remote_path),
                cache_path=str(cache_path),
                status=status,
                detail=detail,
            )
        )
    return checks


def grouped(items: Iterable[FailedJob]) -> Dict[str, List[FailedJob]]:
    groups: Dict[str, List[FailedJob]] = defaultdict(list)
    for item in items:
        groups[item.category].append(item)
    return dict(groups)


def format_root_check_progress(job_name: str, current: int, total: int) -> str:
    match = re.search(r"run(\d+)", job_name, re.I)
    run_label = match.group(1) if match else "unknown"
    percent = (100.0 * current / total) if total else 100.0
    return f"[ROOT check {current}/{total} | {percent:5.1f}%] run={run_label} job={job_name}"


def emit_root_check_progress(message: str, inline: bool, previous_width: int) -> int:
    if not inline:
        print(message, flush=True)
        return previous_width

    padded = message.ljust(previous_width)
    sys.stdout.write("\r" + padded)
    sys.stdout.flush()
    return max(previous_width, len(message))


def print_group(name: str, jobs: Sequence[FailedJob], workflow: str, swif2_bin: str) -> None:
    print(f"[{name}] count={len(jobs)}")
    for job in jobs:
        print(f"* {job.job_name}")
        print(f"    problem    : {job.problem}")
        if job.details:
            print(f"    details    : {job.details}")
        if job.slurm_state or job.slurm_exitcode:
            print(f"    slurm      : state={job.slurm_state or '?'} exit={job.slurm_exitcode or '?'} id={job.slurm_id or '?'}")
    print(f"  Suggestion  : {jobs[0].suggestion}")
    names = " ".join(job.job_name for job in jobs)
    if name.startswith("resource_"):
        print(f"  Next step   : consider resource bumping before retry")
    elif name == "scheduler_submit_timeout":
        print(f"  Next step   : simple retry first; only escalate if this keeps happening")
    else:
        print(f"  Next step   : inspect logs for one job, then retry if appropriate")
    print(f"  Inspect one : {swif2_bin} show-job -workflow {workflow} -name {jobs[0].job_name}")
    print(f"  Retry group : {swif2_bin} retry-jobs {workflow} -names {names}")
    print()


def print_success_root_checks(checks: Sequence[RootOutputCheck]) -> None:
    print("Completed ROOT Cache Checks")
    print("-" * 80)
    if not checks:
        print("No completed ROOT outputs were found to inspect.")
        return

    grouped_checks: Dict[str, List[RootOutputCheck]] = defaultdict(list)
    for check in checks:
        grouped_checks[check.status].append(check)

    total = len(checks)
    healthy = len(grouped_checks.get("healthy_root_file", []))
    suspicious = total - healthy
    print(f"ROOT outputs inspected : {total}")
    print(f"Healthy cache files    : {healthy}")
    print(f"Suspicious cache files : {suspicious}")
    print()

    if suspicious == 0:
        print("All inspected cache ROOT files opened cleanly.")
        return

    for status in sorted(grouped_checks):
        if status == "healthy_root_file":
            continue
        bucket = grouped_checks[status]
        print(f"[{status}] count={len(bucket)}")
        for check in bucket:
            print(f"* {check.job_name}")
            print(f"    local      : {check.local_name}")
            print(f"    remote     : {check.remote_path}")
            print(f"    cache      : {check.cache_path}")
            print(f"    detail     : {check.detail}")
        print()


def main() -> int:
    args = parse_args()
    problem_names = load_problem_job_names(args.swif2_bin, args.workflow)
    done_names = load_done_job_names(args.swif2_bin, args.workflow)

    if args.limit > 0:
        problem_names = problem_names[: args.limit]
        done_names = done_names[: args.limit]

    failed_jobs = [inspect_failed_job(args.swif2_bin, args.workflow, job_name) for job_name in problem_names]
    groups = grouped(failed_jobs) if failed_jobs else {}
    root_checks: List[RootOutputCheck] = []
    completed_total = len(done_names)
    inline_progress = sys.stdout.isatty()
    progress_width = 0
    if completed_total:
        print(f"Checking completed ROOT outputs in cache ({completed_total} jobs)...", flush=True)
    for index, job_name in enumerate(done_names, start=1):
        progress_width = emit_root_check_progress(
            format_root_check_progress(job_name, index, completed_total),
            inline_progress,
            progress_width,
        )
        root_checks.extend(inspect_success_root_outputs(args.swif2_bin, args.workflow, job_name))
    if completed_total and inline_progress:
        print()

    print(f"Workflow: {args.workflow}")
    print(f"Failed jobs inspected: {len(failed_jobs)}")
    print(f"Completed jobs inspected: {len(done_names)}")
    print()

    print("Failure groups")
    print("-" * 80)
    if not groups:
        print("No current failed/problem jobs found.")
        print()
    else:
        for group_name in sorted(groups):
            print_group(group_name, groups[group_name], args.workflow, args.swif2_bin)

    print_success_root_checks(root_checks)
    print()

    resource_groups = [name for name in groups if name.startswith("resource_")]
    non_resource_groups = [name for name in groups if not name.startswith("resource_")]

    print("Recommended flow")
    print("-" * 80)
    if non_resource_groups:
        print("* For script/path/cache-like failures, fix the underlying issue first, then retry only those failed jobs.")
    if resource_groups:
        print("* For memory/time/disk failures, rebalance resources before retrying.")
    if "scheduler_submit_timeout" in groups:
        print("* For scheduler submit timeouts, retry the failed jobs first; these are not wall-time failures inside your script.")
    if "output_exists_on_tape" in groups:
        print("* For output-exists-on-tape jobs, treat the MSS copy as already present; skip/recreate those jobs only if you need a different destination or a clean workflow state.")
    suspicious_root_checks = [check for check in root_checks if check.status != "healthy_root_file"]
    if suspicious_root_checks:
        print("* For suspicious cache ROOT files, inspect or restage those cache entries before trusting the completed job.")
    print("* After retries or resource changes, run:")
    print(f"    {args.swif2_bin} run {args.workflow}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
