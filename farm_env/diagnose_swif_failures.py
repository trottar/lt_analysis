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
from typing import Dict, Iterable, List, Sequence

DEFAULT_SWIF2 = os.environ.get("SWIF2_BIN", "swif2")

MEMORY_PAT = re.compile(r"(out[_ -]?of[_ -]?memory|oom|memory limit|maxrss|exceeded memory)", re.I)
TIME_PAT = re.compile(r"(time[_ -]?limit|timed out|due to time limit|wall ?time)", re.I)
DISK_PAT = re.compile(r"(no space left on device|disk full|disk quota|scratch|no space)", re.I)
SITE_PAT = re.compile(r"(site[_ -]?launch[_ -]?fail|invalid account|partition)", re.I)
CACHE_PAT = re.compile(r"(jcache|cache|coin_all_|raw file|missing raw|staging)", re.I)
PATH_PAT = re.compile(r"(not found|no such file|does not exist|permission denied|macro .* not found)", re.I)


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


def parse_show_job_output(text: str) -> Dict[str, str]:
    data: Dict[str, str] = {}
    for raw in text.splitlines():
        if "=" not in raw:
            continue
        key, value = raw.split("=", 1)
        data[key.strip()] = value.strip()
    return data


def classify_failure(problem: str, details: str) -> tuple[str, str]:
    haystack = f"{problem} {details}"
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


def grouped(items: Iterable[FailedJob]) -> Dict[str, List[FailedJob]]:
    groups: Dict[str, List[FailedJob]] = defaultdict(list)
    for item in items:
        groups[item.category].append(item)
    return dict(groups)


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
    else:
        print(f"  Next step   : inspect logs for one job, then retry if appropriate")
    print(f"  Inspect one : {swif2_bin} show-job -workflow {workflow} -name {jobs[0].job_name}")
    print(f"  Retry group : {swif2_bin} retry-jobs {workflow} -names {names}")
    print()


def main() -> int:
    args = parse_args()
    problem_names = load_problem_job_names(args.swif2_bin, args.workflow)
    if not problem_names:
        print("No current failed/problem jobs found.")
        return 0

    if args.limit > 0:
        problem_names = problem_names[: args.limit]

    failed_jobs = [inspect_failed_job(args.swif2_bin, args.workflow, job_name) for job_name in problem_names]
    groups = grouped(failed_jobs)

    print(f"Workflow: {args.workflow}")
    print(f"Failed jobs inspected: {len(failed_jobs)}")
    print()

    print("Failure groups")
    print("-" * 80)
    for group_name in sorted(groups):
        print_group(group_name, groups[group_name], args.workflow, args.swif2_bin)

    resource_groups = [name for name in groups if name.startswith("resource_")]
    non_resource_groups = [name for name in groups if not name.startswith("resource_")]

    print("Recommended flow")
    print("-" * 80)
    if non_resource_groups:
        print("* For script/path/cache-like failures, fix the underlying issue first, then retry only those failed jobs.")
    if resource_groups:
        print("* For memory/time/disk failures, rebalance resources before retrying.")
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
