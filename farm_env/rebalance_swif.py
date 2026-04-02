#!/usr/bin/env python3
"""
Inspect a SWIF2 workflow, identify likely resource-related problem jobs,
and prepare or apply resource bumps using swif2 modify-jobs.

Dry-run is the default. Use --apply to execute the modify-jobs commands and
optionally 'swif2 run' afterward.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Sequence, Set

DEFAULT_SWIF2 = os.environ.get("SWIF2_BIN", "swif2")

MEMORY_PAT = re.compile(r"(out[_ -]?of[_ -]?memory|oom|memory limit|maxrss|exceeded memory)", re.I)
TIME_PAT = re.compile(r"(time[_ -]?limit|timed out|due to time limit|wall ?time)", re.I)
DISK_PAT = re.compile(r"(no space left on device|disk|scratch)", re.I)


@dataclass(frozen=True)
class ProblemJob:
    job_name: str
    problem: str
    details: str
    categories: frozenset[str]
    ram_bytes: int = 0
    disk_bytes: int = 0
    time_secs: int = 0



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Bump SWIF2 job resources for KaonLT workflows.")
    parser.add_argument("workflow", help="SWIF2 workflow name")
    parser.add_argument("--swif2-bin", default=DEFAULT_SWIF2, help="SWIF2 client executable")
    parser.add_argument("--ram-mult", type=float, default=1.5, help="RAM multiplier for OOM-like failures")
    parser.add_argument("--disk-mult", type=float, default=1.5, help="Disk multiplier for disk-like failures")
    parser.add_argument("--time-mult", type=float, default=1.5, help="Time multiplier for timeout-like failures")
    parser.add_argument("--apply", action="store_true", help="Execute modify-jobs and rerun the workflow")
    parser.add_argument("--no-run", action="store_true", help="With --apply, do not call 'swif2 run' afterwards")
    return parser.parse_args()



def run_command(cmd: Sequence[str], capture: bool = False) -> subprocess.CompletedProcess:
    return subprocess.run(list(cmd), text=True, capture_output=capture, check=False)



def load_problem_jobs(swif2_bin: str, workflow: str) -> List[ProblemJob]:
    cmd = [swif2_bin, "status", workflow, "-jobs", "-display", "json"]
    result = run_command(cmd, capture=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to query workflow {workflow}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

    payload = json.loads(result.stdout)
    jobs = payload.get("jobs", [])
    problems: List[ProblemJob] = []

    for job in jobs:
        if not isinstance(job, dict):
            continue
        if str(job.get("job_attempt_status", "")).lower() != "problem":
            continue
        job_name = str(job.get("job_name", ""))
        if not job_name:
            continue
        problem = str(job.get("job_attempt_problem", ""))
        details = str(job.get("job_attempt_problem_details", ""))
        ram_bytes = int(job.get("site_job_ram_bytes") or 0)
        disk_bytes = int(job.get("site_job_disk_bytes") or 0)
        time_secs = int(job.get("site_job_time_secs") or 0)
        haystack = f"{problem} {details}"
        categories: Set[str] = set()
        if MEMORY_PAT.search(haystack):
            categories.add("memory")
        if TIME_PAT.search(haystack):
            categories.add("time")
        if DISK_PAT.search(haystack):
            categories.add("disk")
        if categories:
            problems.append(
                ProblemJob(
                    job_name=job_name,
                    problem=problem,
                    details=details,
                    categories=frozenset(categories),
                    ram_bytes=ram_bytes,
                    disk_bytes=disk_bytes,
                    time_secs=time_secs,
                )
            )
    return problems



def build_modify_commands(args: argparse.Namespace, problems: Sequence[ProblemJob]) -> List[List[str]]:
    grouped: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    for problem in problems:
        if "memory" in problem.categories and problem.ram_bytes > 0:
            new_ram = str(int(math.ceil(problem.ram_bytes * args.ram_mult)))
            grouped["memory"][new_ram].append(problem.job_name)
        if "time" in problem.categories and problem.time_secs > 0:
            new_time = str(int(math.ceil(problem.time_secs * args.time_mult)))
            grouped["time"][new_time].append(problem.job_name)
        if "disk" in problem.categories and problem.disk_bytes > 0:
            new_disk = str(int(math.ceil(problem.disk_bytes * args.disk_mult)))
            grouped["disk"][new_disk].append(problem.job_name)

    commands: List[List[str]] = []
    for new_ram, names in sorted(grouped.get("memory", {}).items()):
        commands.append(
            [
                args.swif2_bin,
                "modify-jobs",
                args.workflow,
                "-names",
                "-ram",
                "set",
                new_ram,
                *sorted(set(names)),
            ]
        )
    for new_time, names in sorted(grouped.get("time", {}).items()):
        commands.append(
            [
                args.swif2_bin,
                "modify-jobs",
                args.workflow,
                "-names",
                "-time",
                "set",
                new_time,
                *sorted(set(names)),
            ]
        )
    for new_disk, names in sorted(grouped.get("disk", {}).items()):
        commands.append(
            [
                args.swif2_bin,
                "modify-jobs",
                args.workflow,
                "-names",
                "-disk",
                "set",
                new_disk,
                *sorted(set(names)),
            ]
        )
    return commands



def main() -> int:
    args = parse_args()
    problems = load_problem_jobs(args.swif2_bin, args.workflow)

    if not problems:
        print("No problem jobs with obvious resource-related signatures were found.")
        return 0

    print("Detected resource-related problem jobs")
    print("-" * 80)
    for problem in problems:
        print(f"* {problem.job_name}: categories={','.join(sorted(problem.categories))}")
        print(f"    problem: {problem.problem}")
        if problem.details:
            print(f"    details: {problem.details}")
    print()

    commands = build_modify_commands(args, problems)
    if not commands:
        print("No modify-jobs commands were generated.")
        return 0

    print("Planned commands")
    print("-" * 80)
    for cmd in commands:
        print(" ".join(shlex.quote(x) for x in cmd))
    if args.apply and not args.no_run:
        print(" ".join(shlex.quote(x) for x in [args.swif2_bin, "run", args.workflow]))
    print()

    if not args.apply:
        print("Dry run only. Use --apply to execute these modifications.")
        return 0

    for cmd in commands:
        result = run_command(cmd, capture=True)
        if result.stdout:
            print(result.stdout.rstrip())
        if result.stderr:
            print(result.stderr.rstrip(), file=sys.stderr)
        if result.returncode != 0:
            return result.returncode

    if not args.no_run:
        result = run_command([args.swif2_bin, "run", args.workflow], capture=True)
        if result.stdout:
            print(result.stdout.rstrip())
        if result.stderr:
            print(result.stderr.rstrip(), file=sys.stderr)
        return result.returncode

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
