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
SUBMIT_TIMEOUT_PAT = re.compile(r"Command timed out:\s*sbatch\b", re.I)


@dataclass(frozen=True)
class ProblemJob:
    job_name: str
    problem: str
    details: str
    categories: frozenset[str]
    num_attempts: int = 0
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


def compute_targets(args: argparse.Namespace, problem: ProblemJob) -> Dict[str, int]:
    targets: Dict[str, int] = {}
    if "memory" in problem.categories and problem.ram_bytes > 0:
        targets["ram"] = int(math.ceil(problem.ram_bytes * args.ram_mult))
    if "time" in problem.categories and problem.time_secs > 0:
        targets["time"] = int(math.ceil(problem.time_secs * args.time_mult))
    if "disk" in problem.categories and problem.disk_bytes > 0:
        targets["disk"] = int(math.ceil(problem.disk_bytes * args.disk_mult))
    return targets



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
        num_attempts = int(job.get("num_attempts") or 0)
        ram_bytes = int(job.get("site_job_ram_bytes") or 0)
        disk_bytes = int(job.get("site_job_disk_bytes") or 0)
        time_secs = int(job.get("site_job_time_secs") or 0)
        if SUBMIT_TIMEOUT_PAT.search(details) or (
            problem.strip().upper() == "SWIF_SYSTEM_ERROR"
            and "sbatch" in details
            and "timed out" in details.lower()
        ):
            # This is a scheduler submit timeout before the job actually starts.
            # Retrying is the right default; resource bumping is not.
            continue
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
                    num_attempts=num_attempts,
                    ram_bytes=ram_bytes,
                    disk_bytes=disk_bytes,
                    time_secs=time_secs,
                )
            )
    return problems



def build_modify_commands(args: argparse.Namespace, problems: Sequence[ProblemJob]) -> List[List[str]]:
    grouped: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    for problem in problems:
        targets = compute_targets(args, problem)
        if "ram" in targets:
            new_ram = str(targets["ram"])
            grouped["memory"][new_ram].append(problem.job_name)
        if "time" in targets:
            new_time = str(targets["time"])
            grouped["time"][new_time].append(problem.job_name)
        if "disk" in targets:
            new_disk = str(targets["disk"])
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
        targets = compute_targets(args, problem)
        print(f"* {problem.job_name}: categories={','.join(sorted(problem.categories))}")
        print(f"    problem: {problem.problem}")
        if problem.num_attempts:
            print(f"    attempts: {problem.num_attempts}")
        if "ram" in targets:
            print(f"    ram: {problem.ram_bytes} -> {targets['ram']}")
        if "time" in targets:
            print(f"    time: {problem.time_secs} -> {targets['time']}")
        if "disk" in targets:
            print(f"    disk: {problem.disk_bytes} -> {targets['disk']}")
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
