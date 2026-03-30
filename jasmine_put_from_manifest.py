#!/usr/bin/env python3
"""
Create tar archives for small files and submit write requests to Jefferson Lab
Jasmine tape storage via jput.

Behavior:
- Reads a JSON manifest describing source paths and /mss destinations.
- Files >= min_size_mb are submitted directly to jput.
- Files < min_size_mb are packed into tar archives.
- Tar archives are packed toward target_archive_size_gb, capped by
  max_archive_size_gb.
- Dry-run by default. Use --submit to actually invoke jput.

The script assumes the archival jput usage forms documented by JLab:
    jput file /mss/dir
    jput file1 file2 ... /mss/dir
    jput -m file1 /mss/dir1 file2 /mss/dir2 ...

If your local jput wrapper differs, adjust build_jput_command().
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tarfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

MB = 1024 ** 2
GB = 1024 ** 3
DEFAULT_MIN_SIZE_MB = 200
DEFAULT_TARGET_ARCHIVE_SIZE_GB = 4
DEFAULT_MAX_ARCHIVE_SIZE_GB = 20
DEFAULT_STAGE_DIR = "/tmp/jasmine_stage"
DEFAULT_JPUT = shutil.which("jput") or "/site/bin/jput"


@dataclass
class SourceFile:
    path: Path
    size_bytes: int
    arcname: str


@dataclass
class PlannedPut:
    local_path: Path
    mss_dir: Path
    reason: str
    size_bytes: int
    members: List[SourceFile] = field(default_factory=list)


@dataclass
class JobConfig:
    source: Path
    destination: Path
    recursive: bool = True
    archive_prefix: Optional[str] = None
    preserve_tree: bool = True


@dataclass
class Settings:
    min_size_bytes: int
    target_archive_size_bytes: int
    max_archive_size_bytes: int
    stage_dir: Path
    jput_bin: str
    submit: bool
    jput_n_flag: bool
    max_pairs_per_request: int
    cleanup_archives: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read a JSON manifest, package sub-threshold files into tar archives, "
            "and submit or print jput requests for /mss tape storage."
        )
    )
    parser.add_argument("manifest", help="Path to JSON manifest")
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Actually call jput. Default is dry-run.",
    )
    parser.add_argument(
        "--keep-archives",
        action="store_true",
        help="Keep staged tar archives after submission.",
    )
    parser.add_argument(
        "--stage-dir",
        default=None,
        help=f"Override staging directory (default: {DEFAULT_STAGE_DIR})",
    )
    return parser.parse_args()


def load_manifest(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def bool_from_manifest(value, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "y", "on"}
    return bool(value)


def resolve_settings(manifest: Dict, cli: argparse.Namespace) -> Settings:
    defaults = manifest.get("defaults", {})
    min_size_mb = float(defaults.get("min_size_mb", DEFAULT_MIN_SIZE_MB))
    target_archive_size_gb = float(
        defaults.get("target_archive_size_gb", DEFAULT_TARGET_ARCHIVE_SIZE_GB)
    )
    max_archive_size_gb = float(
        defaults.get("max_archive_size_gb", DEFAULT_MAX_ARCHIVE_SIZE_GB)
    )

    if target_archive_size_gb <= 0:
        raise ValueError("target_archive_size_gb must be > 0")
    if max_archive_size_gb <= 0:
        raise ValueError("max_archive_size_gb must be > 0")
    if target_archive_size_gb > max_archive_size_gb:
        raise ValueError("target_archive_size_gb cannot exceed max_archive_size_gb")

    stage_dir = Path(cli.stage_dir or defaults.get("stage_dir", DEFAULT_STAGE_DIR)).expanduser()
    jput_bin = str(defaults.get("jput_bin", DEFAULT_JPUT))
    jput_n_flag = bool_from_manifest(defaults.get("jput_n_flag"), False)
    max_pairs_per_request = int(defaults.get("max_pairs_per_request", 32))
    cleanup_archives = not cli.keep_archives

    return Settings(
        min_size_bytes=int(min_size_mb * MB),
        target_archive_size_bytes=int(target_archive_size_gb * GB),
        max_archive_size_bytes=int(max_archive_size_gb * GB),
        stage_dir=stage_dir,
        jput_bin=jput_bin,
        submit=cli.submit,
        jput_n_flag=jput_n_flag,
        max_pairs_per_request=max_pairs_per_request,
        cleanup_archives=cleanup_archives,
    )


def resolve_jobs(manifest: Dict) -> List[JobConfig]:
    jobs = manifest.get("jobs")
    if not isinstance(jobs, list) or not jobs:
        raise ValueError("Manifest must contain a non-empty 'jobs' list")

    resolved: List[JobConfig] = []
    for index, entry in enumerate(jobs, start=1):
        if not isinstance(entry, dict):
            raise ValueError(f"jobs[{index}] must be an object")
        try:
            source = Path(entry["source"]).expanduser().resolve()
            destination = Path(entry["destination"]).expanduser()
        except KeyError as exc:
            raise ValueError(f"jobs[{index}] missing key: {exc.args[0]}") from exc

        recursive = bool_from_manifest(entry.get("recursive"), True)
        preserve_tree = bool_from_manifest(entry.get("preserve_tree"), True)
        archive_prefix = entry.get("archive_prefix")
        resolved.append(
            JobConfig(
                source=source,
                destination=destination,
                recursive=recursive,
                archive_prefix=archive_prefix,
                preserve_tree=preserve_tree,
            )
        )
    return resolved


def ensure_supported_destination(destination: Path) -> None:
    if not str(destination).startswith("/mss/"):
        raise ValueError(f"Destination must be under /mss/: {destination}")



def iter_source_files(job: JobConfig) -> Iterator[SourceFile]:
    source = job.source
    if not source.exists():
        raise FileNotFoundError(f"Source path does not exist: {source}")

    if source.is_file():
        yield SourceFile(
            path=source,
            size_bytes=source.stat().st_size,
            arcname=source.name,
        )
        return

    if not source.is_dir():
        raise ValueError(f"Unsupported source type: {source}")

    if not job.recursive:
        children = sorted([p for p in source.iterdir() if p.is_file()])
    else:
        children = sorted([p for p in source.rglob("*") if p.is_file()])

    for child in children:
        arcname = child.relative_to(source).as_posix() if job.preserve_tree else child.name
        yield SourceFile(
            path=child,
            size_bytes=child.stat().st_size,
            arcname=arcname,
        )



def chunk_small_files(
    files: Sequence[SourceFile],
    target_bytes: int,
    max_bytes: int,
) -> List[List[SourceFile]]:
    """
    Greedy first-fit-decreasing binning toward the target size, bounded by max size.
    This is simple and stable for tape-prep purposes.
    """
    chunks: List[List[SourceFile]] = []
    chunk_sizes: List[int] = []

    for item in sorted(files, key=lambda f: f.size_bytes, reverse=True):
        placed = False

        # Prefer a chunk still below target.
        preferred_indices = [
            i for i, size in enumerate(chunk_sizes)
            if size < target_bytes and size + item.size_bytes <= max_bytes
        ]
        preferred_indices.sort(key=lambda i: chunk_sizes[i], reverse=True)

        for idx in preferred_indices:
            chunks[idx].append(item)
            chunk_sizes[idx] += item.size_bytes
            placed = True
            break

        if placed:
            continue

        # Otherwise, fit anywhere under the hard max.
        fallback_indices = [
            i for i, size in enumerate(chunk_sizes)
            if size + item.size_bytes <= max_bytes
        ]
        fallback_indices.sort(key=lambda i: chunk_sizes[i], reverse=True)
        for idx in fallback_indices:
            chunks[idx].append(item)
            chunk_sizes[idx] += item.size_bytes
            placed = True
            break

        if not placed:
            chunks.append([item])
            chunk_sizes.append(item.size_bytes)

    return chunks



def plan_job(job: JobConfig, settings: Settings) -> List[PlannedPut]:
    ensure_supported_destination(job.destination)
    files = list(iter_source_files(job))
    if not files:
        raise ValueError(f"No files found under source: {job.source}")

    direct_files = [f for f in files if f.size_bytes >= settings.min_size_bytes]
    small_files = [f for f in files if f.size_bytes < settings.min_size_bytes]

    planned: List[PlannedPut] = []
    for entry in direct_files:
        planned.append(
            PlannedPut(
                local_path=entry.path,
                mss_dir=job.destination,
                reason="direct>=threshold",
                size_bytes=entry.size_bytes,
                members=[entry],
            )
        )

    if small_files:
        groups = chunk_small_files(
            small_files,
            target_bytes=settings.target_archive_size_bytes,
            max_bytes=settings.max_archive_size_bytes,
        )
        base_prefix = job.archive_prefix or sanitize_name(job.source.name or "archive")
        for idx, group in enumerate(groups, start=1):
            planned.append(
                PlannedPut(
                    local_path=settings.stage_dir / f"{base_prefix}_{idx:04d}.tar",
                    mss_dir=job.destination,
                    reason="tarball<200MB_inputs",
                    size_bytes=sum(member.size_bytes for member in group),
                    members=group,
                )
            )

    return planned



def sanitize_name(name: str) -> str:
    safe = [c if c.isalnum() or c in {"-", "_", "."} else "_" for c in name]
    text = "".join(safe).strip("._")
    return text or "archive"



def create_tarball(plan: PlannedPut) -> None:
    plan.local_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(plan.local_path, mode="w") as tar:
        for member in plan.members:
            tar.add(member.path, arcname=member.arcname, recursive=False)



def maybe_check_existing_stub(plan: PlannedPut) -> Optional[Path]:
    # Best-effort local check only. Jasmine itself remains authoritative.
    candidate = plan.mss_dir / plan.local_path.name
    return candidate if candidate.exists() else None



def build_jput_command(pairs: Sequence[Tuple[Path, Path]], settings: Settings) -> List[str]:
    if not pairs:
        raise ValueError("No source/destination pairs provided")

    cmd = [settings.jput_bin]
    if settings.jput_n_flag:
        cmd.append("-n")

    if len(pairs) == 1:
        src, dst_dir = pairs[0]
        cmd.extend([str(src), with_trailing_slash(dst_dir)])
        return cmd

    cmd.append("-m")
    for src, dst_dir in pairs:
        cmd.extend([str(src), with_trailing_slash(dst_dir)])
    return cmd



def with_trailing_slash(path: Path) -> str:
    text = str(path)
    return text if text.endswith("/") else text + "/"



def batched(seq: Sequence[PlannedPut], batch_size: int) -> Iterator[Sequence[PlannedPut]]:
    for start in range(0, len(seq), batch_size):
        yield seq[start:start + batch_size]



def format_bytes(num_bytes: int) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    value = float(num_bytes)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.2f} {unit}"
        value /= 1024
    return f"{num_bytes} B"



def print_plan(plans: Sequence[PlannedPut], settings: Settings) -> None:
    direct_count = sum(1 for plan in plans if plan.reason == "direct>=threshold")
    tar_count = len(plans) - direct_count
    tar_member_count = sum(len(plan.members) for plan in plans if plan.reason != "direct>=threshold")

    print("\nPlanned submissions")
    print("=" * 80)
    print(f"Direct files : {direct_count}")
    print(f"Tar archives : {tar_count}")
    print(f"Tar members  : {tar_member_count}")
    print(f"Threshold    : {format_bytes(settings.min_size_bytes)}")
    print(f"Tar target   : {format_bytes(settings.target_archive_size_bytes)}")
    print(f"Tar max      : {format_bytes(settings.max_archive_size_bytes)}")
    print()

    for idx, plan in enumerate(plans, start=1):
        existing = maybe_check_existing_stub(plan)
        print(f"[{idx:03d}] {plan.reason}")
        print(f"  local      : {plan.local_path}")
        print(f"  mss dir    : {plan.mss_dir}")
        print(f"  payload    : {format_bytes(plan.size_bytes)}")
        if existing is not None:
            print(f"  warning    : stub already exists locally: {existing}")
        if plan.reason != "direct>=threshold":
            member_sizes = sum(member.size_bytes for member in plan.members)
            print(f"  members    : {len(plan.members)} files, {format_bytes(member_sizes)} total")
            if member_sizes < settings.min_size_bytes:
                print("  note       : final tarball is below threshold because there were not enough small files to aggregate further")
        print()



def run_command(cmd: Sequence[str]) -> int:
    print("Executing:")
    print("  " + " ".join(cmd))
    print()
    proc = subprocess.Popen(
        list(cmd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        print(line.rstrip())
    return proc.wait()



def submit(plans: Sequence[PlannedPut], settings: Settings) -> int:
    staged_archives: List[Path] = []
    exit_code = 0
    try:
        for plan in plans:
            if plan.reason != "direct>=threshold":
                create_tarball(plan)
                staged_archives.append(plan.local_path)

        for batch in batched(plans, settings.max_pairs_per_request):
            pairs = [(plan.local_path, plan.mss_dir) for plan in batch]
            cmd = build_jput_command(pairs, settings)
            rc = run_command(cmd)
            if rc != 0:
                exit_code = rc
                print(f"ERROR: jput returned non-zero exit status {rc}", file=sys.stderr)
                break
    finally:
        if settings.cleanup_archives:
            for path in staged_archives:
                try:
                    if path.exists():
                        path.unlink()
                except OSError as exc:
                    print(f"WARN: failed to remove staged archive {path}: {exc}", file=sys.stderr)
    return exit_code



def validate_environment(settings: Settings) -> List[str]:
    warnings: List[str] = []
    if not shutil.which(Path(settings.jput_bin).name) and not Path(settings.jput_bin).exists():
        warnings.append(f"jput binary not found: {settings.jput_bin}")
    home = Path.home()
    probable_cert_paths = [
        home / ".globus" / "usercert.pem",
        home / ".globus" / "userkey.pem",
        home / ".jlab" / "certificate.pem",
    ]
    if not any(path.exists() for path in probable_cert_paths):
        warnings.append(
            "No obvious certificate file found under ~/.globus or ~/.jlab; verify your JLab scientific computing certificate manually."
        )
    return warnings



def main() -> int:
    cli = parse_args()
    manifest = load_manifest(Path(cli.manifest).expanduser())
    settings = resolve_settings(manifest, cli)
    jobs = resolve_jobs(manifest)

    warnings = validate_environment(settings)
    if warnings:
        print("Environment warnings")
        print("-" * 80)
        for msg in warnings:
            print(f"* {msg}")
        print()

    all_plans: List[PlannedPut] = []
    for job in jobs:
        all_plans.extend(plan_job(job, settings))

    print_plan(all_plans, settings)

    if not settings.submit:
        print("Dry run only. No jput requests were submitted.")
        print("Use --submit to actually stage tarballs and invoke jput.")
        return 0

    return submit(all_plans, settings)


if __name__ == "__main__":
    raise SystemExit(main())
