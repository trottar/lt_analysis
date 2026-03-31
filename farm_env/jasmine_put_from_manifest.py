#!/usr/bin/env python3
"""
Prepare files for Jefferson Lab Jasmine tape storage from a JSON manifest.

Grounding:
- This script is based on the user's supplied Jasmine procedure:
  * files written to tape appear in the /mss stub tree,
  * users submit write requests via jput,
  * small files should be packaged into .tar or .zip archives,
  * files smaller than 200 MB should be avoided when possible,
  * optimal processing sizes are typically larger than that.

Conservative design choices:
- Small files are grouped into tar archives before submission.
- Submission is one planned output per jput invocation by default. This avoids
  assuming undocumented multi-file CLI behavior.
- The exact jput command shape is configurable from the manifest if your local
  environment requires a wrapper or different argument order.
- Dry-run is the default.
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
import shutil
import subprocess
import sys
import tarfile
from dataclasses import dataclass, field
from os.path import expandvars
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence

try:
    from ltsep_paths import LtsepPaths, resolve_ltsep_paths
except ModuleNotFoundError:
    from farm_env.ltsep_paths import LtsepPaths, resolve_ltsep_paths

MB = 1024 ** 2
GB = 1024 ** 3
DEFAULT_MIN_SIZE_MB = 200
DEFAULT_TARGET_ARCHIVE_SIZE_GB = 4
DEFAULT_MAX_ARCHIVE_SIZE_GB = 20
DEFAULT_STAGE_DIR = "/scratch/$USER/jasmine_stage"
DEFAULT_JPUT = shutil.which("jput") or "jput"
DEFAULT_RUN_MATCH_REGEX_TEMPLATE = r"(^|[^0-9]){run}([^0-9]|$)"
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_MANIFEST_DIR = REPO_ROOT / "input" / "kaon"
VALID_PHI = {"center", "left", "right"}
VALID_EPSILON = {"high", "low"}
VALID_TARGET = {"lh2", "data", "dummy"}
VALID_PRODUCT_KINDS = {"replay", "skim"}


def render_run_template(value, run):
    if value is None:
        return None
    return expandvars(str(value)).format(run=int(run))


def read_runs_file(path_value) -> List[int]:
    path = Path(expandvars(str(path_value))).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Runs file does not exist: {path}")

    runs: List[int] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            text = raw.strip()
            if not text:
                continue
            if text.startswith("#") or text.startswith(">"):
                continue
            try:
                runs.append(int(text))
            except ValueError:
                # Ignore non-integer header lines such as:
                # Q3p0W2p32center_highe
                continue

    if not runs:
        raise ValueError(f"No run numbers found in runs file: {path}")

    return runs


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

    @property
    def destination_stub(self) -> Path:
        return self.mss_dir / self.local_path.name


@dataclass
class JobConfig:
    source: Path
    destination: Path
    skim_destination: Optional[Path] = None
    recursive: bool = True
    archive_prefix: Optional[str] = None
    preserve_tree: bool = True
    tar_subdir: Optional[str] = None
    tar_destination: Optional[Path] = None
    run: Optional[int] = None
    run_match_regex: Optional[str] = None


@dataclass
class SubmissionConfig:
    jput_bin: str
    command_mode: str
    command_template: Optional[List[str]]
    extra_args: List[str]


@dataclass
class Settings:
    min_size_bytes: int
    target_archive_size_bytes: int
    max_archive_size_bytes: int
    stage_dir: Path
    product_kind: str
    submit: bool
    cleanup_archives: bool
    submission: SubmissionConfig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select a kaon manifest from input/kaon using phi/epsilon/Q2/W/target, "
            "package sub-threshold files into tar archives, and optionally invoke "
            "jput for /mss tape storage."
        ),
        epilog=(
            "Examples:\n"
            "  ./jasmine_put_from_manifest.py center high 3p0 3p14 lh2 --product-kind replay\n"
            "  ./jasmine_put_from_manifest.py left low 4p4 2p74 dummy --product-kind skim --submit"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("phi", nargs="?", help="Phi setting: center, left, or right")
    parser.add_argument("epsilon", nargs="?", help="Epsilon setting: high or low")
    parser.add_argument("q2", nargs="?", help="Q2 token, for example 3p0 or Q3p0")
    parser.add_argument("w", nargs="?", help="W token, for example 3p14 or W3p14")
    parser.add_argument("target", nargs="?", help="Target type: lh2/data or dummy")
    parser.add_argument(
        "--manifest-path",
        default=None,
        help="Exact manifest path to use. Overrides phi/epsilon/Q2/W/target selection.",
    )
    parser.add_argument(
        "--manifest-dir",
        default=str(DEFAULT_MANIFEST_DIR),
        help=f"Directory containing KaonLT manifests (default: {DEFAULT_MANIFEST_DIR})",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="Actually invoke submission commands. Default is dry-run.",
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
    parser.add_argument(
        "--product-kind",
        choices=sorted(VALID_PRODUCT_KINDS),
        default="replay",
        help="Which product family to upload: replay or skim (default: replay)",
    )
    parser.add_argument(
        "--run",
        type=int,
        default=None,
        help="Restrict processing to a single run contained in the selected manifest.",
    )
    return parser.parse_args()


def normalize_named_token(name: str, value: str) -> str:
    text = value.strip().lower()
    if not text:
        raise ValueError(f"{name} must not be empty")
    return text


def normalize_kinematics_token(name: str, value: str, prefix: str) -> str:
    text = value.strip().lower()
    if not text:
        raise ValueError(f"{name} must not be empty")
    if text.startswith(prefix):
        text = text[1:]
    return text


def build_manifest_stem(phi: str, epsilon: str, q2: str, w: str, target: str) -> str:
    stem = f"Q{q2}W{w}{phi}_{epsilon}e"
    if target == "dummy":
        stem += "_dummy"
    return stem


def resolve_manifest_path(cli: argparse.Namespace) -> Path:
    if cli.manifest_path is not None:
        manifest_path = Path(expandvars(str(cli.manifest_path))).expanduser()
        if not manifest_path.exists():
            raise FileNotFoundError(f"Manifest path does not exist: {manifest_path}")
        return manifest_path

    if not all([cli.phi, cli.epsilon, cli.q2, cli.w, cli.target]):
        raise ValueError(
            "Provide either --manifest-path or all of phi, epsilon, Q2, W, and target."
        )

    phi = normalize_named_token("phi", cli.phi)
    epsilon = normalize_named_token("epsilon", cli.epsilon)
    q2 = normalize_kinematics_token("Q2", cli.q2, "q")
    w = normalize_kinematics_token("W", cli.w, "w")
    target = normalize_named_token("target", cli.target)

    if phi not in VALID_PHI:
        raise ValueError(f"phi must be one of: {', '.join(sorted(VALID_PHI))}")
    if epsilon not in VALID_EPSILON:
        raise ValueError(f"epsilon must be one of: {', '.join(sorted(VALID_EPSILON))}")
    if target not in VALID_TARGET:
        raise ValueError(f"target must be one of: {', '.join(sorted(VALID_TARGET))}")
    if epsilon == "low" and phi == "right":
        raise ValueError("No low-epsilon right-setting manifest exists in input/kaon.")

    manifest_dir = Path(expandvars(str(cli.manifest_dir))).expanduser()
    manifest_stem = build_manifest_stem(phi, epsilon, q2, w, target)
    manifest_path = manifest_dir / f"{manifest_stem}.json"

    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Manifest not found for phi={phi}, epsilon={epsilon}, Q2={q2}, W={w}, "
            f"target={target}: {manifest_path}"
        )

    return manifest_path



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



def resolve_submission(defaults: Dict) -> SubmissionConfig:
    submit_cfg = defaults.get("submit", {})
    if not isinstance(submit_cfg, dict):
        raise ValueError("defaults.submit must be an object when provided")

    command_mode = str(submit_cfg.get("command_mode", "src_dest")).strip()
    if command_mode not in {"src_dest", "template"}:
        raise ValueError("defaults.submit.command_mode must be 'src_dest' or 'template'")

    command_template = submit_cfg.get("command_template")
    if command_mode == "template":
        if not isinstance(command_template, list) or not command_template:
            raise ValueError(
                "defaults.submit.command_template must be a non-empty list when command_mode='template'"
            )
        command_template = [str(x) for x in command_template]
    else:
        command_template = None

    extra_args = submit_cfg.get("extra_args", [])
    if not isinstance(extra_args, list):
        raise ValueError("defaults.submit.extra_args must be a list when provided")

    return SubmissionConfig(
        jput_bin=str(submit_cfg.get("jput_bin", DEFAULT_JPUT)),
        command_mode=command_mode,
        command_template=command_template,
        extra_args=[str(x) for x in extra_args],
    )



def resolve_settings(manifest: Dict, cli: argparse.Namespace) -> Settings:
    defaults = manifest.get("defaults", {})
    min_size_mb = float(defaults.get("min_size_mb", DEFAULT_MIN_SIZE_MB))
    target_archive_size_gb = float(defaults.get("target_archive_size_gb", DEFAULT_TARGET_ARCHIVE_SIZE_GB))
    max_archive_size_gb = float(defaults.get("max_archive_size_gb", DEFAULT_MAX_ARCHIVE_SIZE_GB))

    if min_size_mb <= 0:
        raise ValueError("min_size_mb must be > 0")
    if target_archive_size_gb <= 0:
        raise ValueError("target_archive_size_gb must be > 0")
    if max_archive_size_gb <= 0:
        raise ValueError("max_archive_size_gb must be > 0")
    if target_archive_size_gb > max_archive_size_gb:
        raise ValueError("target_archive_size_gb cannot exceed max_archive_size_gb")

    stage_dir = Path(expandvars(cli.stage_dir or defaults.get("stage_dir", DEFAULT_STAGE_DIR))).expanduser()
    submission = resolve_submission(defaults)

    return Settings(
        min_size_bytes=int(min_size_mb * MB),
        target_archive_size_bytes=int(target_archive_size_gb * GB),
        max_archive_size_bytes=int(max_archive_size_gb * GB),
        stage_dir=stage_dir,
        product_kind=str(cli.product_kind),
        submit=cli.submit,
        cleanup_archives=not cli.keep_archives,
        submission=submission,
    )



def resolve_product_source_root(settings: Settings, paths: LtsepPaths) -> Path:
    if settings.product_kind == "replay":
        return paths.replay_source_dir
    return paths.skim_source_dir



def resolve_jobs(
    manifest: Dict, settings: Settings, paths: LtsepPaths, cli: argparse.Namespace
) -> List[JobConfig]:
    jobs = manifest.get("jobs")
    if not isinstance(jobs, list) or not jobs:
        raise ValueError("Manifest must contain a non-empty 'jobs' list")

    resolved: List[JobConfig] = []
    source_root = resolve_product_source_root(settings, paths).resolve()

    for index, entry in enumerate(jobs, start=1):
        if not isinstance(entry, dict):
            raise ValueError(f"jobs[{index}] must be an object")

        recursive = bool_from_manifest(entry.get("recursive"), True)
        preserve_tree = bool_from_manifest(entry.get("preserve_tree"), True)
        runs = entry.get("runs")
        runs_file = entry.get("runs_file")

        if runs is not None and runs_file is not None:
            raise ValueError(f"jobs[{index}] cannot define both 'runs' and 'runs_file'")

        archive_prefix_default = "FullReplay" if settings.product_kind == "replay" else "ApplyCuts"

        # Single fixed job with no run expansion
        if runs is None and runs_file is None:
            destination_value = entry.get("destination")
            if destination_value is None:
                raise ValueError(f"jobs[{index}] missing 'destination'")
            source = source_root
            destination = Path(expandvars(str(destination_value))).expanduser()

            tar_destination_value = entry.get("tar_destination")
            tar_destination = None
            if tar_destination_value is not None:
                tar_destination = Path(expandvars(str(tar_destination_value))).expanduser()
            skim_destination_value = entry.get("skim_destination")
            skim_destination = None
            if skim_destination_value is not None:
                skim_destination = Path(expandvars(str(skim_destination_value))).expanduser()

            resolved.append(
                JobConfig(
                    source=source,
                    destination=destination,
                    skim_destination=skim_destination,
                    recursive=recursive,
                    archive_prefix=entry.get("archive_prefix") or f"{archive_prefix_default}_{paths.anatype}LT",
                    preserve_tree=preserve_tree,
                    tar_subdir=entry.get("tar_subdir"),
                    tar_destination=tar_destination,
                )
            )
            continue

        if runs_file is not None:
            runs = read_runs_file(runs_file)

        if not isinstance(runs, list) or not runs:
            raise ValueError(f"jobs[{index}] 'runs' must be a non-empty list")

        if cli.run is not None:
            runs = [run for run in runs if int(run) == int(cli.run)]
            if not runs:
                continue

        destination_template = entry.get("destination_template", entry.get("destination"))
        archive_prefix_template = entry.get("archive_prefix_template", entry.get("archive_prefix"))
        skim_archive_prefix_template = entry.get("skim_archive_prefix_template")
        tar_destination_template = entry.get("tar_destination_template", entry.get("tar_destination"))
        skim_destination_template = entry.get("skim_destination_template", entry.get("skim_destination"))
        run_match_regex_template = entry.get("run_match_regex_template")

        if destination_template is None:
            raise ValueError(f"jobs[{index}] missing 'destination_template' (or 'destination')")
        if run_match_regex_template is None:
            run_match_regex_template = DEFAULT_RUN_MATCH_REGEX_TEMPLATE

        for run in runs:
            source = source_root
            destination = Path(render_run_template(destination_template, run)).expanduser()
            if settings.product_kind == "skim":
                if skim_archive_prefix_template is not None:
                    archive_prefix = render_run_template(skim_archive_prefix_template, run)
                else:
                    archive_prefix = f"ApplyCuts_run_{int(run)}_{paths.anatype}LT"
            elif archive_prefix_template is not None:
                archive_prefix = render_run_template(archive_prefix_template, run)
            else:
                archive_prefix = f"{archive_prefix_default}_run_{int(run)}_{paths.anatype}LT"
            tar_destination = None
            if tar_destination_template is not None:
                tar_destination = Path(render_run_template(tar_destination_template, run)).expanduser()
            skim_destination = None
            if skim_destination_template is not None:
                skim_destination = Path(render_run_template(skim_destination_template, run)).expanduser()
            run_match_regex = None
            if run_match_regex_template is not None:
                run_match_regex = render_run_template(run_match_regex_template, run)

            resolved.append(
                JobConfig(
                    source=source,
                    destination=destination,
                    skim_destination=skim_destination,
                    recursive=recursive,
                    archive_prefix=archive_prefix,
                    preserve_tree=preserve_tree,
                    tar_subdir=render_run_template(entry.get("tar_subdir"), run),
                    tar_destination=tar_destination,
                    run=int(run),
                    run_match_regex=run_match_regex,
                )
            )

    if cli.run is not None and not resolved:
        raise ValueError(f"Run {cli.run} was not found in the selected manifest.")

    return resolved



def ensure_supported_destination(destination: Path) -> None:
    if not str(destination).startswith("/mss/"):
        raise ValueError(f"Destination must be under /mss/: {destination}")


def resolve_job_destination(job: JobConfig, settings: Settings) -> Path:
    if settings.product_kind == "skim":
        if job.skim_destination is None:
            raise ValueError(
                f"Manifest job for source {job.source} is missing 'skim_destination' "
                f"but --product-kind=skim was requested."
            )
        return job.skim_destination
    return job.destination



def iter_source_files(job: JobConfig) -> Iterator[SourceFile]:
    source = job.source
    if not source.exists():
        raise FileNotFoundError(f"Source path does not exist: {source}")

    if source.is_file():
        name_text = source.name
        if job.run_match_regex is not None and re.search(job.run_match_regex, name_text) is None:
            return
        yield SourceFile(path=source, size_bytes=source.stat().st_size, arcname=source.name)
        return

    if not source.is_dir():
        raise ValueError(f"Unsupported source type: {source}")

    if job.recursive:
        children = sorted(p for p in source.rglob("*") if p.is_file())
    else:
        children = sorted(p for p in source.iterdir() if p.is_file())

    for child in children:
        rel_text = child.relative_to(source).as_posix()
        name_text = child.name

        if job.run_match_regex is not None and re.search(job.run_match_regex, name_text) is None:
            continue

        arcname = rel_text if job.preserve_tree else child.name
        yield SourceFile(path=child, size_bytes=child.stat().st_size, arcname=arcname)



def chunk_small_files(files: Sequence[SourceFile], target_bytes: int, max_bytes: int) -> List[List[SourceFile]]:
    """Simple greedy packing for small-file tar archives."""
    chunks: List[List[SourceFile]] = []
    chunk_sizes: List[int] = []

    for item in sorted(files, key=lambda f: f.size_bytes, reverse=True):
        placed = False

        preferred = [
            i for i, size in enumerate(chunk_sizes)
            if size < target_bytes and size + item.size_bytes <= max_bytes
        ]
        preferred.sort(key=lambda i: chunk_sizes[i], reverse=True)
        for idx in preferred:
            chunks[idx].append(item)
            chunk_sizes[idx] += item.size_bytes
            placed = True
            break

        if placed:
            continue

        fallback = [i for i, size in enumerate(chunk_sizes) if size + item.size_bytes <= max_bytes]
        fallback.sort(key=lambda i: chunk_sizes[i], reverse=True)
        for idx in fallback:
            chunks[idx].append(item)
            chunk_sizes[idx] += item.size_bytes
            placed = True
            break

        if not placed:
            chunks.append([item])
            chunk_sizes.append(item.size_bytes)

    return chunks



def sanitize_name(name: str) -> str:
    safe = [c if c.isalnum() or c in {"-", "_", "."} else "_" for c in name]
    text = "".join(safe).strip("._")
    return text or "archive"



def resolve_tar_destination(job: JobConfig) -> Path:
    if job.tar_destination is not None:
        return job.tar_destination
    if job.tar_subdir:
        return job.destination.parent / job.tar_subdir
    return job.destination



def plan_job(job: JobConfig, settings: Settings) -> List[PlannedPut]:
    job_destination = resolve_job_destination(job, settings)
    ensure_supported_destination(job_destination)
    tar_mss_dir = resolve_tar_destination(job) if settings.product_kind == "replay" else job_destination
    ensure_supported_destination(tar_mss_dir)

    files = list(iter_source_files(job))
    if not files:
        if job.run is not None:
            print(f"WARN: no files matched for run {job.run} under {job.source}", file=sys.stderr)
            return []
        raise ValueError(f"No files found under source: {job.source}")

    direct_files = [f for f in files if f.size_bytes >= settings.min_size_bytes]
    small_files = [f for f in files if f.size_bytes < settings.min_size_bytes]

    planned: List[PlannedPut] = []
    for entry in direct_files:
        planned.append(
            PlannedPut(
                local_path=entry.path,
                mss_dir=job_destination,
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
        prefix = job.archive_prefix or sanitize_name(job.source.name or "archive")

        for idx, group in enumerate(groups, start=1):
            planned.append(
                PlannedPut(
                    local_path=settings.stage_dir / f"{prefix}_{idx:04d}.tar",
                    mss_dir=tar_mss_dir,
                    reason="tarball<200MB_inputs",
                    size_bytes=sum(member.size_bytes for member in group),
                    members=group,
                )
            )

    return planned



def create_tarball(plan: PlannedPut) -> None:
    plan.local_path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(plan.local_path, mode="w") as tar:
        for member in plan.members:
            tar.add(member.path, arcname=member.arcname, recursive=False)



def format_bytes(num_bytes: int) -> str:
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    value = float(num_bytes)
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return f"{num_bytes} B"



def maybe_check_existing_stub(plan: PlannedPut) -> bool:
    return plan.destination_stub.exists()



def build_submit_command(plan: PlannedPut, settings: Settings) -> List[str]:
    sub = settings.submission
    if sub.command_mode == "src_dest":
        return [sub.jput_bin, *sub.extra_args, str(plan.local_path), with_trailing_slash(plan.mss_dir)]

    if sub.command_mode == "template":
        replacements = {
            "{jput_bin}": sub.jput_bin,
            "{src}": str(plan.local_path),
            "{dest}": with_trailing_slash(plan.mss_dir),
            "{dest_stub}": str(plan.destination_stub),
            "{name}": plan.local_path.name,
        }
        cmd: List[str] = []
        assert sub.command_template is not None
        for token in sub.command_template:
            expanded = str(token)
            for key, value in replacements.items():
                expanded = expanded.replace(key, value)
            cmd.append(expanded)
        return cmd

    raise ValueError(f"Unsupported command mode: {sub.command_mode}")



def with_trailing_slash(path: Path) -> str:
    text = str(path)
    return text if text.endswith("/") else text + "/"



def print_plan(plans: Sequence[PlannedPut], settings: Settings) -> None:
    direct_count = sum(1 for p in plans if p.reason == "direct>=threshold")
    tar_count = len(plans) - direct_count
    tar_members = sum(len(p.members) for p in plans if p.reason != "direct>=threshold")

    print("\nPlanned submissions")
    print("=" * 80)
    print(f"Direct files : {direct_count}")
    print(f"Tar archives : {tar_count}")
    print(f"Tar members  : {tar_members}")
    print(f"Threshold    : {format_bytes(settings.min_size_bytes)}")
    print(f"Tar target   : {format_bytes(settings.target_archive_size_bytes)}")
    print(f"Tar max      : {format_bytes(settings.max_archive_size_bytes)}")
    print(f"Submit mode  : {settings.submission.command_mode}")
    print()

    for idx, plan in enumerate(plans, start=1):
        print(f"[{idx:03d}] {plan.reason}")
        print(f"  local        : {plan.local_path}")
        print(f"  mss dir      : {plan.mss_dir}")
        print(f"  destination  : {plan.destination_stub}")
        print(f"  payload      : {format_bytes(plan.size_bytes)}")
        if maybe_check_existing_stub(plan):
            print("  warning      : destination stub already exists locally")
        if plan.reason != "direct>=threshold":
            member_total = sum(member.size_bytes for member in plan.members)
            print(f"  members      : {len(plan.members)} files, {format_bytes(member_total)} total")
            if member_total < settings.min_size_bytes:
                print("  note         : final tarball is still below threshold because not enough small files were available to combine further")
        cmd_preview = build_submit_command(plan, settings)
        print("  submit cmd   : " + " ".join(shlex.quote(x) for x in cmd_preview))
        print()



def run_command(cmd: Sequence[str]) -> int:
    print("Executing:")
    print("  " + " ".join(shlex.quote(x) for x in cmd))
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

        for plan in plans:
            cmd = build_submit_command(plan, settings)
            rc = run_command(cmd)
            if rc != 0:
                exit_code = rc
                print(f"ERROR: submission command returned non-zero exit status {rc}", file=sys.stderr)
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
    sub = settings.submission
    if shutil.which(Path(sub.jput_bin).name) is None and not Path(sub.jput_bin).exists():
        warnings.append(f"submission binary not found in PATH and not found as a file: {sub.jput_bin}")

    home = Path.home()
    probable_cert_paths = [
        home / ".globus" / "usercert.pem",
        home / ".globus" / "userkey.pem",
        home / ".jlab" / "certificate.pem",
    ]
    if not any(path.exists() for path in probable_cert_paths):
        warnings.append(
            "No obvious certificate file found under ~/.globus or ~/.jlab; verify certificate setup manually before submission."
        )
    return warnings



def main() -> int:
    cli = parse_args()
    try:
        manifest_path = resolve_manifest_path(cli)
        manifest = load_manifest(manifest_path)
    except (FileNotFoundError, ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    settings = resolve_settings(manifest, cli)
    paths = resolve_ltsep_paths(__file__)
    jobs = resolve_jobs(manifest, settings, paths, cli)

    print(f"Using manifest: {manifest_path}")
    print(f"Product kind : {settings.product_kind}")
    print(f"Source root  : {resolve_product_source_root(settings, paths)}")

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
        print("Dry run only. No submission commands were executed.")
        print("Use --submit to create tarballs and invoke the configured submission command.")
        return 0

    return submit(all_plans, settings)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
