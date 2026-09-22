#!/usr/bin/env python3
"""Print a narrow KaonLT memory startup observation without active-state synthesis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Any, Callable


CORE_RECORDS = (
    "docs/memory/AGENTS.md",
    "docs/memory/CURRENT.md",
    "docs/memory/MEMORY.md",
    "docs/memory/handoffs/CURRENT_HANDOFF.md",
    "docs/memory/USER.md",
)
STABLE_HANDOFF = "No exceptional transfer state is recorded."
REFERENCE = re.compile(r"(?<![A-Za-z0-9_])(?:(?:docs/memory/)?[A-Za-z0-9_.-]+/)*[A-Za-z0-9_.-]+\.md")


def default_root() -> Path:
    return Path(__file__).resolve().parents[1]


def git_command(root: Path, *arguments: str) -> tuple[bool, str | None]:
    try:
        result = subprocess.run(
            ["git", *arguments], cwd=root, check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.CalledProcessError):
        return False, None
    return True, result.stdout


def observe_git(root: Path) -> dict[str, Any]:
    """Collect runtime Git facts, preserving clean output versus unavailability."""
    branch_available, branch = git_command(root, "branch", "--show-current")
    head_available, head = git_command(root, "rev-parse", "HEAD")
    status_available, status = git_command(root, "status", "--short", "--untracked-files=all")
    status_short = status.splitlines() if status_available and status is not None else None
    return {
        "branch": branch.strip() if branch_available and branch is not None else None,
        "head": head.strip() if head_available and head is not None else None,
        "worktree": {
            "available": status_available,
            "clean": not status_short if status_short is not None else None,
            "status_short": status_short,
        },
    }


def direct_memory_references(root: Path, current: Path) -> list[str]:
    if not current.is_file():
        return []
    memory_root = root / "docs" / "memory"
    found: set[str] = set()
    for match in REFERENCE.findall(current.read_text(encoding="utf-8")):
        relative = Path(match)
        candidate = root / relative if match.startswith("docs/memory/") else memory_root / relative
        try:
            resolved = candidate.resolve()
            resolved.relative_to(memory_root.resolve())
        except ValueError:
            continue
        if resolved.is_file():
            found.add(resolved.relative_to(root).as_posix())
    return sorted(found)


def run_health(root: Path) -> tuple[int, str]:
    script = root / "tools" / "check_memory_health.py"
    if not script.is_file():
        return 1, "MEMORY HEALTH: FAIL: health checker is missing"
    result = subprocess.run(
        [sys.executable, str(script), "--root", str(root)],
        cwd=root,
        capture_output=True,
        text=True,
    )
    return result.returncode, (result.stdout + result.stderr).strip()


def health_status(returncode: int, output: str) -> str:
    if returncode:
        return "fail"
    if "MEMORY HEALTH: WARN" in output:
        return "warn"
    return "pass"


def exceptional_handoff(path: Path) -> dict[str, bool | str | None]:
    """Report only whether canonical handoff transfer state is exceptional."""
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        return {"present": None, "error": f"cannot read handoff: {error}"}
    match = re.search(r"^## Transfer State\s*$\n(?P<body>.*?)(?=^##\s|\Z)", text, re.MULTILINE | re.DOTALL)
    if match is None:
        return {"present": None, "error": "cannot locate Transfer State"}
    transfer = match.group("body").strip()
    if not transfer:
        return {"present": None, "error": "Transfer State is empty"}
    return {"present": transfer != STABLE_HANDOFF}


def collect_summary(
    root: Path,
    health_runner: Callable[[Path], tuple[int, str]] = run_health,
    git_observer: Callable[[Path], dict[str, Any]] = observe_git,
) -> dict[str, Any]:
    core_records = []
    for relative in CORE_RECORDS:
        path = root / relative
        core_records.append({"path": relative, "bytes": path.stat().st_size if path.is_file() else None})
    current = root / CORE_RECORDS[1]
    returncode, health_output = health_runner(root)
    return {
        "core_records": core_records,
        "current_references": direct_memory_references(root, current),
        "exceptional_handoff": exceptional_handoff(root / CORE_RECORDS[3]),
        "git": git_observer(root),
        "memory_health": {"status": health_status(returncode, health_output), "returncode": returncode},
    }


def worktree_label(worktree: dict[str, Any]) -> str:
    if not worktree.get("available"):
        return "UNAVAILABLE"
    return "CLEAN" if worktree.get("clean") else "DIRTY"


def handoff_label(handoff: dict[str, bool | str | None]) -> str:
    present = handoff.get("present")
    if present is None:
        return "UNAVAILABLE"
    return "YES" if present else "NO"


def print_human(summary: dict[str, Any]) -> None:
    git = summary["git"]
    print("KaonLT memory bootstrap")
    print(f"git branch: {git['branch'] or 'unavailable'}")
    print(f"git HEAD: {git['head'] or 'unavailable'}")
    print(f"worktree: {worktree_label(git['worktree'])}")
    print("startup core:")
    for record in summary["core_records"]:
        size = "missing" if record["bytes"] is None else str(record["bytes"])
        print(f"  {record['path']}: {size} bytes")
    print("CURRENT.md references:")
    for reference in summary["current_references"]:
        print(f"  {reference}")
    if not summary["current_references"]:
        print("  none resolved")
    print(f"exceptional handoff: {handoff_label(summary['exceptional_handoff'])}")
    print(f"memory health: {summary['memory_health']['status'].upper()}")


def fixture_current() -> str:
    return "# Current\n\nSee evidence/example.md and nowhere else.\n"


def fixture_git(*, clean: bool, status_short: list[str] | None = None) -> dict[str, Any]:
    return {
        "branch": "fixture",
        "head": "a" * 40,
        "worktree": {
            "available": True,
            "clean": clean,
            "status_short": [] if status_short is None else status_short,
        },
    }


def self_test() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        memory = root / "docs" / "memory"
        evidence = memory / "evidence"
        evidence.mkdir(parents=True)
        (evidence / "example.md").write_text("evidence\n", encoding="utf-8")
        for relative, text in (
            ("AGENTS.md", "agents\n"),
            ("CURRENT.md", fixture_current()),
            ("MEMORY.md", "memory\n"),
            ("handoffs/CURRENT_HANDOFF.md", "# Handoff\n\n## Transfer State\n\nNo exceptional transfer state is recorded.\n\n## Resume\n\nCURRENT.md\n"),
            ("USER.md", "user\n"),
        ):
            path = memory / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8")
        summary = collect_summary(
            root,
            lambda _: (0, "MEMORY HEALTH: PASS"),
            lambda _: fixture_git(clean=True),
        )
        assert [record["path"] for record in summary["core_records"]] == list(CORE_RECORDS), summary
        assert all(isinstance(record["bytes"], int) for record in summary["core_records"]), summary
        assert summary["current_references"] == ["docs/memory/evidence/example.md"], summary
        assert summary["exceptional_handoff"] == {"present": False}, summary
        assert summary["git"]["worktree"] == {"available": True, "clean": True, "status_short": []}, summary
        assert "active_state" not in summary, summary
        assert summary["memory_health"]["status"] == "pass", summary

        dirty = collect_summary(
            root,
            lambda _: (0, "MEMORY HEALTH: WARN: fixture"),
            lambda _: fixture_git(clean=False, status_short=[" M docs/memory/CURRENT.md", "?? scratch.txt"]),
        )
        assert dirty["git"]["worktree"]["clean"] is False, dirty
        assert dirty["git"]["worktree"]["status_short"] == [" M docs/memory/CURRENT.md", "?? scratch.txt"], dirty
        assert dirty["memory_health"]["status"] == "warn", dirty
        (memory / "handoffs/CURRENT_HANDOFF.md").write_text("# Handoff\n\n## Transfer State\n\nExceptional transfer.\n\n## Resume\n\nCURRENT.md\n", encoding="utf-8")
        exceptional = collect_summary(root, lambda _: (1, "MEMORY HEALTH: FAIL"), lambda _: fixture_git(clean=True))
        assert exceptional["exceptional_handoff"] == {"present": True}, exceptional
        assert exceptional["memory_health"]["status"] == "fail", exceptional


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root(), help="repository root")
    parser.add_argument("--json", action="store_true", help="emit machine-readable summary")
    parser.add_argument("--self-test", action="store_true", help="run deterministic temporary-fixture checks")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        self_test()
        print("SELF-TEST: PASS")
        return 0
    summary = collect_summary(args.root.resolve())
    if args.json:
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print_human(summary)
    return 0 if summary["memory_health"]["returncode"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
