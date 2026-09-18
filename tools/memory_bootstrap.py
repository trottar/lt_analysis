#!/usr/bin/env python3
"""Print a narrow, non-authoritative KaonLT memory startup summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Any, Callable

import update_memory_manifest as manifest_tool


ACTIVE_FILES = (
    "docs/memory/CURRENT.md",
    "docs/memory/MEMORY.md",
    "docs/memory/handoffs/CURRENT_HANDOFF.md",
)
REFERENCE = re.compile(r"(?<![A-Za-z0-9_])(?:(?:docs/memory/)?[A-Za-z0-9_.-]+/)*[A-Za-z0-9_.-]+\.md")


def default_root() -> Path:
    return Path(__file__).resolve().parents[1]


def git_value(root: Path, *arguments: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", *arguments], cwd=root, check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip() or None


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


def active_state(current: Path) -> dict[str, str | int] | None:
    try:
        return manifest_tool.parse_active_state(current)
    except manifest_tool.ActiveStateError:
        return None


def collect_summary(
    root: Path,
    health_runner: Callable[[Path], tuple[int, str]] = run_health,
) -> dict[str, Any]:
    records = []
    for relative in ACTIVE_FILES:
        path = root / relative
        records.append({"path": relative, "bytes": path.stat().st_size if path.is_file() else None})
    current = root / ACTIVE_FILES[0]
    returncode, health_output = health_runner(root)
    return {
        "active_records": records,
        "active_state": active_state(current),
        "current_references": direct_memory_references(root, current),
        "git": {"branch": git_value(root, "branch", "--show-current"), "head": git_value(root, "rev-parse", "HEAD")},
        "memory_health": {"status": health_status(returncode, health_output), "returncode": returncode},
    }


def print_human(summary: dict[str, Any]) -> None:
    print("KaonLT memory bootstrap")
    print(f"git branch: {summary['git']['branch'] or 'unavailable'}")
    print(f"git HEAD: {summary['git']['head'] or 'unavailable'}")
    print("active state:")
    if summary["active_state"] is None:
        print("  unavailable: CURRENT.md frontmatter is invalid")
    else:
        for key, value in summary["active_state"].items():
            print(f"  {key}: {value}")
    print("active records:")
    for record in summary["active_records"]:
        size = "missing" if record["bytes"] is None else str(record["bytes"])
        print(f"  {record['path']}: {size} bytes")
    print("CURRENT.md references:")
    for reference in summary["current_references"]:
        print(f"  {reference}")
    if not summary["current_references"]:
        print("  none resolved")
    print(f"memory health: {summary['memory_health']['status'].upper()}")


def fixture_current() -> str:
    return """---
memory_schema: 2
active_objective: Fixture objective
current_work_item: Fixture work item
active_status: ACTIVE
next_action: Fixture next action
scientific_source_commit: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
bundle_profile_commit: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
---
# Current

See evidence/example.md and nowhere else.
"""


def self_test() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        evidence = root / "docs" / "memory" / "evidence"
        evidence.mkdir(parents=True)
        (evidence / "example.md").write_text("evidence\n", encoding="utf-8")
        current = root / "docs" / "memory" / "CURRENT.md"
        current.write_text(fixture_current(), encoding="utf-8")
        (root / "docs" / "memory" / "MEMORY.md").write_text("memory\n", encoding="utf-8")
        handoff = root / "docs" / "memory" / "handoffs"
        handoff.mkdir()
        (handoff / "CURRENT_HANDOFF.md").write_text("handoff\n", encoding="utf-8")
        summary = collect_summary(root, lambda _: (0, "MEMORY HEALTH: PASS"))
        assert summary["current_references"] == ["docs/memory/evidence/example.md"], summary
        assert summary["active_state"] is not None, summary
        assert summary["active_state"]["active_status"] == "ACTIVE", summary
        assert summary["memory_health"]["status"] == "pass", summary


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
