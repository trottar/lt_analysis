#!/usr/bin/env python3
"""Check deterministic structural and integrity properties of KaonLT memory."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Iterable


SOFT_LIMIT_BYTES = 16 * 1024
HARD_LIMIT_BYTES = 32 * 1024
REQUIRED_FILES = (
    "docs/memory/README.md",
    "docs/memory/AGENTS.md",
    "docs/memory/MAINTENANCE.md",
    "docs/memory/COMMUNICATION.md",
    "docs/memory/CODEX.md",
    "docs/memory/chats/CHAT_INDEX.md",
    "docs/memory/templates/CODEX_CONTRACT.md",
    "docs/memory/templates/MEMORY_UPDATE.md",
    "docs/memory/templates/RUNTIME_EVIDENCE.md",
    "docs/memory/manifest.json",
    "tools/check_memory_health.py",
    "tools/update_memory_manifest.py",
    "tools/memory_bootstrap.py",
)
ACTIVE_FILES = (
    "docs/memory/CURRENT.md",
    "docs/memory/MEMORY.md",
    "docs/memory/handoffs/CURRENT_HANDOFF.md",
)
REQUIRED_README_LINKS = (
    "AGENTS.md",
    "MAINTENANCE.md",
    "COMMUNICATION.md",
    "CODEX.md",
    "chats/CHAT_INDEX.md",
    "templates/CODEX_CONTRACT.md",
    "templates/MEMORY_UPDATE.md",
    "templates/RUNTIME_EVIDENCE.md",
    "manifest.json",
)
TOP_LEVEL_ROLES = ("README", "CURRENT", "MEMORY", "AGENTS", "MAINTENANCE", "COMMUNICATION", "CODEX")
MARKDOWN_LINK = re.compile(r"(?<!!)\[[^\]]+\]\(([^)#]+)(?:#[^)]*)?\)")
ROLE_SUFFIX = re.compile(r"(?:[._-](?:copy|backup|old|v\d+))+$", re.IGNORECASE)


def default_root() -> Path:
    return Path(__file__).resolve().parents[1]


def markdown_links(readme: Path) -> set[str]:
    links: set[str] = set()
    for target in MARKDOWN_LINK.findall(readme.read_text(encoding="utf-8")):
        target = target.strip()
        if "://" not in target and not target.startswith("mailto:"):
            links.add(target)
    return links


def role_base(stem: str) -> str:
    return ROLE_SUFFIX.sub("", stem).upper()


def run_manifest_check(root: Path) -> tuple[bool, str]:
    script = root / "tools" / "update_memory_manifest.py"
    if not script.is_file():
        return False, "manifest checker is missing"
    result = subprocess.run(
        [sys.executable, str(script), "--root", str(root), "--check"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    message = (result.stdout + result.stderr).strip()
    return result.returncode == 0, message


def run_checks(
    root: Path,
    soft_limit: int,
    hard_limit: int,
    *,
    check_manifest: bool = True,
) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    if soft_limit <= 0 or hard_limit < soft_limit:
        return ["size limits must be positive and hard limit must be at least soft limit"], warnings

    for relative in REQUIRED_FILES:
        if not (root / relative).is_file():
            errors.append(f"missing required file: {relative}")

    for relative in ACTIVE_FILES:
        path = root / relative
        if not path.is_file():
            errors.append(f"missing active file: {relative}")
            continue
        size = path.stat().st_size
        if size > hard_limit:
            errors.append(f"active file exceeds hard limit ({size} > {hard_limit}): {relative}")
        elif size > soft_limit:
            warnings.append(f"active file exceeds soft limit ({size} > {soft_limit}): {relative}")

    readme = root / "docs" / "memory" / "README.md"
    if readme.is_file():
        links = markdown_links(readme)
        for target in REQUIRED_README_LINKS:
            if target not in links:
                errors.append(f"README missing required control link: {target}")
        for target in links:
            target_path = (readme.parent / target).resolve()
            try:
                target_path.relative_to(root.resolve())
            except ValueError:
                continue
            if not target_path.exists():
                errors.append(f"README link target does not exist: {target}")

    memory_root = root / "docs" / "memory"
    if memory_root.is_dir():
        role_files: dict[str, list[str]] = {role: [] for role in TOP_LEVEL_ROLES}
        for path in memory_root.iterdir():
            if path.is_file() and path.suffix.lower() == ".md":
                role = role_base(path.stem)
                if role in role_files:
                    role_files[role].append(path.name)
        for role, names in role_files.items():
            canonical = f"{role}.md"
            if names != [canonical]:
                errors.append(f"duplicate or missing top-level {role} role: {', '.join(sorted(names)) or 'none'}")

    if check_manifest and not errors:
        passed, message = run_manifest_check(root)
        if not passed:
            errors.append(f"manifest freshness check failed: {message or 'no diagnostic'}")
    return errors, warnings


def print_result(errors: Iterable[str], warnings: Iterable[str]) -> int:
    errors = list(errors)
    warnings = list(warnings)
    for warning in warnings:
        print(f"MEMORY HEALTH: WARN: {warning}")
    for error in errors:
        print(f"MEMORY HEALTH: FAIL: {error}")
    if errors:
        return 1
    print("MEMORY HEALTH: WARN" if warnings else "MEMORY HEALTH: PASS")
    return 0


def self_test() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        memory = root / "docs" / "memory"
        (memory / "handoffs").mkdir(parents=True)
        (memory / "chats").mkdir()
        (memory / "templates").mkdir()
        (root / "tools").mkdir()
        readme_links = "\n".join(f"[x]({target})" for target in REQUIRED_README_LINKS)
        (memory / "README.md").write_text(readme_links, encoding="utf-8")
        for name in ("CURRENT.md", "MEMORY.md", "AGENTS.md", "MAINTENANCE.md", "COMMUNICATION.md", "CODEX.md"):
            (memory / name).write_text("ok\n", encoding="utf-8")
        (memory / "handoffs" / "CURRENT_HANDOFF.md").write_text("ok\n", encoding="utf-8")
        (memory / "chats" / "CHAT_INDEX.md").write_text("ok\n", encoding="utf-8")
        for name in ("CODEX_CONTRACT.md", "MEMORY_UPDATE.md", "RUNTIME_EVIDENCE.md"):
            (memory / "templates" / name).write_text("ok\n", encoding="utf-8")
        (memory / "manifest.json").write_text("{}\n", encoding="utf-8")
        for name in ("check_memory_health.py", "update_memory_manifest.py", "memory_bootstrap.py"):
            (root / "tools" / name).write_text("# fixture\n", encoding="utf-8")

        errors, warnings = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert not errors and not warnings, (errors, warnings)
        (memory / "CURRENT_copy.md").write_text("duplicate\n", encoding="utf-8")
        errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert any("top-level CURRENT role" in error for error in errors), errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root(), help="repository root")
    parser.add_argument("--soft-limit-bytes", type=int, default=SOFT_LIMIT_BYTES)
    parser.add_argument("--hard-limit-bytes", type=int, default=HARD_LIMIT_BYTES)
    parser.add_argument("--self-test", action="store_true", help="run deterministic temporary-fixture checks")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        self_test()
        print("SELF-TEST: PASS")
        return 0
    errors, warnings = run_checks(
        args.root.resolve(), args.soft_limit_bytes, args.hard_limit_bytes
    )
    return print_result(errors, warnings)


if __name__ == "__main__":
    raise SystemExit(main())
