#!/usr/bin/env python3
"""Check deterministic structural and integrity properties of KaonLT memory."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Callable, Iterable

import update_memory_manifest as manifest_tool


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
ACTIVE_SURFACES = (
    "docs/memory/CURRENT.md",
    "docs/memory/handoffs/CURRENT_HANDOFF.md",
    "docs/memory/roadmap/CURRENT.md",
)
CURRENT_SECTIONS = (
    "Active Objective",
    "Current Work Item",
    "Verified State",
    "Source / Evidence Identity",
    "Blockers",
    "Next Action",
    "Success Criteria",
    "Do Not Reopen Without New Evidence",
    "Relevant References",
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
HEADING = re.compile(r"^##\s+(.+?)\s*$", re.MULTILINE)
ROLE_SUFFIX = re.compile(r"(?:[._-](?:copy|backup|old|v\d+))+$", re.IGNORECASE)
COMMIT = re.compile(r"^[0-9a-f]{40}$")
AUTHORITY_MARKERS = {
    "docs/memory/AGENTS.md": (
        "## Execution authority",
        "Codex must not commit, push, update remote refs",
        "The user alone commits/pushes accepted changes and runs farm",
        "Workflow: Codex local changes -> ChatGPT audit -> user commit/push",
    ),
    "docs/memory/CODEX.md": (
        "Codex must not commit, push, update remote refs",
        "The user alone commits/pushes accepted changes and runs farm",
        "Workflow: Codex local changes -> ChatGPT audit -> user commit/push",
    ),
    "docs/memory/COMMUNICATION.md": (
        "## Execution authority",
        "must not commit, push, update remote refs",
        "commits/pushes accepted changes and runs farm validation.",
        "changes -> ChatGPT audit -> user commit/push",
    ),
    "docs/memory/templates/CODEX_CONTRACT.md": (
        "**Execution authority:**",
        "the user alone commits/pushes and runs the farm",
        "Workflow: Codex local changes -> ChatGPT audit -> user commit/push",
    ),
}


def default_root() -> Path:
    return Path(__file__).resolve().parents[1]


def markdown_links_text(text: str) -> set[str]:
    links: set[str] = set()
    for target in MARKDOWN_LINK.findall(text):
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


def git_commit_resolver(root: Path, commit: str) -> bool | None:
    try:
        repository = subprocess.run(
            ["git", "rev-parse", "--is-inside-work-tree"],
            cwd=root,
            capture_output=True,
            text=True,
        )
    except OSError:
        return None
    if repository.returncode != 0:
        return None
    result = subprocess.run(
        ["git", "cat-file", "-e", f"{commit}^{{commit}}"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    return result.returncode == 0


def current_sections(text: str) -> tuple[list[str], dict[str, str]]:
    matches = list(HEADING.finditer(text))
    headings = [match.group(1) for match in matches]
    contents: dict[str, str] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        contents[match.group(1)] = text[match.end():end].strip()
    return headings, contents


def check_local_links(root: Path, document: Path, text: str, label: str) -> list[str]:
    errors: list[str] = []
    for target in markdown_links_text(text):
        target_path = (document.parent / target).resolve()
        try:
            target_path.relative_to(root.resolve())
        except ValueError:
            continue
        if not target_path.exists():
            errors.append(f"{label} link target does not exist: {target}")
    return errors


def check_current_structure(root: Path, text: str) -> list[str]:
    errors: list[str] = []
    headings, contents = current_sections(text)
    if headings != list(CURRENT_SECTIONS):
        errors.append("CURRENT.md must contain exactly the ordered M.2 active-state sections")
    for heading in CURRENT_SECTIONS:
        if not contents.get(heading):
            errors.append(f"CURRENT.md section is missing or empty: {heading}")
    errors.extend(check_local_links(root, root / "docs/memory/CURRENT.md", text, "CURRENT.md"))
    for block in re.split(r"\n\s*\n", text):
        if "CLOSED / RUNTIME VALIDATED" in block:
            evidence_links = [target for target in markdown_links_text(block) if target.startswith("evidence/")]
            if not evidence_links:
                errors.append("CURRENT.md runtime-closed claim lacks an evidence/ link")
    return errors


def check_active_states(
    root: Path,
    resolver: Callable[[Path, str], bool | None],
) -> list[str]:
    errors: list[str] = []
    states: dict[str, dict[str, str | int]] = {}
    for relative in ACTIVE_SURFACES:
        path = root / relative
        if not path.is_file():
            continue
        try:
            states[relative] = manifest_tool.parse_active_state(path)
        except manifest_tool.ActiveStateError as error:
            errors.append(f"invalid active-state frontmatter in {relative}: {error}")
    if len(states) == len(ACTIVE_SURFACES):
        canonical = states[ACTIVE_SURFACES[0]]
        for relative, state in states.items():
            if state != canonical:
                errors.append(f"active-state metadata differs from CURRENT.md: {relative}")
        for key, value in canonical.items():
            if key.endswith("_commit"):
                if not isinstance(value, str) or not COMMIT.fullmatch(value):
                    errors.append(f"active-state {key} is not a 40-character lowercase SHA")
                    continue
                resolved = resolver(root, value)
                if resolved is False:
                    errors.append(f"active-state {key} does not resolve to a commit: {value}")
    return errors


def run_checks(
    root: Path,
    soft_limit: int,
    hard_limit: int,
    *,
    check_manifest: bool = True,
    commit_resolver: Callable[[Path, str], bool | None] = git_commit_resolver,
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
        text = readme.read_text(encoding="utf-8")
        links = markdown_links_text(text)
        for target in REQUIRED_README_LINKS:
            if target not in links:
                errors.append(f"README missing required control link: {target}")
        errors.extend(check_local_links(root, readme, text, "README"))

    current = root / "docs" / "memory" / "CURRENT.md"
    if current.is_file():
        errors.extend(check_current_structure(root, current.read_text(encoding="utf-8")))
    errors.extend(check_active_states(root, commit_resolver))

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

    for relative, markers in AUTHORITY_MARKERS.items():
        path = root / relative
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        for marker in markers:
            if marker not in text:
                errors.append(f"authority boundary marker is missing from {relative}: {marker}")

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


def fixture_frontmatter() -> str:
    return """---
memory_schema: 2
active_objective: Fixture objective
current_work_item: Fixture work item
active_status: ACTIVE
next_action: Fixture next action
baseline_commit: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
source_commit: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
bundle_profile_commit: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb
---
"""


def fixture_current() -> str:
    sections = []
    for heading in CURRENT_SECTIONS:
        if heading == "Verified State":
            body = "CLOSED / RUNTIME VALIDATED — [evidence](evidence/accepted.md)"
        elif heading == "Relevant References":
            body = "[evidence](evidence/accepted.md)"
        else:
            body = "fixture"
        sections.append(f"## {heading}\n\n{body}")
    return fixture_frontmatter() + "# Fixture CURRENT\n\n" + "\n\n".join(sections) + "\n"


def write_fixture(root: Path) -> None:
    memory = root / "docs" / "memory"
    (memory / "handoffs").mkdir(parents=True)
    (memory / "roadmap").mkdir()
    (memory / "evidence").mkdir()
    (memory / "chats").mkdir()
    (memory / "templates").mkdir()
    (root / "tools").mkdir()
    readme_links = "\n".join(f"[x]({target})" for target in REQUIRED_README_LINKS)
    (memory / "README.md").write_text(readme_links, encoding="utf-8")
    (memory / "CURRENT.md").write_text(fixture_current(), encoding="utf-8")
    for name in ("MEMORY.md", "MAINTENANCE.md"):
        (memory / name).write_text("ok\n", encoding="utf-8")
    for relative, markers in AUTHORITY_MARKERS.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(markers) + "\n", encoding="utf-8")
    (memory / "handoffs" / "CURRENT_HANDOFF.md").write_text(fixture_frontmatter() + "# Handoff\n", encoding="utf-8")
    (memory / "roadmap" / "CURRENT.md").write_text(fixture_frontmatter() + "# Roadmap\n", encoding="utf-8")
    (memory / "evidence" / "accepted.md").write_text("evidence\n", encoding="utf-8")
    (memory / "chats" / "CHAT_INDEX.md").write_text("ok\n", encoding="utf-8")
    for name in ("MEMORY_UPDATE.md", "RUNTIME_EVIDENCE.md"):
        (memory / "templates" / name).write_text("ok\n", encoding="utf-8")
    (memory / "manifest.json").write_text("{}\n", encoding="utf-8")
    for name in ("check_memory_health.py", "update_memory_manifest.py", "memory_bootstrap.py"):
        (root / "tools" / name).write_text("# fixture\n", encoding="utf-8")


def self_test() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        write_fixture(root)
        errors, warnings = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert not errors and not warnings, (errors, warnings)

        current = root / "docs/memory/CURRENT.md"
        current.write_text(current.read_text(encoding="utf-8") + "\n## Next Action\n\nduplicate\n", encoding="utf-8")
        errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert any("exactly the ordered" in error for error in errors), errors
        current.write_text(fixture_current().replace("evidence/accepted.md", "evidence/missing.md"), encoding="utf-8")
        errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert any("link target does not exist" in error for error in errors), errors
        current.write_text(fixture_current().replace("[evidence](evidence/accepted.md)", "no evidence link"), encoding="utf-8")
        errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert any("runtime-closed claim lacks" in error for error in errors), errors
        current.write_text(fixture_current(), encoding="utf-8")

        handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
        handoff.write_text(handoff.read_text(encoding="utf-8").replace("Fixture work item", "different"), encoding="utf-8")
        errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert any("metadata differs" in error for error in errors), errors
        handoff.write_text(fixture_frontmatter() + "# Handoff\n", encoding="utf-8")

        for field in ("baseline_commit", "source_commit"):
            current.write_text(
                fixture_current().replace(f"{field}: " + ("a" if field == "baseline_commit" else "b") * 40 + "\n", ""),
                encoding="utf-8",
            )
            errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
            assert any("invalid active-state frontmatter" in error for error in errors), errors
        current.write_text(
            fixture_current().replace("source_commit:", "scientific_source_commit:"),
            encoding="utf-8",
        )
        errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert any("invalid active-state frontmatter" in error for error in errors), errors

        for field in ("baseline_commit", "source_commit", "bundle_profile_commit"):
            current.write_text(
                fixture_current().replace(f"{field}: " + ("a" if field == "baseline_commit" else "b") * 40, f"{field}: invalid"),
                encoding="utf-8",
            )
            errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
            assert any("not a 40-character" in error for error in errors), errors
        current.write_text(fixture_current(), encoding="utf-8")

        for field in ("baseline_commit", "source_commit"):
            replacement = f"{field}: {'c' * 40}"
            handoff.write_text(
                fixture_frontmatter().replace(
                    f"{field}: " + ("a" if field == "baseline_commit" else "b") * 40,
                    replacement,
                ),
                encoding="utf-8",
            )
            errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
            assert any("metadata differs" in error for error in errors), errors
        handoff.write_text(fixture_frontmatter() + "# Handoff\n", encoding="utf-8")
        errors, _ = run_checks(
            root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False,
            commit_resolver=lambda _root, _commit: False,
        )
        assert any("does not resolve" in error for error in errors), errors
        errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert not errors, errors

        agents = root / "docs/memory/AGENTS.md"
        agents.write_text("missing authority boundary\n", encoding="utf-8")
        errors, _ = run_checks(root, SOFT_LIMIT_BYTES, HARD_LIMIT_BYTES, check_manifest=False)
        assert any("authority boundary marker" in error for error in errors), errors


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
