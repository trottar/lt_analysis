#!/usr/bin/env python3
"""Check deterministic schema-3 structural and semantic KaonLT memory health."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from typing import Iterable


REQUIRED_FILES = (
    "docs/memory/README.md", "docs/memory/AGENTS.md", "docs/memory/CURRENT.md",
    "docs/memory/MEMORY.md", "docs/memory/handoffs/CURRENT_HANDOFF.md",
    "docs/memory/USER.md", "docs/memory/TOOLS.md", "docs/memory/COMMUNICATION.md",
    "docs/memory/CODEX.md", "docs/memory/MAINTENANCE.md", "docs/memory/LEARNINGS.md",
    "docs/memory/roadmap/STATUS.md", "docs/memory/sources/SOURCE_INDEX.md",
    "docs/memory/sources/ARTIFACT_INDEX.md", "docs/memory/chats/CHAT_INDEX.md",
    "docs/memory/history/CHAT_INDEX.md", "docs/memory/history/PROJECT_HISTORY.md",
    "docs/memory/import/README.md", "docs/memory/templates/CODEX_CONTRACT.md",
    "docs/memory/templates/MEMORY_UPDATE.md", "docs/memory/templates/RUNTIME_EVIDENCE.md",
    "docs/memory/manifest.json", "tools/check_memory_health.py",
    "tools/update_memory_manifest.py", "tools/memory_bootstrap.py",
)
LIVE_SURFACES = {
    "docs/memory/README.md": "KaonLT development memory",
    "docs/memory/AGENTS.md": "KaonLT memory operating rules",
    "docs/memory/CURRENT.md": "Current KaonLT development state",
    "docs/memory/MEMORY.md": "Durable KaonLT project knowledge",
    "docs/memory/handoffs/CURRENT_HANDOFF.md": "Current KaonLT handoff",
    "docs/memory/USER.md": "KaonLT user collaboration context",
    "docs/memory/TOOLS.md": "KaonLT tools and operational commands",
    "docs/memory/COMMUNICATION.md": "KaonLT farm communication",
    "docs/memory/CODEX.md": "KaonLT Codex workflow",
    "docs/memory/MAINTENANCE.md": "KaonLT memory maintenance",
    "docs/memory/LEARNINGS.md": "KaonLT durable learnings",
    "docs/memory/roadmap/STATUS.md": "Approved KaonLT roadmap status",
    "docs/memory/sources/SOURCE_INDEX.md": "KaonLT source and evidence index",
    "docs/memory/sources/ARTIFACT_INDEX.md": "KaonLT artifact index",
    "docs/memory/chats/CHAT_INDEX.md": "KaonLT chat-era index",
    "docs/memory/history/CHAT_INDEX.md": "KaonLT chat and handoff index",
    "docs/memory/history/PROJECT_HISTORY.md": "KaonLT project history",
    "docs/memory/import/README.md": "KaonLT imported historical provenance",
}
NO_NEXT_SURFACES = tuple(relative for relative in LIVE_SURFACES if relative != "docs/memory/CURRENT.md")
SIZE_LIMITS = {
    "docs/memory/CURRENT.md": (8 * 1024, 16 * 1024),
    "docs/memory/handoffs/CURRENT_HANDOFF.md": (6 * 1024, 10 * 1024),
    "docs/memory/MEMORY.md": (30 * 1024, 50 * 1024),
}
CURRENT_SECTIONS = (
    "Active Objective", "Current Work Item", "Verified State", "Source / Evidence Identity",
    "Blockers", "Next Action", "Success Criteria", "Do Not Reopen Without New Evidence",
    "Relevant References",
)
MEMORY_SECTIONS = (
    "Authority and evidence", "Source and provenance semantics",
    "Scientific ownership and production ordering", "HGCer and Method-A/Method-B boundaries",
    "Diagnostics and runtime validation", "Canonical record ownership",
)
HANDOFF_SECTIONS = ("Transfer State", "Resume")
BOOTSTRAP_ORDER = ("AGENTS.md", "CURRENT.md", "MEMORY.md", "handoffs/CURRENT_HANDOFF.md", "USER.md")
REQUIRED_README_LINKS = (
    "AGENTS.md", "CURRENT.md", "MEMORY.md", "handoffs/CURRENT_HANDOFF.md", "USER.md",
    "TOOLS.md", "COMMUNICATION.md", "CODEX.md", "MAINTENANCE.md", "LEARNINGS.md",
    "chats/CHAT_INDEX.md", "sources/SOURCE_INDEX.md", "sources/ARTIFACT_INDEX.md",
    "history/PROJECT_HISTORY.md", "import/README.md", "roadmap/STATUS.md",
    "templates/CODEX_CONTRACT.md", "templates/MEMORY_UPDATE.md",
    "templates/RUNTIME_EVIDENCE.md", "manifest.json",
)
TOP_LEVEL_ROLES = (
    "README", "CURRENT", "MEMORY", "AGENTS", "MAINTENANCE", "COMMUNICATION", "CODEX",
    "USER", "TOOLS", "LEARNINGS",
)
APPROVED_STATUS_LABELS = (
    "CLOSED / RUNTIME VALIDATED", "SOURCE REVIEWED",
    "DEVELOPMENT COMPLETE, FARM VALIDATION PENDING", "ACTIVE", "DEFERRED", "BLOCKED", "NEXT",
)
RUNTIME_CLOSED = "CLOSED / RUNTIME VALIDATED"
MARKDOWN_LINK = re.compile(r"(?<!!)\[[^\]]+\]\(([^)\s]+)(?:\s+[^)]*)?\)")
H1 = re.compile(r"^#\s+(.+?)\s*$", re.MULTILINE)
H2 = re.compile(r"^##\s+(.+?)\s*$", re.MULTILINE)
PHASE_HEADING = re.compile(r"^#{3,4}\s+.+?$", re.MULTILINE)
FENCE = re.compile(r"^\s*```", re.MULTILINE)
ROLE_SUFFIX = re.compile(r"(?:[._-](?:copy|backup|old|v\d+))+$", re.IGNORECASE)
LOWER_SHA40 = re.compile(r"(?<![0-9a-f])[0-9a-f]{40}(?![0-9a-f])")
LOWER_HASH64 = re.compile(r"(?<![0-9a-f])[0-9a-f]{64}(?![0-9a-f])")
STATUS_CODE = re.compile(r"`([^`\n]+)`")
AUTHORITY_MARKERS = {
    "docs/memory/AGENTS.md": (
        "## Execution authority", "Codex must not commit, push, update remote refs",
        "The user alone commits/pushes accepted changes and runs farm",
        "Workflow: Codex local changes -> ChatGPT audit -> user commit/push",
    ),
    "docs/memory/CODEX.md": (
        "Codex must not commit, push, update remote refs",
        "The user alone commits/pushes accepted changes and runs farm",
        "Workflow: Codex local changes -> ChatGPT audit -> user commit/push",
    ),
    "docs/memory/COMMUNICATION.md": (
        "## Execution authority", "must not commit, push, update remote refs",
        "commits/pushes accepted changes and runs farm validation.",
        "changes -> ChatGPT audit -> user commit/push",
    ),
    "docs/memory/templates/CODEX_CONTRACT.md": (
        "**Execution authority:**", "the user alone commits/pushes and runs the farm",
        "Workflow: Codex local changes -> ChatGPT audit -> user commit/push",
    ),
}
MAINTENANCE_THRESHOLD_ROWS = (
    "| `docs/memory/CURRENT.md` | 8 KiB | 16 KiB |",
    "| `docs/memory/handoffs/CURRENT_HANDOFF.md` | 6 KiB | 10 KiB |",
    "| `docs/memory/MEMORY.md` | 30 KiB | 50 KiB |",
)


class CurrentFrontmatterError(ValueError):
    """Raised when CURRENT.md does not use the minimal schema-3 frontmatter."""


def default_root() -> Path:
    return Path(__file__).resolve().parents[1]


def parse_current_schema3(path: Path) -> dict[str, int]:
    """Parse the health-owned, deliberately minimal CURRENT frontmatter."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise CurrentFrontmatterError(f"cannot read {path}: {error}") from error
    if not lines or lines[0] != "---":
        raise CurrentFrontmatterError(f"{path} must begin with frontmatter")
    try:
        closing = lines.index("---", 1)
    except ValueError as error:
        raise CurrentFrontmatterError(f"{path} frontmatter is not closed") from error
    parsed: dict[str, str] = {}
    for line in lines[1:closing]:
        if not line or line.lstrip().startswith("#") or ":" not in line:
            raise CurrentFrontmatterError(f"{path} has invalid frontmatter line: {line!r}")
        key, value = (part.strip() for part in line.split(":", 1))
        if not key or not value or key in parsed:
            raise CurrentFrontmatterError(f"{path} has missing, empty, or duplicate frontmatter field: {line!r}")
        parsed[key] = value
    if set(parsed) != {"memory_schema"}:
        missing = sorted({"memory_schema"} - set(parsed))
        extra = sorted(set(parsed) - {"memory_schema"})
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unexpected " + ", ".join(extra))
        raise CurrentFrontmatterError(f"{path} schema-3 frontmatter fields differ: {'; '.join(details)}")
    if parsed["memory_schema"] != "3":
        raise CurrentFrontmatterError(f"{path} memory_schema must be 3")
    return {"memory_schema": 3}


def markdown_links_text(text: str) -> set[str]:
    links: set[str] = set()
    for target in MARKDOWN_LINK.findall(text):
        target = target.strip().split("#", 1)[0]
        if target and "://" not in target and not target.startswith("mailto:"):
            links.add(target)
    return links


def role_base(stem: str) -> str:
    return ROLE_SUFFIX.sub("", stem).upper()


def read_document(path: Path, label: str) -> tuple[str | None, list[str]]:
    try:
        data = path.read_bytes()
    except OSError as error:
        return None, [f"cannot read {label}: {error}"]
    if b"\0" in data:
        return None, [f"{label} contains a NUL byte"]
    try:
        text = data.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    except UnicodeDecodeError as error:
        return None, [f"{label} is not readable UTF-8: {error}"]
    if not data.endswith(b"\n"):
        return text, [f"{label} must end with a newline"]
    return text, []


def current_sections(text: str) -> tuple[list[str], dict[str, str]]:
    matches = list(H2.finditer(text))
    headings = [match.group(1) for match in matches]
    contents: dict[str, str] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        contents[match.group(1)] = text[match.end():end].strip()
    return headings, contents


def resolve_local_link(document: Path, target: str) -> Path | None:
    if "://" in target or target.startswith("mailto:"):
        return None
    target = target.split("#", 1)[0]
    return (document.parent / target).resolve() if target else None


def check_local_links(root: Path, document: Path, text: str, label: str) -> list[str]:
    errors: list[str] = []
    for target in markdown_links_text(text):
        resolved = resolve_local_link(document, target)
        if resolved is None:
            continue
        try:
            resolved.relative_to(root.resolve())
        except ValueError:
            errors.append(f"{label} link target escapes repository: {target}")
        else:
            if not resolved.exists():
                errors.append(f"{label} link target does not exist: {target}")
    return errors


def check_representation(root: Path, relative: str, title: str) -> tuple[str | None, list[str]]:
    text, errors = read_document(root / relative, relative)
    if text is None:
        return None, errors
    if H1.findall(text) != [title]:
        errors.append(f"{relative} must contain exactly one H1: # {title}")
    if len(FENCE.findall(text)) % 2:
        errors.append(f"{relative} has an unbalanced fenced triple-backtick block")
    errors.extend(check_local_links(root, root / relative, text, relative))
    return text, errors


def check_size_limit(relative: str, size: int, soft_limit: int, hard_limit: int) -> tuple[list[str], list[str]]:
    if size > hard_limit:
        return [f"{relative} exceeds hard limit ({size} > {hard_limit})"], []
    if size > soft_limit:
        return [], [f"{relative} exceeds soft limit ({size} > {soft_limit})"]
    return [], []


def check_size_limits(root: Path) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    for relative, (soft_limit, hard_limit) in SIZE_LIMITS.items():
        path = root / relative
        if path.is_file():
            role_errors, role_warnings = check_size_limit(relative, path.stat().st_size, soft_limit, hard_limit)
            errors.extend(role_errors)
            warnings.extend(role_warnings)
    return errors, warnings


def check_startup_contract(text: str, label: str) -> list[str]:
    start = text.lower().find("read these five files")
    startup = text[start:] if start >= 0 else text
    positions = [startup.find(relative) for relative in BOOTSTRAP_ORDER]
    errors = []
    if any(position < 0 for position in positions) or positions != sorted(positions):
        errors.append(f"{label} five-file startup order differs from the required order")
    normalized_startup = re.sub(r"\s+", " ", startup)
    for marker in ("CURRENT direct references", "exact active task", "required canonical evidence/decision/phase records"):
        if marker not in normalized_startup:
            errors.append(f"{label} selective-expansion sources are incomplete")
            break
    return errors


def check_bootstrap(texts: dict[str, str]) -> list[str]:
    errors: list[str] = []
    agents = texts.get("docs/memory/AGENTS.md", "")
    if "Read these five files in full, in this exact order:" not in agents:
        errors.append("AGENTS.md must require the full five-file core read")
    if "Only after those five are read in full may a session expand selectively from:" not in agents:
        errors.append("AGENTS.md must require full-read-before-selective-expansion")
    if "Do not eagerly load the whole memory hierarchy." not in agents:
        errors.append("AGENTS.md must prohibit eager hierarchy expansion")
    errors.extend(check_startup_contract(agents, "AGENTS.md"))

    readme = texts.get("docs/memory/README.md", "")
    if "read these five files in full and in this exact order" not in readme:
        errors.append("README.md must require the full five-file core read")
    if "Only then expand selectively" not in readme or "Do not load" not in readme or "whole hierarchy" not in readme:
        errors.append("README.md must preserve full-read-before-selective-expansion meaning")
    errors.extend(check_startup_contract(readme, "README.md"))

    maintenance = texts.get("docs/memory/MAINTENANCE.md", "")
    if "read these five files in full, in this exact order" not in maintenance:
        errors.append("MAINTENANCE.md must require the full five-file core read")
    if "Only after the five-file core is read" not in maintenance or "Do not eagerly load" not in maintenance or "whole" not in maintenance or "hierarchy" not in maintenance:
        errors.append("MAINTENANCE.md must preserve full-read-before-selective-expansion meaning")
    errors.extend(check_startup_contract(maintenance, "MAINTENANCE.md"))
    return errors


def check_runtime_evidence_blocks(root: Path, document: Path, text: str, label: str) -> list[str]:
    errors: list[str] = []
    evidence_root = (root / "docs/memory/evidence").resolve()
    for block in re.split(r"\n\s*\n", text):
        if RUNTIME_CLOSED not in block:
            continue
        valid = False
        for target in markdown_links_text(block):
            resolved = resolve_local_link(document, target)
            if resolved is None:
                continue
            try:
                resolved.relative_to(evidence_root)
            except ValueError:
                continue
            if resolved.is_file():
                valid = True
                break
        if not valid:
            errors.append(f"{label} runtime-closed claim lacks a direct canonical evidence link")
    return errors


def check_current(root: Path, text: str) -> list[str]:
    errors: list[str] = []
    try:
        parsed = parse_current_schema3(root / "docs/memory/CURRENT.md")
    except CurrentFrontmatterError as error:
        errors.append(f"invalid schema-3 CURRENT frontmatter: {error}")
    else:
        if parsed != {"memory_schema": 3}:
            errors.append("CURRENT.md must use only minimal schema-3 frontmatter")
    headings, contents = current_sections(text)
    if headings != list(CURRENT_SECTIONS):
        errors.append("CURRENT.md must contain exactly the ordered active-state sections")
    for heading in CURRENT_SECTIONS:
        if not contents.get(heading):
            errors.append(f"CURRENT.md section is missing or empty: {heading}")
    next_count = text.count("NEXT —")
    if next_count != 1:
        errors.append("CURRENT.md must contain exactly one active NEXT — statement")
    elif contents.get("Next Action", "").count("NEXT —") != 1:
        errors.append("CURRENT.md NEXT — statement must occur only in ## Next Action")
    errors.extend(check_runtime_evidence_blocks(root, root / "docs/memory/CURRENT.md", text, "CURRENT.md"))
    for label in STATUS_CODE.findall(text):
        if re.fullmatch(r"[A-Z][A-Z ,/]+", label) and label not in APPROVED_STATUS_LABELS:
            errors.append(f"CURRENT.md has unapproved uppercase status label: {label}")
    return errors


def check_memory(text: str) -> list[str]:
    errors: list[str] = []
    headings, _ = current_sections(text)
    if headings != list(MEMORY_SECTIONS):
        errors.append("MEMORY.md must contain exactly the ordered durable-knowledge sections")
    if text.startswith("---\n"):
        errors.append("MEMORY.md must not contain YAML/frontmatter")
    if "NEXT —" in text:
        errors.append("MEMORY.md must not contain active NEXT —")
    if (set(CURRENT_SECTIONS) | set(HANDOFF_SECTIONS)).intersection(headings):
        errors.append("MEMORY.md must not contain active-state or handoff headings")
    if LOWER_SHA40.search(text):
        errors.append("MEMORY.md must not contain a full 40-character lowercase SHA")
    if LOWER_HASH64.search(text):
        errors.append("MEMORY.md must not contain a full 64-character lowercase hash")
    for status in (RUNTIME_CLOSED, "DEVELOPMENT COMPLETE, FARM VALIDATION PENDING", "SOURCE REVIEWED"):
        if status in text:
            errors.append(f"MEMORY.md must not contain phase-ledger status label: {status}")
    return errors


def check_handoff(text: str) -> list[str]:
    errors: list[str] = []
    headings, contents = current_sections(text)
    if text.startswith("---\n"):
        errors.append("handoff must not contain active-state frontmatter")
    if headings != list(HANDOFF_SECTIONS):
        errors.append("handoff must contain exactly the Transfer State and Resume sections")
    transfer = contents.get("Transfer State", "")
    if not transfer:
        errors.append("handoff Transfer State must not be empty")
    if "No exceptional transfer state" in transfer and transfer.strip() != "No exceptional transfer state is recorded.":
        errors.append("handoff stable Transfer State must use the exact no-transfer sentence")
    normalized = text.replace("`", "")
    if "CURRENT.md is the sole authoritative resumable state." not in normalized:
        errors.append("handoff must state CURRENT.md sole authority")
    if "handoff cannot override CURRENT.md" not in normalized:
        errors.append("handoff must state it cannot override CURRENT.md")
    if "CURRENT.md" not in contents.get("Resume", ""):
        errors.append("handoff Resume must point to CURRENT.md")
    if "NEXT —" in text:
        errors.append("handoff must not contain active NEXT —")
    return errors


def phase_sections(text: str) -> list[str]:
    matches = list(PHASE_HEADING.finditer(text))
    sections: list[str] = []
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        sections.append(text[match.end():end].strip())
    return sections


def check_roadmap(root: Path, text: str) -> list[str]:
    errors: list[str] = []
    if (root / "docs/memory/roadmap/CURRENT.md").exists():
        errors.append("legacy roadmap/CURRENT.md must not exist")
    if text.startswith("---\n"):
        errors.append("roadmap STATUS must not contain active-state frontmatter")
    if re.search(r"^## NEXT\s*$", text, re.MULTILINE):
        errors.append("roadmap STATUS must not contain ## NEXT")
    if "NEXT —" in text:
        errors.append("roadmap STATUS must not contain active NEXT —")
    normalized = text.replace("`", "")
    if any(marker not in normalized for marker in ("CURRENT.md", "active objective", "blockers", "next action")):
        errors.append("roadmap STATUS must point active objective, blockers, and next action ownership to CURRENT.md")
    for section in phase_sections(text):
        first_line = section.splitlines()[0] if section else ""
        match = re.match(r"^`([^`]+)`\s+—", first_line)
        if not match or match.group(1) not in APPROVED_STATUS_LABELS:
            errors.append("roadmap phase-state claim must begin with an exact backticked approved label")
            break
    errors.extend(check_runtime_evidence_blocks(root, root / "docs/memory/roadmap/STATUS.md", text, "roadmap STATUS"))
    return errors


def check_roles(root: Path) -> list[str]:
    errors: list[str] = []
    memory = root / "docs/memory"
    roles: dict[str, list[str]] = {role: [] for role in TOP_LEVEL_ROLES}
    if memory.is_dir():
        for path in memory.iterdir():
            if path.is_file() and path.suffix.lower() == ".md" and role_base(path.stem) in roles:
                roles[role_base(path.stem)].append(path.name)
    for role, names in roles.items():
        if names != [f"{role}.md"]:
            errors.append(f"duplicate or missing top-level {role} role: {', '.join(sorted(names)) or 'none'}")
    for directory, canonical, label in (
        (memory / "handoffs", "CURRENT_HANDOFF", "handoff CURRENT_HANDOFF"),
        (memory / "roadmap", "STATUS", "roadmap STATUS"),
    ):
        names = sorted(
            path.name for path in directory.glob("*.md") if role_base(path.stem) == canonical
        ) if directory.is_dir() else []
        if names != [f"{canonical}.md"]:
            errors.append(f"duplicate or missing {label} role: {', '.join(names) or 'none'}")
    return errors


def check_readme_navigation(text: str) -> list[str]:
    links = markdown_links_text(text)
    return [f"README missing required control link: {target}" for target in REQUIRED_README_LINKS if target not in links]


def check_maintenance_thresholds(text: str) -> list[str]:
    return [f"MAINTENANCE.md is missing final role-specific threshold row: {row}" for row in MAINTENANCE_THRESHOLD_ROWS if row not in text]


def run_manifest_check(root: Path) -> tuple[bool, str]:
    script = root / "tools/update_memory_manifest.py"
    if not script.is_file():
        return False, "manifest checker is missing"
    result = subprocess.run([sys.executable, str(script), "--root", str(root), "--check"], cwd=root, capture_output=True, text=True)
    return result.returncode == 0, (result.stdout + result.stderr).strip()


def run_checks(root: Path, *, check_manifest: bool = True) -> tuple[list[str], list[str]]:
    """Return deterministic schema-3 health errors and size-policy warnings."""
    errors: list[str] = []
    warnings: list[str] = []
    for relative in REQUIRED_FILES:
        if not (root / relative).is_file():
            errors.append(f"missing required file: {relative}")
    texts: dict[str, str] = {}
    for relative, title in LIVE_SURFACES.items():
        if not (root / relative).is_file():
            continue
        text, representation_errors = check_representation(root, relative, title)
        errors.extend(representation_errors)
        if text is not None:
            texts[relative] = text
            if relative != "docs/memory/CURRENT.md" and text.startswith("---\n"):
                errors.append(f"only CURRENT.md may contain active-state frontmatter: {relative}")
    size_errors, size_warnings = check_size_limits(root)
    errors.extend(size_errors)
    warnings.extend(size_warnings)
    if "docs/memory/CURRENT.md" in texts:
        errors.extend(check_current(root, texts["docs/memory/CURRENT.md"]))
    if "docs/memory/MEMORY.md" in texts:
        errors.extend(check_memory(texts["docs/memory/MEMORY.md"]))
    if "docs/memory/handoffs/CURRENT_HANDOFF.md" in texts:
        errors.extend(check_handoff(texts["docs/memory/handoffs/CURRENT_HANDOFF.md"]))
    if "docs/memory/roadmap/STATUS.md" in texts:
        errors.extend(check_roadmap(root, texts["docs/memory/roadmap/STATUS.md"]))
    errors.extend(check_bootstrap(texts))
    for relative in NO_NEXT_SURFACES:
        if "NEXT —" in texts.get(relative, ""):
            errors.append(f"non-CURRENT live/control surface must not contain active NEXT —: {relative}")
    if "docs/memory/README.md" in texts:
        errors.extend(check_readme_navigation(texts["docs/memory/README.md"]))
    if "docs/memory/MAINTENANCE.md" in texts:
        errors.extend(check_maintenance_thresholds(texts["docs/memory/MAINTENANCE.md"]))
    errors.extend(check_roles(root))
    for relative, markers in AUTHORITY_MARKERS.items():
        text = texts.get(relative)
        if text is None:
            text, read_errors = read_document(root / relative, relative)
            errors.extend(read_errors)
        if text is not None:
            for marker in markers:
                if marker not in text:
                    errors.append(f"authority boundary marker is missing from {relative}: {marker}")
    if check_manifest and not errors:
        passed, message = run_manifest_check(root)
        if not passed:
            errors.append(f"manifest freshness check failed: {message or 'no diagnostic'}")
    return errors, warnings


def print_result(errors: Iterable[str], warnings: Iterable[str]) -> int:
    errors, warnings = list(errors), list(warnings)
    for warning in warnings:
        print(f"MEMORY HEALTH: WARN: {warning}")
    for error in errors:
        print(f"MEMORY HEALTH: FAIL: {error}")
    if errors:
        return 1
    print("MEMORY HEALTH: WARN" if warnings else "MEMORY HEALTH: PASS")
    return 0


def fixture_current() -> str:
    sections = []
    for heading in CURRENT_SECTIONS:
        if heading == "Verified State":
            body = "`CLOSED / RUNTIME VALIDATED` — [accepted evidence](evidence/accepted.md)"
        elif heading == "Next Action":
            body = "NEXT — Complete the one fixture action."
        elif heading == "Relevant References":
            body = "[accepted evidence](evidence/accepted.md)"
        else:
            body = "fixture"
        sections.append(f"## {heading}\n\n{body}")
    return "---\nmemory_schema: 3\n---\n# Current KaonLT development state\n\n" + "\n\n".join(sections) + "\n"


def write_fixture(root: Path) -> None:
    """Create a complete, valid, schema-3-only fixture tree."""
    memory = root / "docs/memory"
    for directory in (memory / "handoffs", memory / "roadmap", memory / "evidence", memory / "sources", memory / "chats", memory / "history", memory / "import", memory / "templates", root / "tools"):
        directory.mkdir(parents=True, exist_ok=True)
    startup = (
        "Read these five files in full, in this exact order:\n\n1. `AGENTS.md`\n2. `CURRENT.md`\n3. `MEMORY.md`\n4. `handoffs/CURRENT_HANDOFF.md`\n5. `USER.md`\n\n"
        "Only after those five are read in full may a session expand selectively from:\n\n- CURRENT direct references\n- exact active task\n- required canonical evidence/decision/phase records\n\nDo not eagerly load the whole memory hierarchy.\n"
    )
    readme_startup = (
        "read these five files in full and in this exact order:\n\n1. [AGENTS.md](AGENTS.md)\n2. [CURRENT.md](CURRENT.md)\n3. [MEMORY.md](MEMORY.md)\n4. [handoffs/CURRENT_HANDOFF.md](handoffs/CURRENT_HANDOFF.md)\n5. [USER.md](USER.md)\n\n"
        "Only then expand selectively from CURRENT direct references, the exact active task, and required canonical evidence/decision/phase records. Do not load the whole hierarchy without task-specific need.\n"
    )
    maintenance_startup = (
        "read these five files in full, in this exact order:\n\n1. `AGENTS.md`\n2. `CURRENT.md`\n3. `MEMORY.md`\n4. `handoffs/CURRENT_HANDOFF.md`\n5. `USER.md`\n\n"
        "Only after the five-file core is read may task-directed expansion use CURRENT direct references, the exact active task, and required canonical evidence/decision/phase records. Do not eagerly load the whole memory hierarchy.\n"
    )
    navigation = "\n".join(f"[x]({target})" for target in REQUIRED_README_LINKS)
    documents = {
        "README.md": "# KaonLT development memory\n\n" + readme_startup + "\n" + navigation + "\n",
        "AGENTS.md": "# KaonLT memory operating rules\n\n" + startup + "\n## Execution authority\n\nCodex must not commit, push, update remote refs. The user alone commits/pushes accepted changes and runs farm.\nWorkflow: Codex local changes -> ChatGPT audit -> user commit/push\n",
        "CURRENT.md": fixture_current(),
        "MEMORY.md": "# Durable KaonLT project knowledge\n\n" + "\n\n".join(f"## {heading}\n\nfixture durable knowledge" for heading in MEMORY_SECTIONS) + "\n",
        "USER.md": "# KaonLT user collaboration context\n\nfixture\n",
        "TOOLS.md": "# KaonLT tools and operational commands\n\nfixture\n",
        "COMMUNICATION.md": "# KaonLT farm communication\n\n## Execution authority\n\nCodex must not commit, push, update remote refs; the user commits/pushes accepted changes and runs farm validation.\nCodex local changes -> ChatGPT audit -> user commit/push\n",
        "CODEX.md": "# KaonLT Codex workflow\n\nCodex must not commit, push, update remote refs. The user alone commits/pushes accepted changes and runs farm.\nWorkflow: Codex local changes -> ChatGPT audit -> user commit/push\n",
        "MAINTENANCE.md": "# KaonLT memory maintenance\n\n" + maintenance_startup + "\n" + "\n".join(MAINTENANCE_THRESHOLD_ROWS) + "\n",
        "LEARNINGS.md": "# KaonLT durable learnings\n\nfixture\n",
        "handoffs/CURRENT_HANDOFF.md": "# Current KaonLT handoff\n\n## Transfer State\n\nNo exceptional transfer state is recorded.\n\n## Resume\n\nCURRENT.md is the sole authoritative resumable state. The handoff cannot override CURRENT.md.\n",
        "roadmap/STATUS.md": "# Approved KaonLT roadmap status\n\n### Fixture phase\n\n`ACTIVE` — fixture dependency/status structure.\n\nCURRENT.md owns the active objective, blockers, and next action.\n",
        "sources/SOURCE_INDEX.md": "# KaonLT source and evidence index\n\nfixture\n",
        "sources/ARTIFACT_INDEX.md": "# KaonLT artifact index\n\nfixture\n",
        "chats/CHAT_INDEX.md": "# KaonLT chat-era index\n\nfixture\n",
        "history/CHAT_INDEX.md": "# KaonLT chat and handoff index\n\nfixture\n",
        "history/PROJECT_HISTORY.md": "# KaonLT project history\n\nfixture\n",
        "import/README.md": "# KaonLT imported historical provenance\n\nfixture\n",
        "templates/CODEX_CONTRACT.md": "# Fixture contract\n\n**Execution authority:** the user alone commits/pushes and runs the farm.\nWorkflow: Codex local changes -> ChatGPT audit -> user commit/push\n",
        "templates/MEMORY_UPDATE.md": "# Fixture memory update\n\nfixture\n",
        "templates/RUNTIME_EVIDENCE.md": "# Fixture runtime evidence\n\nfixture\n",
    }
    for relative, text in documents.items():
        (memory / relative).write_text(text, encoding="utf-8")
    (memory / "evidence/accepted.md").write_text("# Evidence\n\nfixture\n", encoding="utf-8")
    (memory / "manifest.json").write_text("{}\n", encoding="utf-8")
    for name in ("check_memory_health.py", "update_memory_manifest.py", "memory_bootstrap.py"):
        (root / "tools" / name).write_text("# fixture\n", encoding="utf-8")


def self_test() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        write_fixture(root)
        errors, warnings = run_checks(root, check_manifest=False)
        assert not errors and not warnings, (errors, warnings)
        current = root / "docs/memory/CURRENT.md"
        current.write_text(current.read_text(encoding="utf-8") + "\n## Next Action\n\nduplicate\n", encoding="utf-8")
        errors, _ = run_checks(root, check_manifest=False)
        assert any("exactly the ordered" in error for error in errors), errors
        current.write_text(fixture_current().replace("[accepted evidence](evidence/accepted.md)", "missing evidence"), encoding="utf-8")
        errors, _ = run_checks(root, check_manifest=False)
        assert any("runtime-closed claim lacks" in error for error in errors), errors
        current.write_text(fixture_current().replace("evidence/accepted.md", "evidence/missing.md"), encoding="utf-8")
        errors, _ = run_checks(root, check_manifest=False)
        assert any("link target does not exist" in error for error in errors), errors
        write_fixture(root)
        handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
        handoff.write_text(handoff.read_text(encoding="utf-8").replace("sole authoritative", "ordinary"), encoding="utf-8")
        errors, _ = run_checks(root, check_manifest=False)
        assert any("sole authority" in error for error in errors), errors
        write_fixture(root)
        memory = root / "docs/memory/MEMORY.md"
        memory.write_text(memory.read_text(encoding="utf-8") + "\n## Next Action\n\nfixture\n", encoding="utf-8")
        errors, _ = run_checks(root, check_manifest=False)
        assert any("durable-knowledge" in error for error in errors), errors
        write_fixture(root)
        legacy = root / "docs/memory/roadmap/CURRENT.md"
        legacy.write_text("# Legacy\n", encoding="utf-8")
        errors, _ = run_checks(root, check_manifest=False)
        assert any("legacy roadmap" in error for error in errors), errors
        legacy.unlink()
        roadmap = root / "docs/memory/roadmap/STATUS.md"
        roadmap.write_text(roadmap.read_text(encoding="utf-8") + "\nNEXT — fixture\n", encoding="utf-8")
        errors, _ = run_checks(root, check_manifest=False)
        assert any("must not contain active NEXT" in error for error in errors), errors
        write_fixture(root)
        duplicate = root / "docs/memory/CURRENT_copy.md"
        duplicate.write_text("# Copy\n", encoding="utf-8")
        errors, _ = run_checks(root, check_manifest=False)
        assert any("duplicate or missing top-level CURRENT role" in error for error in errors), errors
        duplicate.unlink()
        current = root / "docs/memory/CURRENT.md"
        current.write_bytes(current.read_bytes() + b"x" * (SIZE_LIMITS["docs/memory/CURRENT.md"][1] + 1))
        errors, _ = run_checks(root, check_manifest=False)
        assert any("CURRENT.md exceeds hard limit" in error for error in errors), errors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root(), help="repository root")
    parser.add_argument("--self-test", action="store_true", help="run deterministic temporary-fixture checks")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        self_test()
        print("SELF-TEST: PASS")
        return 0
    errors, warnings = run_checks(args.root.resolve())
    return print_result(errors, warnings)


if __name__ == "__main__":
    raise SystemExit(main())
