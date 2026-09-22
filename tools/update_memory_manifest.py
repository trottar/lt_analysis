#!/usr/bin/env python3
"""Generate and verify the non-authoritative KaonLT memory integrity index."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import tempfile
from typing import Any


SCHEMA_VERSION = 2
MANIFEST_RELATIVE = Path("docs") / "memory" / "manifest.json"
CURRENT_RELATIVE = Path("docs") / "memory" / "CURRENT.md"
ACTIVE_STATE_KEYS = (
    "memory_schema",
    "active_objective",
    "current_work_item",
    "active_status",
    "next_action",
    "baseline_commit",
    "source_commit",
    "bundle_profile_commit",
)


class ActiveStateError(ValueError):
    """Raised when the compact CURRENT.md metadata is structurally invalid."""


def default_root() -> Path:
    return Path(__file__).resolve().parents[1]


def git_output(root: Path, *arguments: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", *arguments],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def parse_active_state(path: Path) -> dict[str, str | int]:
    """Parse the deliberately small, dependency-free YAML-like frontmatter."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise ActiveStateError(f"cannot read {path}: {error}") from error
    if not lines or lines[0] != "---":
        raise ActiveStateError(f"{path} must begin with frontmatter")
    try:
        closing = lines.index("---", 1)
    except ValueError as error:
        raise ActiveStateError(f"{path} frontmatter is not closed") from error

    parsed: dict[str, str] = {}
    for line in lines[1:closing]:
        if not line or line.lstrip().startswith("#") or ":" not in line:
            raise ActiveStateError(f"{path} has invalid frontmatter line: {line!r}")
        key, value = line.split(":", 1)
        key, value = key.strip(), value.strip()
        if not key or not value or key in parsed:
            raise ActiveStateError(f"{path} has missing, empty, or duplicate frontmatter field: {line!r}")
        parsed[key] = value

    schema = parsed.get("memory_schema")
    if schema == "3":
        expected = {"memory_schema"}
        if set(parsed) != expected:
            missing = sorted(expected - set(parsed))
            extra = sorted(set(parsed) - expected)
            details = []
            if missing:
                details.append("missing " + ", ".join(missing))
            if extra:
                details.append("unexpected " + ", ".join(extra))
            raise ActiveStateError(f"{path} schema-3 frontmatter fields differ: {'; '.join(details)}")
        return {"memory_schema": 3}

    if schema is not None and schema != "2":
        raise ActiveStateError(f"{path} memory_schema must be 2 or 3")

    expected = set(ACTIVE_STATE_KEYS)
    if set(parsed) != expected:
        missing = sorted(expected - set(parsed))
        extra = sorted(set(parsed) - expected)
        details = []
        if missing:
            details.append("missing " + ", ".join(missing))
        if extra:
            details.append("unexpected " + ", ".join(extra))
        raise ActiveStateError(f"{path} frontmatter fields differ: {'; '.join(details)}")
    result: dict[str, str | int] = {key: parsed[key] for key in ACTIVE_STATE_KEYS}
    result["memory_schema"] = 2
    return result


def versionable_memory_paths(root: Path) -> list[Path]:
    """Return tracked plus nonignored versionable memory files, excluding manifest."""
    memory_root = root / "docs" / "memory"
    manifest = root / MANIFEST_RELATIVE
    listing = git_output(
        root, "ls-files", "-z", "--cached", "--others", "--exclude-standard", "--", "docs/memory"
    )
    if listing is None:
        candidates = memory_root.rglob("*") if memory_root.exists() else []
        paths = [path for path in candidates if path.is_file()]
    else:
        paths = [root / Path(item) for item in listing.split("\0") if item]
    return sorted(
        {path.resolve() for path in paths if path.is_file() and path.resolve() != manifest.resolve()},
        key=lambda path: path.relative_to(root).as_posix(),
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_manifest(root: Path) -> dict[str, Any]:
    entries = []
    for path in versionable_memory_paths(root):
        entries.append(
            {
                "bytes": path.stat().st_size,
                "path": path.relative_to(root).as_posix(),
                "sha256": sha256(path),
            }
        )
    return {
        "active_state": parse_active_state(root / CURRENT_RELATIVE),
        "files": entries,
        "generated_date_utc": datetime.now(timezone.utc).date().isoformat(),
        "observed_git_head": git_output(root, "rev-parse", "HEAD"),
        "schema_version": SCHEMA_VERSION,
    }


def encoded_manifest(manifest: dict[str, Any]) -> bytes:
    return (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")


def write_manifest(root: Path) -> None:
    manifest_path = root / MANIFEST_RELATIVE
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    payload = encoded_manifest(build_manifest(root))
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=manifest_path.parent, prefix=".manifest-", delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(payload)
    try:
        os.replace(temporary, manifest_path)
    finally:
        if temporary.exists():
            temporary.unlink()


def check_manifest(root: Path) -> list[str]:
    manifest_path = root / MANIFEST_RELATIVE
    if not manifest_path.is_file():
        return [f"missing manifest: {MANIFEST_RELATIVE.as_posix()}"]
    try:
        actual = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        return [f"cannot read manifest: {error}"]
    try:
        expected = build_manifest(root)
    except ActiveStateError as error:
        return [f"CURRENT active state is invalid: {error}"]

    problems: list[str] = []
    if actual.get("schema_version") != SCHEMA_VERSION:
        problems.append("manifest schema_version is unsupported")
    if actual.get("active_state") != expected["active_state"]:
        problems.append("manifest active_state differs from CURRENT.md")
    if actual.get("files") != expected["files"]:
        problems.append("memory file inventory, byte count, or SHA-256 differs from manifest")
    if not isinstance(actual.get("generated_date_utc"), str):
        problems.append("manifest generated_date_utc is missing or invalid")
    if actual.get("observed_git_head") is not None and not isinstance(actual["observed_git_head"], str):
        problems.append("manifest observed_git_head is invalid")
    return problems


def fixture_current() -> str:
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
# Current
"""


def self_test() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        memory = root / "docs" / "memory"
        memory.mkdir(parents=True)
        current = memory / "CURRENT.md"
        current.write_text(fixture_current(), encoding="utf-8")
        write_manifest(root)
        assert not check_manifest(root), "fresh manifest did not verify"
        current.write_text(fixture_current().replace("Fixture objective", "Changed objective"), encoding="utf-8")
        assert check_manifest(root), "changed bytes or active state were not detected"
        current.write_text("# No frontmatter\n", encoding="utf-8")
        try:
            parse_active_state(current)
        except ActiveStateError:
            pass
        else:
            raise AssertionError("missing frontmatter was accepted")
        for field in ("baseline_commit", "source_commit"):
            current.write_text(
                fixture_current().replace(f"{field}: " + ("a" if field == "baseline_commit" else "b") * 40 + "\n", ""),
                encoding="utf-8",
            )
            try:
                parse_active_state(current)
            except ActiveStateError:
                pass
            else:
                raise AssertionError(f"missing {field} was accepted")
        current.write_text(
            fixture_current().replace("source_commit:", "scientific_source_commit:"),
            encoding="utf-8",
        )
        try:
            parse_active_state(current)
        except ActiveStateError:
            pass
        else:
            raise AssertionError("obsolete scientific_source_commit was accepted")
        current.write_text("---\nmemory_schema: 3\n---\n# Current\n", encoding="utf-8")
        assert parse_active_state(current) == {"memory_schema": 3}
        write_manifest(root)
        assert not check_manifest(root), "schema-3 manifest did not verify"
        current.write_text("---\nmemory_schema: 3\nactive_status: forbidden\n---\n", encoding="utf-8")
        try:
            parse_active_state(current)
        except ActiveStateError:
            pass
        else:
            raise AssertionError("schema-3 extra frontmatter was accepted")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true", help="generate docs/memory/manifest.json")
    mode.add_argument("--check", action="store_true", help="verify current memory bytes against manifest")
    mode.add_argument("--self-test", action="store_true", help="run deterministic temporary-fixture checks")
    parser.add_argument("--root", type=Path, default=default_root(), help="repository root")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        self_test()
        print("SELF-TEST: PASS")
        return 0
    root = args.root.resolve()
    try:
        if args.write:
            write_manifest(root)
            print(f"MANIFEST: WROTE {MANIFEST_RELATIVE.as_posix()}")
            return 0
        problems = check_manifest(root)
    except ActiveStateError as error:
        print(f"MANIFEST: FAIL: CURRENT active state is invalid: {error}")
        return 1
    if problems:
        for problem in problems:
            print(f"MANIFEST: FAIL: {problem}")
        return 1
    print("MANIFEST: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
