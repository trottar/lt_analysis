#!/usr/bin/env python3
"""Generate and verify the KaonLT memory file-integrity index."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
from typing import Any


SCHEMA_VERSION = 3
MANIFEST_RELATIVE = Path("docs") / "memory" / "manifest.json"
MANIFEST_KEYS = {"files", "schema_version"}
ENTRY_KEYS = {"bytes", "path", "sha256"}


def default_root() -> Path:
    return Path(__file__).resolve().parents[1]


def git_output(root: Path, *arguments: str) -> str | None:
    try:
        result = subprocess.run(
            ["git", *arguments], cwd=root, check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


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
    return {
        "files": [
            {
                "bytes": path.stat().st_size,
                "path": path.relative_to(root).as_posix(),
                "sha256": sha256(path),
            }
            for path in versionable_memory_paths(root)
        ],
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
    if not isinstance(actual, dict):
        return ["manifest top level must be an object"]

    problems: list[str] = []
    if set(actual) != MANIFEST_KEYS:
        problems.append("manifest top-level keys must be exactly files and schema_version")
        return problems
    if actual["schema_version"] != SCHEMA_VERSION:
        problems.append("manifest schema_version is unsupported")
    if not isinstance(actual["files"], list):
        problems.append("manifest files must be a list")
        return problems
    for entry in actual["files"]:
        if not isinstance(entry, dict) or set(entry) != ENTRY_KEYS:
            problems.append("manifest file entry keys must be exactly bytes, path, and sha256")
            break
    expected = build_manifest(root)
    if actual["files"] != expected["files"]:
        problems.append("memory file inventory, byte count, or SHA-256 differs from manifest")
    return problems


def self_test() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        memory = root / "docs" / "memory"
        memory.mkdir(parents=True)
        source = memory / "SOURCE.md"
        source.write_text("source\n", encoding="utf-8")

        write_manifest(root)
        manifest_path = root / MANIFEST_RELATIVE
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert set(payload) == MANIFEST_KEYS, payload
        assert payload["schema_version"] == 3, payload
        assert all(entry["path"] != MANIFEST_RELATIVE.as_posix() for entry in payload["files"]), payload
        assert not check_manifest(root), "fresh integrity manifest did not verify"

        source.write_text("changed source\n", encoding="utf-8")
        assert check_manifest(root), "file-byte change was not detected"
        write_manifest(root)
        added = memory / "ADDED.md"
        added.write_text("added\n", encoding="utf-8")
        assert check_manifest(root), "file addition was not detected"
        added.unlink()
        assert not check_manifest(root), "file removal back to indexed inventory was not accepted"

        payload = build_manifest(root)
        for key, value in (
            ("active_state", {"memory_schema": 3}),
            ("generated_date_utc", "2026-09-22"),
            ("observed_git_head", "a" * 40),
            ("unexpected", True),
        ):
            legacy = dict(payload)
            legacy[key] = value
            manifest_path.write_bytes(encoded_manifest(legacy))
            assert any("top-level keys" in problem for problem in check_manifest(root)), key
        unsupported = dict(payload)
        unsupported["schema_version"] = 2
        manifest_path.write_bytes(encoded_manifest(unsupported))
        assert any("schema_version" in problem for problem in check_manifest(root))


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
    if args.write:
        write_manifest(root)
        print(f"MANIFEST: WROTE {MANIFEST_RELATIVE.as_posix()}")
        return 0
    problems = check_manifest(root)
    if problems:
        for problem in problems:
            print(f"MANIFEST: FAIL: {problem}")
        return 1
    print("MANIFEST: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
