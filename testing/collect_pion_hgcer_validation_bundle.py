"""Collect a narrow, read-only Phase-C HGCer validation review bundle.

This utility never imports the analysis runtime.  It only locates frozen
artifacts, validates checkpoint metadata, extracts selected PDF pages, hashes
the source bytes, and writes a new ZIP archive.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import datetime as _datetime
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Any, Callable, Iterable, Mapping, Optional, Sequence, Union
import zipfile


COLLECTOR_SCHEMA_VERSION = "pion_hgcer_validation_bundle/v1"
NORMAL_Q4P4W2P74_SETTINGS = (
    ("Left", "lowe"),
    ("Left", "highe"),
    ("Center", "lowe"),
    ("Center", "highe"),
    ("Right", "highe"),
)
_VALID_PHI_SETTINGS = frozenset({"Left", "Center", "Right"})
_VALID_EPSILON_TOKENS = frozenset({"lowe", "highe"})
_CHECKPOINT_SUFFIX = "_kaon_pion-background_hgcer_refinement_checkpoint_"
_HGCER_DEBUG_SUFFIX = "_kaon_rand_sub_"


@dataclass(frozen=True)
class PdfBackend:
    """One optional, non-raster PDF page-copy backend."""

    name: str
    identity: str
    kind: str
    module: Any = None
    executable: Optional[str] = None


CommandRunner = Callable[[Sequence[str], Optional[Path]], Mapping[str, Any]]


def _safe_token(value: object, field: str) -> str:
    if not isinstance(value, str):
        raise ValueError("{}_invalid".format(field))
    token = value.strip()
    if not token or token in {".", ".."} or any(character in token for character in "\\/:"):
        raise ValueError("{}_invalid".format(field))
    return token


def resolve_settings(
    phi: Optional[str] = None, epsilon: Optional[str] = None
) -> tuple[tuple[str, str], ...]:
    """Return one explicit setting or the frozen normal five-setting list."""
    if (phi is None) != (epsilon is None):
        raise ValueError("phi_and_epsilon_must_be_supplied_together")
    if phi is None:
        return NORMAL_Q4P4W2P74_SETTINGS
    normalized_phi = _safe_token(phi, "phi")
    normalized_epsilon = _safe_token(epsilon, "epsilon").lower()
    if normalized_phi not in _VALID_PHI_SETTINGS:
        raise ValueError("phi_invalid")
    if normalized_epsilon not in _VALID_EPSILON_TOKENS:
        raise ValueError("epsilon_invalid")
    return ((normalized_phi, normalized_epsilon),)


def checkpoint_basename(phi: str, kinematic: str, epsilon: str) -> str:
    """Return the existing deterministic Phase-C checkpoint basename."""
    phi, epsilon = resolve_settings(phi, epsilon)[0]
    kinematic = _safe_token(kinematic, "kinematic")
    return "{}{}{}_{}.json".format(phi, _CHECKPOINT_SUFFIX, kinematic, epsilon)


def hgcer_debug_basename(phi: str, kinematic: str, epsilon: str) -> str:
    """Return the existing deterministic HGCer-debug PDF basename."""
    phi, epsilon = resolve_settings(phi, epsilon)[0]
    kinematic = _safe_token(kinematic, "kinematic")
    return "{}{}{}_{}_hgcer-debug.pdf".format(phi, _HGCER_DEBUG_SUFFIX, kinematic, epsilon)


def hgcer_validation_basename(phi: str, kinematic: str, epsilon: str) -> str:
    """Return the slim validation-PDF basename derived from the debug basename."""
    debug_name = hgcer_debug_basename(phi, kinematic, epsilon)
    return "{}_hgcer-validation.pdf".format(debug_name[:-len("_hgcer-debug.pdf")])


def select_validation_pages(page_count: int) -> list[int]:
    """Return the exact, unique 1-based page selection for a slim review PDF."""
    if isinstance(page_count, bool) or not isinstance(page_count, int) or page_count < 1:
        raise ValueError("pdf_page_count_invalid")
    if page_count <= 8:
        return list(range(1, page_count + 1))
    return [1] + list(range(page_count - 6, page_count + 1))


def sha256_file(path: Path) -> str:
    """Return the SHA-256 of file bytes without loading a large PDF at once."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _command_text(command: Sequence[str]) -> str:
    return " ".join(str(item) for item in command)


def run_command(
    command: Sequence[str], cwd: Optional[Path] = None
) -> dict[str, Any]:
    """Run a read-only helper command and retain complete diagnostic output."""
    try:
        completed = subprocess.run(
            [str(item) for item in command],
            cwd=None if cwd is None else os.fspath(cwd),
            capture_output=True,
            text=True,
            errors="replace",
            check=False,
        )
    except OSError as exc:
        return {
            "command": [str(item) for item in command],
            "returncode": 127,
            "stdout": "",
            "stderr": "{}: {}".format(type(exc).__name__, exc),
        }
    return {
        "command": [str(item) for item in command],
        "returncode": int(completed.returncode),
        "stdout": completed.stdout or "",
        "stderr": completed.stderr or "",
    }


def _normalized_command_result(result: Mapping[str, Any], command: Sequence[str]) -> dict[str, Any]:
    return {
        "command": [str(item) for item in result.get("command", command)],
        "returncode": int(result.get("returncode", 127)),
        "stdout": str(result.get("stdout", "") or ""),
        "stderr": str(result.get("stderr", "") or ""),
    }


def _run(
    command_runner: CommandRunner, command: Sequence[str], cwd: Optional[Path]
) -> dict[str, Any]:
    return _normalized_command_result(command_runner(command, cwd), command)


def _module_identity(module_name: str, module: Any) -> str:
    version = getattr(module, "__version__", None)
    return "{}{}".format(module_name, " {}".format(version) if version else "")


def discover_pdf_backends(
    *,
    importer: Callable[[str], Any] = importlib.import_module,
    which: Callable[[str], Optional[str]] = shutil.which,
) -> list[PdfBackend]:
    """Discover only supported non-raster PDF page-copy tools in priority order."""
    backends: list[PdfBackend] = []
    for module_name in ("pypdf", "PyPDF2"):
        try:
            module = importer(module_name)
        except Exception:
            continue
        if hasattr(module, "PdfReader") and hasattr(module, "PdfWriter"):
            backends.append(PdfBackend(
                name=module_name,
                identity=_module_identity(module_name, module),
                kind="python",
                module=module,
            ))
    qpdf = which("qpdf")
    if qpdf:
        backends.append(PdfBackend("qpdf", str(qpdf), "qpdf", executable=str(qpdf)))
    pdfseparate = which("pdfseparate")
    pdfunite = which("pdfunite")
    if pdfseparate and pdfunite:
        backends.append(PdfBackend(
            "pdfseparate_pdfunite",
            "pdfseparate={} pdfunite={}".format(pdfseparate, pdfunite),
            "pdfseparate_pdfunite",
            executable="{}\n{}".format(pdfseparate, pdfunite),
        ))
    mutool = which("mutool")
    if mutool:
        backends.append(PdfBackend("mutool", str(mutool), "mutool", executable=str(mutool)))
    return backends


def _page_range_spec(pages: Iterable[int]) -> str:
    ordered = list(pages)
    if not ordered:
        raise ValueError("pdf_pages_empty")
    ranges: list[str] = []
    start = end = ordered[0]
    for page in ordered[1:]:
        if page == end + 1:
            end = page
            continue
        ranges.append(str(start) if start == end else "{}-{}".format(start, end))
        start = end = page
    ranges.append(str(start) if start == end else "{}-{}".format(start, end))
    return ",".join(ranges)


def _successful_command(
    command_runner: CommandRunner, command: Sequence[str], cwd: Optional[Path]
) -> dict[str, Any]:
    result = _run(command_runner, command, cwd)
    if result["returncode"] != 0:
        raise RuntimeError(
            "{} returned {}: {}".format(
                _command_text(command), result["returncode"], result["stderr"].strip(),
            )
        )
    return result


def _extract_with_python(source: Path, pages: Sequence[int], backend: PdfBackend) -> tuple[int, bytes]:
    reader = backend.module.PdfReader(os.fspath(source))
    page_count = len(reader.pages)
    selected = select_validation_pages(page_count)
    if list(pages) != selected:
        raise ValueError("pdf_page_selection_invalid")
    writer = backend.module.PdfWriter()
    for page_number in selected:
        writer.add_page(reader.pages[page_number - 1])
    output = io.BytesIO()
    writer.write(output)
    return page_count, output.getvalue()


def _extract_with_qpdf(
    source: Path, backend: PdfBackend, command_runner: CommandRunner, temporary: Path,
) -> tuple[int, list[int], bytes]:
    count_result = _successful_command(
        command_runner, [backend.executable, "--show-npages", os.fspath(source)], None,
    )
    try:
        page_count = int(count_result["stdout"].strip())
    except ValueError as exc:
        raise RuntimeError("qpdf_page_count_invalid") from exc
    pages = select_validation_pages(page_count)
    destination = temporary / "slim-qpdf.pdf"
    _successful_command(
        command_runner,
        [backend.executable, "--empty", "--pages", os.fspath(source), _page_range_spec(pages), "--", os.fspath(destination)],
        None,
    )
    return page_count, pages, destination.read_bytes()


def _extract_with_pdfseparate_pdfunite(
    source: Path, backend: PdfBackend, command_runner: CommandRunner, temporary: Path,
) -> tuple[int, list[int], bytes]:
    separator, unite = backend.executable.split("\n", 1)
    pattern = temporary / "page-%d.pdf"
    _successful_command(command_runner, [separator, os.fspath(source), os.fspath(pattern)], None)
    extracted = sorted(
        temporary.glob("page-*.pdf"),
        key=lambda path: int(re.search(r"(\d+)(?=\.pdf$)", path.name).group(1)),
    )
    page_count = len(extracted)
    pages = select_validation_pages(page_count)
    destination = temporary / "slim-pdfunite.pdf"
    _successful_command(
        command_runner,
        [unite] + [os.fspath(extracted[index - 1]) for index in pages] + [os.fspath(destination)],
        None,
    )
    return page_count, pages, destination.read_bytes()


def _extract_with_mutool(
    source: Path, backend: PdfBackend, command_runner: CommandRunner, temporary: Path,
) -> tuple[int, list[int], bytes]:
    info_result = _successful_command(
        command_runner, [backend.executable, "info", os.fspath(source)], None,
    )
    match = re.search(r"^Pages:\s*(\d+)\s*$", info_result["stdout"], re.MULTILINE | re.IGNORECASE)
    if match is None:
        raise RuntimeError("mutool_page_count_invalid")
    page_count = int(match.group(1))
    pages = select_validation_pages(page_count)
    destination = temporary / "slim-mutool.pdf"
    _successful_command(
        command_runner,
        [backend.executable, "merge", "-o", os.fspath(destination), os.fspath(source), _page_range_spec(pages)],
        None,
    )
    return page_count, pages, destination.read_bytes()


def extract_pdf_pages(
    source: Path,
    *,
    backends: Sequence[PdfBackend],
    command_runner: CommandRunner = run_command,
) -> dict[str, Any]:
    """Extract page 1 plus the final seven pages with the first working backend."""
    if not backends:
        return {
            "available": False,
            "error_code": "pdf_extraction_dependency_unavailable",
            "attempts": [],
        }
    attempts: list[dict[str, str]] = []
    for backend in backends:
        try:
            if backend.kind == "python":
                reader = backend.module.PdfReader(os.fspath(source))
                page_count = len(reader.pages)
                pages = select_validation_pages(page_count)
                writer = backend.module.PdfWriter()
                for page_number in pages:
                    writer.add_page(reader.pages[page_number - 1])
                output = io.BytesIO()
                writer.write(output)
                payload = output.getvalue()
            else:
                with tempfile.TemporaryDirectory(prefix="kaonlt-hgcer-bundle-") as temporary_name:
                    temporary = Path(temporary_name)
                    if backend.kind == "qpdf":
                        page_count, pages, payload = _extract_with_qpdf(
                            source, backend, command_runner, temporary
                        )
                    elif backend.kind == "pdfseparate_pdfunite":
                        page_count, pages, payload = _extract_with_pdfseparate_pdfunite(
                            source, backend, command_runner, temporary
                        )
                    elif backend.kind == "mutool":
                        page_count, pages, payload = _extract_with_mutool(
                            source, backend, command_runner, temporary
                        )
                    else:
                        raise RuntimeError("pdf_backend_kind_invalid")
            return {
                "available": True,
                "backend": backend.name,
                "backend_identity": backend.identity,
                "page_count": page_count,
                "pages": pages,
                "payload": payload,
                "attempts": attempts,
            }
        except Exception as exc:
            attempts.append({
                "backend": backend.name,
                "backend_identity": backend.identity,
                "error": "{}: {}".format(type(exc).__name__, exc),
            })
    return {
        "available": False,
        "error_code": "pdf_extraction_failed",
        "attempts": attempts,
    }


def _issue(
    issues: list[dict[str, Any]],
    code: str,
    *,
    setting: Optional[Mapping[str, str]] = None,
    artifact: Optional[str] = None,
    detail: Optional[str] = None,
) -> None:
    entry: dict[str, Any] = {"code": code}
    if setting is not None:
        entry["setting"] = dict(setting)
    if artifact is not None:
        entry["artifact"] = artifact
    if detail:
        entry["detail"] = detail
    issues.append(entry)


def _source_artifact(path: Path, archive_path: str, artifact: str) -> dict[str, Any]:
    result: dict[str, Any] = {
        "artifact": artifact,
        "source_path": os.fspath(path.resolve()),
        "archive_path": archive_path,
        "status": "missing",
        "byte_size": None,
        "sha256": None,
    }
    if path.is_file():
        try:
            result.update({
                "status": "exists",
                "byte_size": path.stat().st_size,
                "sha256": sha256_file(path),
            })
        except OSError as exc:
            result.update({
                "status": "unreadable",
                "error": "{}: {}".format(type(exc).__name__, exc),
            })
    return result


def _checkpoint_metadata(path: Path, setting: Mapping[str, str]) -> tuple[dict[str, Any], list[dict[str, str]]]:
    metadata: dict[str, Any] = {
        "schema_version": None,
        "setting": None,
        "non_authoritative": None,
        "production_objects_mutated": None,
        "refinement_applied": None,
        "metadata_status": "invalid",
    }
    errors: list[dict[str, str]] = []
    try:
        payload = json.loads(
            path.read_text(encoding="utf-8"),
            parse_constant=lambda constant: (_ for _ in ()).throw(
                ValueError("nonstandard JSON constant: {}".format(constant))
            ),
        )
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        errors.append({"code": "checkpoint_json_invalid", "detail": "{}: {}".format(type(exc).__name__, exc)})
        return metadata, errors
    if not isinstance(payload, dict) or not isinstance(payload.get("setting"), dict):
        errors.append({"code": "checkpoint_json_invalid", "detail": "checkpoint setting metadata is absent"})
        return metadata, errors
    checkpoint_setting = payload["setting"]
    metadata.update({
        "schema_version": payload.get("schema_version"),
        "setting": checkpoint_setting,
        "non_authoritative": payload.get("non_authoritative"),
        "production_objects_mutated": payload.get("production_objects_mutated"),
        "refinement_applied": payload.get("refinement_applied"),
        "metadata_status": "match",
    })
    expected = {
        "phi_setting": setting["phi"],
        "epsilon_filename_token": setting["epsilon"],
        "kinematic_token": setting["kinematic"],
        "particle_type": "kaon",
    }
    mismatches = {
        field: {"expected": expected_value, "actual": checkpoint_setting.get(field)}
        for field, expected_value in expected.items()
        if checkpoint_setting.get(field) != expected_value
    }
    if mismatches:
        metadata["metadata_status"] = "mismatch"
        metadata["mismatches"] = mismatches
        errors.append({"code": "checkpoint_metadata_mismatch", "detail": json.dumps(mismatches, sort_keys=True)})
    return metadata, errors


def _format_command_records(records: Sequence[Mapping[str, Any]]) -> str:
    blocks: list[str] = []
    for record in records:
        blocks.extend((
            "$ {}".format(_command_text(record["command"])),
            "returncode: {}".format(record["returncode"]),
            "stdout:",
            record["stdout"].rstrip(),
            "stderr:",
            record["stderr"].rstrip(),
            "",
        ))
    return "\n".join(blocks).rstrip() + "\n"


def collect_source_state(
    repo_root: Path, command_runner: CommandRunner = run_command
) -> tuple[str, Optional[str]]:
    """Capture source identity without changing the repository."""
    records = [
        _run(command_runner, ["git", "rev-parse", "HEAD"], repo_root),
        _run(command_runner, ["git", "status", "--short"], repo_root),
        _run(command_runner, [sys.executable, "--version"], repo_root),
    ]
    head = records[0]["stdout"].strip() if records[0]["returncode"] == 0 else None
    return _format_command_records(records), head or None


def collect_source_checks(repo_root: Path, command_runner: CommandRunner = run_command) -> tuple[str, list[dict[str, Any]]]:
    """Run the required local source checks and retain failures instead of raising."""
    checks = (
        ("py_compile", [sys.executable, "-m", "py_compile", "src/cuts/pion_hgcer_refinement_plots.py", "testing/test_pion_hgcer_refinement_plots.py"]),
        ("hgcer_refinement_plots_unittest", [sys.executable, "-m", "unittest", "testing.test_pion_hgcer_refinement_plots"]),
        ("git_diff_check", ["git", "diff", "--check"]),
    )
    records: list[dict[str, Any]] = []
    for name, command in checks:
        record = _run(command_runner, command, repo_root)
        record["name"] = name
        records.append(record)
    return _format_command_records(records), records


def _archive_json(archive: zipfile.ZipFile, archive_path: str, payload: Mapping[str, Any]) -> None:
    archive.writestr(
        archive_path,
        json.dumps(payload, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n",
    )


def _validate_output_path(output_path: Path, source_paths: Sequence[Path]) -> None:
    if output_path.exists():
        raise ValueError("output_path_already_exists")
    if not output_path.parent.is_dir():
        raise ValueError("output_parent_missing")
    resolved_output = output_path.resolve()
    if resolved_output in {source.resolve() for source in source_paths}:
        raise ValueError("output_path_is_source_artifact")


def collect_validation_bundle(
    *,
    outdir: Union[Path, str],
    kinematic: str,
    output: Union[Path, str],
    phi: Optional[str] = None,
    epsilon: Optional[str] = None,
    repo_root: Optional[Union[Path, str]] = None,
    pdf_backends: Optional[Sequence[PdfBackend]] = None,
    command_runner: CommandRunner = run_command,
) -> dict[str, Any]:
    """Create a best-effort bundle and return its manifest/result status."""
    selected = resolve_settings(phi, epsilon)
    kinematic = _safe_token(kinematic, "kinematic")
    source_root = Path(outdir).expanduser()
    output_path = Path(output).expanduser()
    if not output_path.is_absolute():
        output_path = Path.cwd() / output_path
    output_path = output_path.resolve()
    source_paths = []
    for selected_phi, selected_epsilon in selected:
        source_paths.extend((
            source_root / checkpoint_basename(selected_phi, kinematic, selected_epsilon),
            source_root / hgcer_debug_basename(selected_phi, kinematic, selected_epsilon),
        ))
    _validate_output_path(output_path, source_paths)
    repository = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[1]
    repository = repository.resolve()
    backends = list(discover_pdf_backends() if pdf_backends is None else pdf_backends)
    source_state_text, git_head = collect_source_state(repository, command_runner)
    source_checks_text, source_checks = collect_source_checks(repository, command_runner)
    issues: list[dict[str, Any]] = []
    for check in source_checks:
        if check["returncode"] != 0:
            _issue(issues, "source_check_failed", artifact=check["name"], detail="returncode={}".format(check["returncode"]))

    manifest: dict[str, Any] = {
        "schema_version": COLLECTOR_SCHEMA_VERSION,
        "generated_at_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
        "git_head": git_head,
        "requested_kinematic": kinematic,
        "requested_settings": [
            {"phi": selected_phi, "epsilon": selected_epsilon}
            for selected_phi, selected_epsilon in selected
        ],
        "settings": [],
        "source_checks": source_checks,
        "errors": issues,
        "complete": False,
    }

    staging_descriptor, staging_name = tempfile.mkstemp(
        prefix=".kaonlt-hgcer-validation-", suffix=".zip", dir=output_path.parent
    )
    os.close(staging_descriptor)
    staging_path = Path(staging_name)
    with ExitStack() as cleanup:
        cleanup.callback(staging_path.unlink, missing_ok=True)
        archive = zipfile.ZipFile(staging_path, "w", compression=zipfile.ZIP_DEFLATED)
        cleanup.callback(archive.close)
        for selected_phi, selected_epsilon in selected:
            setting = {"phi": selected_phi, "epsilon": selected_epsilon, "kinematic": kinematic}
            directory = "{}_{}".format(selected_phi, selected_epsilon)
            archive.writestr(directory + "/", b"")
            checkpoint_path = source_root / checkpoint_basename(selected_phi, kinematic, selected_epsilon)
            debug_path = source_root / hgcer_debug_basename(selected_phi, kinematic, selected_epsilon)
            checkpoint_archive = "{}/{}".format(directory, checkpoint_path.name)
            slim_archive = "{}/{}".format(directory, hgcer_validation_basename(selected_phi, kinematic, selected_epsilon))
            setting_manifest: dict[str, Any] = {
                "phi": selected_phi,
                "epsilon": selected_epsilon,
                "artifacts": {},
                "errors": [],
            }

            checkpoint = _source_artifact(checkpoint_path, checkpoint_archive, "phase_c_checkpoint")
            setting_manifest["artifacts"]["phase_c_checkpoint"] = checkpoint
            if checkpoint["status"] != "exists":
                code = "missing_source_artifact" if checkpoint["status"] == "missing" else "source_artifact_unreadable"
                _issue(issues, code, setting=setting, artifact="phase_c_checkpoint", detail=checkpoint.get("error"))
                setting_manifest["errors"].append(code)
            else:
                archive.writestr(checkpoint_archive, checkpoint_path.read_bytes())
                metadata, metadata_errors = _checkpoint_metadata(checkpoint_path, setting)
                checkpoint["checkpoint_metadata"] = metadata
                for error in metadata_errors:
                    _issue(issues, error["code"], setting=setting, artifact="phase_c_checkpoint", detail=error["detail"])
                    setting_manifest["errors"].append(error["code"])

            debug = _source_artifact(debug_path, slim_archive, "hgcer_debug_pdf")
            setting_manifest["artifacts"]["hgcer_debug_pdf"] = debug
            if debug["status"] != "exists":
                code = "missing_source_artifact" if debug["status"] == "missing" else "source_artifact_unreadable"
                _issue(issues, code, setting=setting, artifact="hgcer_debug_pdf", detail=debug.get("error"))
                setting_manifest["errors"].append(code)
            else:
                extraction = extract_pdf_pages(
                    debug_path, backends=backends, command_runner=command_runner
                )
                debug["extraction_attempts"] = extraction.get("attempts", [])
                if extraction["available"]:
                    slim_payload = extraction["payload"]
                    archive.writestr(slim_archive, slim_payload)
                    debug.update({
                        "original_page_count": extraction["page_count"],
                        "extracted_pages": extraction["pages"],
                        "slim_pdf": {
                            "archive_path": slim_archive,
                            "byte_size": len(slim_payload),
                            "sha256": sha256_bytes(slim_payload),
                            "backend": extraction["backend"],
                            "backend_identity": extraction["backend_identity"],
                        },
                    })
                else:
                    debug["extraction_error"] = extraction["error_code"]
                    _issue(issues, extraction["error_code"], setting=setting, artifact="hgcer_debug_pdf")
                    setting_manifest["errors"].append(extraction["error_code"])
            manifest["settings"].append(setting_manifest)

        manifest["complete"] = not issues
        _archive_json(archive, "manifest.json", manifest)
        archive.writestr("source_state.txt", source_state_text)
        archive.writestr("source_checks.txt", source_checks_text)
        archive.close()
        if output_path.exists():
            raise ValueError("output_path_already_exists")
        os.replace(staging_path, output_path)
        cleanup.pop_all()

    return {
        "output_path": os.fspath(output_path),
        "returncode": 0 if manifest["complete"] else 1,
        "manifest": manifest,
    }


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, help="directory containing existing analysis artifacts")
    parser.add_argument("--kinematic", required=True, help="kinematic filename token, for example Q4p4W2p74")
    parser.add_argument("--output", required=True, help="new validation ZIP path")
    parser.add_argument("--phi", help="one phi setting: Left, Center, or Right")
    parser.add_argument("--epsilon", help="one epsilon filename token: lowe or highe")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = build_argument_parser().parse_args(argv)
    try:
        result = collect_validation_bundle(
            outdir=arguments.outdir,
            kinematic=arguments.kinematic,
            output=arguments.output,
            phi=arguments.phi,
            epsilon=arguments.epsilon,
        )
    except ValueError as exc:
        print("collector configuration error: {}".format(exc), file=sys.stderr)
        return 2
    print("wrote {}".format(result["output_path"]))
    if result["returncode"]:
        codes = sorted({entry["code"] for entry in result["manifest"]["errors"]})
        print(
            "bundle completed with validation/collection issues: {}; see manifest.json".format(
                ", ".join(codes)
            ),
            file=sys.stderr,
        )
    return int(result["returncode"])


if __name__ == "__main__":
    raise SystemExit(main())
