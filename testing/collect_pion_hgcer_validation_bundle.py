"""Collect a narrow, read-only Phase-E.7.2 farm-review bundle.

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
import math
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


COLLECTOR_SCHEMA_VERSION = "pion_hgcer_validation_bundle/v2"
PROFILE_SCHEMA_VERSION = "pion_hgcer_validation_bundle_profile/v2"
DEFAULT_PROFILE_PATH = Path(__file__).with_name(
    "pion_hgcer_validation_bundle_profile.json"
)


class PdfPageSelectionError(ValueError):
    """A source PDF cannot contain the required meeting-summary pages."""

    def __init__(self, page_count: int):
        self.page_count = int(page_count)
        super().__init__("pdf_page_count_too_short_for_meeting_summary")


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


def load_validation_profile(
    profile_path: Optional[Union[Path, str]] = None,
) -> dict[str, Any]:
    """Read the declarative bundle profile without importing analysis code."""
    path = DEFAULT_PROFILE_PATH if profile_path is None else Path(profile_path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("validation_bundle_profile_invalid: {}".format(exc)) from exc
    if not isinstance(payload, dict) or payload.get("schema_version") != PROFILE_SCHEMA_VERSION:
        raise ValueError("validation_bundle_profile_invalid")
    if not isinstance(payload.get("validation_profile"), str):
        raise ValueError("validation_bundle_profile_invalid")
    settings = payload.get("settings")
    if not isinstance(settings, list) or not settings:
        raise ValueError("validation_bundle_profile_invalid")
    for setting in settings:
        if not isinstance(setting, dict):
            raise ValueError("validation_bundle_profile_invalid")
        _safe_token(setting.get("phi"), "profile_phi")
        _safe_token(setting.get("epsilon"), "profile_epsilon")
    artifacts = payload.get("artifacts")
    if not isinstance(artifacts, dict):
        raise ValueError("validation_bundle_profile_invalid")
    checkpoint = artifacts.get("phase_c_checkpoint")
    phase_d = artifacts.get("phase_d_checkpoint")
    correction = artifacts.get("parent_preserving_correction")
    procedure = artifacts.get("procedure_pdf")
    if (
        not isinstance(checkpoint, dict)
        or not isinstance(phase_d, dict)
        or not isinstance(correction, dict)
        or not isinstance(procedure, dict)
    ):
        raise ValueError("validation_bundle_profile_invalid")
    for artifact, field in (
        (checkpoint, "basename_template"),
        (phase_d, "basename_template"),
        (correction, "basename_template"),
        (procedure, "source_basename_template"),
        (procedure, "page_manifest_basename_template"),
        (procedure, "slim_basename_template"),
    ):
        if not isinstance(artifact.get(field), str):
            raise ValueError("validation_bundle_profile_invalid")
    selection = procedure.get("page_selection")
    if (
        not isinstance(selection, dict)
        or selection.get("kind") != "final_pages"
        or isinstance(selection.get("count"), bool)
        or not isinstance(selection.get("count"), int)
        or selection["count"] < 1
    ):
        raise ValueError("validation_bundle_profile_invalid")
    source_identity = payload.get("source_identity")
    if not isinstance(source_identity, dict):
        raise ValueError("validation_bundle_profile_invalid")
    required_commit = source_identity.get("required_analysis_commit")
    allowed_files = source_identity.get("allowed_committed_files")
    if (
        not isinstance(required_commit, str)
        or re.fullmatch(r"[0-9a-f]{40}", required_commit) is None
        or not isinstance(allowed_files, list)
        or not all(isinstance(path, str) and path for path in allowed_files)
    ):
        raise ValueError("validation_bundle_profile_invalid")
    return json.loads(json.dumps(payload, ensure_ascii=True, allow_nan=False))


_DEFAULT_PROFILE = load_validation_profile()
VALIDATION_PROFILE = _DEFAULT_PROFILE["validation_profile"]
REQUIRED_ANALYSIS_COMMIT = _DEFAULT_PROFILE["source_identity"]["required_analysis_commit"]
AUTHORIZED_SETTINGS = tuple(
    (setting["phi"], setting["epsilon"])
    for setting in _DEFAULT_PROFILE["settings"]
)
ALLOWED_COMMITTED_FILES = frozenset(
    _DEFAULT_PROFILE["source_identity"]["allowed_committed_files"]
)


def _format_profile_basename(
    template: str, phi: str, kinematic: str, epsilon: str,
) -> str:
    """Format one profile-owned basename and reject path traversal."""
    try:
        basename = template.format(phi=phi, kinematic=kinematic, epsilon=epsilon)
    except (KeyError, ValueError) as exc:
        raise ValueError("validation_bundle_profile_invalid") from exc
    if Path(basename).name != basename or not basename:
        raise ValueError("validation_bundle_profile_invalid")
    return basename


def resolve_settings(
    phi: Optional[str] = None,
    epsilon: Optional[str] = None,
    profile: Optional[Mapping[str, Any]] = None,
) -> tuple[tuple[str, str], ...]:
    """Return the settings declared by the active bundle profile."""
    active = _DEFAULT_PROFILE if profile is None else profile
    authorized = tuple(
        (str(setting["phi"]), str(setting["epsilon"]))
        for setting in active["settings"]
    )
    if (phi is None) != (epsilon is None):
        raise ValueError("phi_and_epsilon_must_be_supplied_together")
    if phi is None:
        return authorized
    normalized_phi = _safe_token(phi, "phi")
    normalized_epsilon = _safe_token(epsilon, "epsilon").lower()
    if (normalized_phi, normalized_epsilon) not in authorized:
        raise ValueError("validation_bundle_profile_setting_not_authorized")
    return ((normalized_phi, normalized_epsilon),)


def checkpoint_basename(
    phi: str, kinematic: str, epsilon: str,
    profile: Optional[Mapping[str, Any]] = None,
) -> str:
    """Return the existing deterministic Phase-C checkpoint basename."""
    active = _DEFAULT_PROFILE if profile is None else profile
    phi, epsilon = resolve_settings(phi, epsilon, active)[0]
    kinematic = _safe_token(kinematic, "kinematic")
    return _format_profile_basename(
        active["artifacts"]["phase_c_checkpoint"]["basename_template"],
        phi, kinematic, epsilon,
    )


def full_background_subtraction_basename(
    phi: str, kinematic: str, epsilon: str,
    profile: Optional[Mapping[str, Any]] = None,
) -> str:
    """Return the deterministic full-background procedure-PDF basename."""
    active = _DEFAULT_PROFILE if profile is None else profile
    phi, epsilon = resolve_settings(phi, epsilon, active)[0]
    kinematic = _safe_token(kinematic, "kinematic")
    return _format_profile_basename(
        active["artifacts"]["procedure_pdf"]["source_basename_template"],
        phi, kinematic, epsilon,
    )


def phase_d_checkpoint_basename(
    phi: str, kinematic: str, epsilon: str,
    profile: Optional[Mapping[str, Any]] = None,
) -> str:
    """Return the deterministic Phase-D checkpoint basename."""
    active = _DEFAULT_PROFILE if profile is None else profile
    phi, epsilon = resolve_settings(phi, epsilon, active)[0]
    return _format_profile_basename(
        active["artifacts"]["phase_d_checkpoint"]["basename_template"],
        phi, _safe_token(kinematic, "kinematic"), epsilon,
    )


def parent_preserving_correction_basename(
    phi: str, kinematic: str, epsilon: str,
    profile: Optional[Mapping[str, Any]] = None,
) -> str:
    """Return the deterministic E.7.1 correction-artifact basename."""
    active = _DEFAULT_PROFILE if profile is None else profile
    phi, epsilon = resolve_settings(phi, epsilon, active)[0]
    return _format_profile_basename(
        active["artifacts"]["parent_preserving_correction"]["basename_template"],
        phi, _safe_token(kinematic, "kinematic"), epsilon,
    )


def full_background_page_manifest_basename(
    phi: str, kinematic: str, epsilon: str,
    profile: Optional[Mapping[str, Any]] = None,
) -> str:
    """Return the renderer-owned full-background page-manifest basename."""
    active = _DEFAULT_PROFILE if profile is None else profile
    phi, epsilon = resolve_settings(phi, epsilon, active)[0]
    return _format_profile_basename(
        active["artifacts"]["procedure_pdf"]["page_manifest_basename_template"],
        phi, _safe_token(kinematic, "kinematic"), epsilon,
    )


def meeting_summary_basename(
    phi: str, kinematic: str, epsilon: str,
    profile: Optional[Mapping[str, Any]] = None,
) -> str:
    """Return the deterministic four-page meeting-summary PDF basename."""
    active = _DEFAULT_PROFILE if profile is None else profile
    phi, epsilon = resolve_settings(phi, epsilon, active)[0]
    kinematic = _safe_token(kinematic, "kinematic")
    return _format_profile_basename(
        active["artifacts"]["procedure_pdf"]["slim_basename_template"],
        phi, kinematic, epsilon,
    )


def select_validation_pages(
    page_count: int,
    page_selection: Optional[Mapping[str, Any]] = None,
) -> list[int]:
    """Return the required final 1-based meeting-summary PDF pages."""
    if isinstance(page_count, bool) or not isinstance(page_count, int) or page_count < 1:
        raise ValueError("pdf_page_count_invalid")
    selection = (
        _DEFAULT_PROFILE["artifacts"]["procedure_pdf"]["page_selection"]
        if page_selection is None else page_selection
    )
    if selection.get("kind") != "final_pages" or not isinstance(selection.get("count"), int):
        raise ValueError("pdf_page_selection_invalid")
    count = selection["count"]
    if page_count < count:
        raise PdfPageSelectionError(page_count)
    return list(range(page_count - count + 1, page_count + 1))


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
    page_selection: Mapping[str, Any],
) -> tuple[int, list[int], bytes]:
    count_result = _successful_command(
        command_runner, [backend.executable, "--show-npages", os.fspath(source)], None,
    )
    try:
        page_count = int(count_result["stdout"].strip())
    except ValueError as exc:
        raise RuntimeError("qpdf_page_count_invalid") from exc
    pages = select_validation_pages(page_count, page_selection)
    destination = temporary / "slim-qpdf.pdf"
    _successful_command(
        command_runner,
        [backend.executable, "--empty", "--pages", os.fspath(source), _page_range_spec(pages), "--", os.fspath(destination)],
        None,
    )
    return page_count, pages, destination.read_bytes()


def _extract_with_pdfseparate_pdfunite(
    source: Path, backend: PdfBackend, command_runner: CommandRunner, temporary: Path,
    page_selection: Mapping[str, Any],
) -> tuple[int, list[int], bytes]:
    separator, unite = backend.executable.split("\n", 1)
    pattern = temporary / "page-%d.pdf"
    _successful_command(command_runner, [separator, os.fspath(source), os.fspath(pattern)], None)
    extracted = sorted(
        temporary.glob("page-*.pdf"),
        key=lambda path: int(re.search(r"(\d+)(?=\.pdf$)", path.name).group(1)),
    )
    page_count = len(extracted)
    pages = select_validation_pages(page_count, page_selection)
    destination = temporary / "slim-pdfunite.pdf"
    _successful_command(
        command_runner,
        [unite] + [os.fspath(extracted[index - 1]) for index in pages] + [os.fspath(destination)],
        None,
    )
    return page_count, pages, destination.read_bytes()


def _extract_with_mutool(
    source: Path, backend: PdfBackend, command_runner: CommandRunner, temporary: Path,
    page_selection: Mapping[str, Any],
) -> tuple[int, list[int], bytes]:
    info_result = _successful_command(
        command_runner, [backend.executable, "info", os.fspath(source)], None,
    )
    match = re.search(r"^Pages:\s*(\d+)\s*$", info_result["stdout"], re.MULTILINE | re.IGNORECASE)
    if match is None:
        raise RuntimeError("mutool_page_count_invalid")
    page_count = int(match.group(1))
    pages = select_validation_pages(page_count, page_selection)
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
    page_selection: Optional[Mapping[str, Any]] = None,
) -> dict[str, Any]:
    """Extract the final meeting-summary pages with one non-raster backend."""
    if not backends:
        return {
            "available": False,
            "error_code": "pdf_extraction_dependency_unavailable",
            "attempts": [],
        }
    active_selection = (
        _DEFAULT_PROFILE["artifacts"]["procedure_pdf"]["page_selection"]
        if page_selection is None else page_selection
    )
    attempts: list[dict[str, str]] = []
    for backend in backends:
        try:
            if backend.kind == "python":
                reader = backend.module.PdfReader(os.fspath(source))
                page_count = len(reader.pages)
                pages = select_validation_pages(page_count, active_selection)
                writer = backend.module.PdfWriter()
                for page_number in pages:
                    writer.add_page(reader.pages[page_number - 1])
                output = io.BytesIO()
                writer.write(output)
                payload = output.getvalue()
            else:
                with tempfile.TemporaryDirectory(prefix="kaonlt-e3-bundle-") as temporary_name:
                    temporary = Path(temporary_name)
                    if backend.kind == "qpdf":
                        page_count, pages, payload = _extract_with_qpdf(
                            source, backend, command_runner, temporary, active_selection
                        )
                    elif backend.kind == "pdfseparate_pdfunite":
                        page_count, pages, payload = _extract_with_pdfseparate_pdfunite(
                            source, backend, command_runner, temporary, active_selection
                        )
                    elif backend.kind == "mutool":
                        page_count, pages, payload = _extract_with_mutool(
                            source, backend, command_runner, temporary, active_selection
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
        except PdfPageSelectionError as exc:
            attempts.append({
                "backend": backend.name,
                "backend_identity": backend.identity,
                "error": "{}: {}".format(type(exc).__name__, exc),
            })
            return {
                "available": False,
                "error_code": str(exc),
                "page_count": exc.page_count,
                "pages": [],
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


def _source_artifact(
    path: Path, archive_path: Optional[str], artifact: str,
) -> dict[str, Any]:
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


def _strict_json_payload(path: Path) -> Mapping[str, Any]:
    """Load a JSON artifact while rejecting nonstandard numeric constants."""
    payload = json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=lambda constant: (_ for _ in ()).throw(
            ValueError("nonstandard JSON constant: {}".format(constant))
        ),
    )
    if not isinstance(payload, Mapping):
        raise ValueError("artifact root is not an object")
    return payload


def _setting_mismatches(
    artifact_setting: object, setting: Mapping[str, str],
) -> dict[str, dict[str, object]]:
    if not isinstance(artifact_setting, Mapping):
        return {"setting": {"expected": "mapping", "actual": artifact_setting}}
    expected = {
        "phi_setting": setting["phi"],
        "epsilon_filename_token": setting["epsilon"],
        "kinematic_token": setting["kinematic"],
        "particle_type": "kaon",
    }
    return {
        field: {"expected": expected_value, "actual": artifact_setting.get(field)}
        for field, expected_value in expected.items()
        if artifact_setting.get(field) != expected_value
    }


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
        payload = _strict_json_payload(path)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        errors.append({"code": "checkpoint_json_invalid", "detail": "{}: {}".format(type(exc).__name__, exc)})
        return metadata, errors
    if not isinstance(payload.get("setting"), Mapping):
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
    mismatches = _setting_mismatches(checkpoint_setting, setting)
    if mismatches:
        metadata["metadata_status"] = "mismatch"
        metadata["mismatches"] = mismatches
        errors.append({"code": "checkpoint_metadata_mismatch", "detail": json.dumps(mismatches, sort_keys=True)})
    return metadata, errors


def _phase_d_metadata(
    path: Path, setting: Mapping[str, str],
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Validate the detached Phase-D authority and setting contract."""
    metadata: dict[str, Any] = {"metadata_status": "invalid"}
    errors: list[dict[str, str]] = []
    try:
        payload = _strict_json_payload(path)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        return metadata, [{"code": "checkpoint_json_invalid", "detail": "{}: {}".format(type(exc).__name__, exc)}]
    metadata.update({"schema_version": payload.get("schema_version"), "setting": payload.get("setting")})
    mismatches = _setting_mismatches(payload.get("setting"), setting)
    flags_valid = (
        payload.get("schema_version") == "pion_hgcer_phase_d_checkpoint/v1"
        and payload.get("non_authoritative") is True
        and payload.get("decision_performed") is False
        and payload.get("statistical_compatibility_claimed") is False
        and payload.get("production_objects_mutated") is False
        and payload.get("refinement_applied") is False
    )
    if mismatches or not flags_valid:
        metadata["metadata_status"] = "mismatch"
        if mismatches:
            metadata["mismatches"] = mismatches
        errors.append({
            "code": "checkpoint_metadata_mismatch",
            "detail": "phase_d_contract_invalid" if not mismatches else json.dumps(mismatches, sort_keys=True),
        })
    else:
        metadata["metadata_status"] = "match"
    return metadata, errors


def _parent_preserving_metadata(
    path: Path, setting: Mapping[str, str],
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Validate E.7.1's frozen wrapper, correction, lattice, and closures."""
    metadata: dict[str, Any] = {"metadata_status": "invalid"}
    try:
        payload = _strict_json_payload(path)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        return metadata, [{"code": "checkpoint_json_invalid", "detail": "{}: {}".format(type(exc).__name__, exc)}]
    correction = payload.get("correction")
    mismatches = _setting_mismatches(payload.get("setting"), setting)
    wrapper_flags = all(payload.get(key) is expected for key, expected in (
        ("non_authoritative", True),
        ("production_objects_mutated", False),
        ("refinement_applied", False),
        ("production_application_performed", False),
        ("event_application_performed", False),
    ))
    correction_flags = isinstance(correction, Mapping) and all(correction.get(key) is expected for key, expected in (
        ("available", True),
        ("non_authoritative", True),
        ("production_objects_mutated", False),
        ("refinement_applied", False),
        ("production_application_performed", False),
        ("event_application_performed", False),
    ))
    parents = correction.get("parents") if isinstance(correction, Mapping) else None
    cells = correction.get("cells") if isinstance(correction, Mapping) else None
    t_edges = correction.get("t_edges") if isinstance(correction, Mapping) else None
    delta_edges = correction.get("delta_edges") if isinstance(correction, Mapping) else None
    geometry_valid = all(
        isinstance(edges, list)
        and len(edges) >= 2
        and all(
            isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value))
            for value in edges
        )
        and all(float(edges[index]) < float(edges[index + 1]) for index in range(len(edges) - 1))
        for edges in (t_edges, delta_edges)
    )
    expected_parent_indices = set(range(len(t_edges) - 1)) if geometry_valid else set()
    observed_parent_indices: set[int] = set()
    parents_valid = isinstance(parents, list) and bool(parents) and geometry_valid
    if parents_valid:
        for parent in parents:
            if not isinstance(parent, Mapping):
                parents_valid = False
                break
            t_index = parent.get("t_index")
            refinable_cell_count = parent.get("refinable_cell_count")
            if (
                isinstance(t_index, bool) or not isinstance(t_index, int)
                or t_index not in expected_parent_indices
                or t_index in observed_parent_indices
                or isinstance(parent.get("t_low"), bool)
                or not isinstance(parent.get("t_low"), (int, float))
                or not math.isfinite(float(parent.get("t_low")))
                or isinstance(parent.get("t_high"), bool)
                or not isinstance(parent.get("t_high"), (int, float))
                or not math.isfinite(float(parent.get("t_high")))
                or parent.get("t_low") != t_edges[t_index]
                or parent.get("t_high") != t_edges[t_index + 1]
                or not isinstance(parent.get("parent_status"), str)
                or not parent.get("parent_status")
                or parent.get("closure_passed") is not True
                or isinstance(refinable_cell_count, bool)
                or not isinstance(refinable_cell_count, int)
                or refinable_cell_count < 0
            ):
                parents_valid = False
                break
            observed_parent_indices.add(t_index)
        parents_valid = parents_valid and observed_parent_indices == expected_parent_indices
    expected_keys = {
        (t_index, delta_index)
        for t_index in range(len(t_edges) - 1)
        for delta_index in range(len(delta_edges) - 1)
    } if geometry_valid else set()
    observed_keys: set[tuple[int, int]] = set()
    cells_valid = isinstance(cells, list) and bool(cells) and geometry_valid
    if cells_valid:
        for cell in cells:
            if not isinstance(cell, Mapping):
                cells_valid = False
                break
            t_index = cell.get("t_index")
            delta_index = cell.get("delta_index")
            if (
                isinstance(t_index, bool) or not isinstance(t_index, int)
                or isinstance(delta_index, bool) or not isinstance(delta_index, int)
                or (t_index, delta_index) in observed_keys
            ):
                cells_valid = False
                break
            observed_keys.add((t_index, delta_index))
        cells_valid = cells_valid and observed_keys == expected_keys
    correction_valid = (
        correction_flags
        and correction.get("schema_version") == "pion_hgcer_parent_preserving_correction/v1"
        and correction.get("status") == "available"
        and isinstance(correction.get("fingerprint"), str)
        and bool(correction["fingerprint"])
        and parents_valid
        and cells_valid
    )
    metadata.update({
        "schema_version": payload.get("schema_version"),
        "setting": payload.get("setting"),
        "correction_fingerprint": correction.get("fingerprint") if isinstance(correction, Mapping) else None,
        "parent_statuses": [parent.get("parent_status") for parent in parents] if isinstance(parents, list) else [],
        "parent_closure_passed": [parent.get("closure_passed") for parent in parents] if isinstance(parents, list) else [],
        "refinable_cell_counts": [parent.get("refinable_cell_count") for parent in parents] if isinstance(parents, list) else [],
    })
    if (
        payload.get("schema_version") != "pion_hgcer_parent_preserving_correction_artifact/v1"
        or not wrapper_flags
        or mismatches
        or not correction_valid
    ):
        metadata["metadata_status"] = "mismatch"
        if mismatches:
            metadata["mismatches"] = mismatches
        return metadata, [{"code": "checkpoint_metadata_mismatch", "detail": "parent_preserving_correction_contract_invalid"}]
    metadata["metadata_status"] = "match"
    return metadata, []


_E72_PAGE_IDS = (
    "full_background.e72.meeting_overview",
    "full_background.e72.ab_evidence_summary",
    "full_background.e72.parent_preserving_map_summary",
    "full_background.e72.parent_closure_summary",
)


def _page_manifest_metadata(
    path: Path, setting: Mapping[str, str], pdf_basename: str,
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Validate the renderer-owned sidecar before accepting final-page extraction."""
    metadata: dict[str, Any] = {"metadata_status": "invalid"}
    try:
        payload = _strict_json_payload(path)
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        return metadata, [{"code": "checkpoint_json_invalid", "detail": "{}: {}".format(type(exc).__name__, exc)}]
    pages = payload.get("pages")
    failures = payload.get("renderer_failures")
    mismatches = _setting_mismatches(payload.get("setting"), setting)
    terminal = pages[-4:] if isinstance(pages, list) else []
    terminal_valid = len(terminal) == 4 and all(
        isinstance(entry, Mapping)
        and entry.get("page_id") == page_id
        and entry.get("scope") == "setting"
        and entry.get("authoritative") is False
        for entry, page_id in zip(terminal, _E72_PAGE_IDS)
    )
    e72_failure = isinstance(failures, list) and any(
        str(failure).startswith("E.7.2:") for failure in failures
    )
    valid = (
        payload.get("schema_version") == "full_background_subtraction_page_manifest/v1"
        and not mismatches
        and payload.get("pdf_basename") == pdf_basename
        and isinstance(pages, list)
        and isinstance(failures, list)
        and terminal_valid
        and not e72_failure
    )
    metadata.update({
        "schema_version": payload.get("schema_version"),
        "setting": payload.get("setting"),
        "pdf_basename": payload.get("pdf_basename"),
        "page_count": len(pages) if isinstance(pages, list) else None,
        "renderer_failures": failures if isinstance(failures, list) else None,
        "metadata_status": "match" if valid else "mismatch",
    })
    if not valid:
        return metadata, [{"code": "checkpoint_metadata_mismatch", "detail": "full_background_page_manifest_contract_invalid"}]
    return metadata, []


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
        _run(command_runner, ["git", "rev-parse", "HEAD^"], repo_root),
        _run(command_runner, ["git", "status", "--short"], repo_root),
        _run(command_runner, [sys.executable, "--version"], repo_root),
    ]
    head = records[0]["stdout"].strip() if records[0]["returncode"] == 0 else None
    return _format_command_records(records), head or None


def collect_source_checks(
    repo_root: Path,
    command_runner: CommandRunner = run_command,
    *,
    required_analysis_commit: str = REQUIRED_ANALYSIS_COMMIT,
) -> tuple[str, list[dict[str, Any]]]:
    """Run the detached E.7.2 farm-gate source and identity checks."""
    checks = (
        (
            "required_analysis_commit_ancestor",
            [
                "git", "merge-base", "--is-ancestor",
                required_analysis_commit, "HEAD",
            ],
        ),
        (
            "committed_files_after_required_analysis_commit",
            [
                "git", "diff", "--name-only",
                "{}..HEAD".format(required_analysis_commit),
            ],
        ),
        (
            "py_compile",
            [
                sys.executable, "-m", "py_compile",
                "src/cuts/full_background_subtraction_plots.py",
                "src/cuts/rand_sub.py",
                "src/cuts/pion_hgcer_parent_preserving_correction.py",
                "testing/test_full_background_subtraction_plots.py",
                "testing/test_pion_hgcer_parent_preserving_correction.py",
                "testing/test_pion_hgcer_phase_e_runtime_contract.py",
                "testing/collect_pion_hgcer_validation_bundle.py",
                "testing/test_collect_pion_hgcer_validation_bundle.py",
            ],
        ),
        (
            "full_background_subtraction_plots_unittest",
            [sys.executable, "-m", "unittest", "testing.test_full_background_subtraction_plots"],
        ),
        (
            "parent_preserving_correction_unittest",
            [sys.executable, "-m", "unittest", "testing.test_pion_hgcer_parent_preserving_correction"],
        ),
        (
            "phase_e_runtime_contract_unittest",
            [sys.executable, "-m", "unittest", "testing.test_pion_hgcer_phase_e_runtime_contract"],
        ),
        (
            "ab_combination_prototype_unittest",
            [sys.executable, "-m", "unittest", "testing.test_pion_hgcer_ab_combination_prototype"],
        ),
        (
            "pion_hgcer_event_contract_unittest",
            [sys.executable, "-m", "unittest", "testing.test_pion_hgcer_event_contract"],
        ),
        (
            "validation_bundle_collector_unittest",
            [
                sys.executable, "-m", "unittest",
                "testing.test_collect_pion_hgcer_validation_bundle",
            ],
        ),
        ("git_diff_check", ["git", "diff", "--check"]),
        (
            "git_diff_check_required_analysis_commit_range",
            ["git", "diff", "--check", "{}..HEAD".format(required_analysis_commit)],
        ),
    )
    records: list[dict[str, Any]] = []
    for name, command in checks:
        record = _run(command_runner, command, repo_root)
        record["name"] = name
        records.append(record)
    return _format_command_records(records), records


def _committed_identity(
    source_checks: Sequence[Mapping[str, Any]],
    allowed_committed_files: Iterable[str] = ALLOWED_COMMITTED_FILES,
) -> tuple[bool, list[str], list[str]]:
    """Return the frozen-base ancestry result and committed post-base paths."""
    by_name = {str(record.get("name")): record for record in source_checks}
    ancestry = by_name.get("required_analysis_commit_ancestor", {})
    diff = by_name.get("committed_files_after_required_analysis_commit", {})
    is_ancestor = int(ancestry.get("returncode", 127)) == 0
    if int(diff.get("returncode", 127)) != 0:
        return is_ancestor, [], []
    committed_files = sorted({
        line.strip() for line in str(diff.get("stdout", "")).splitlines()
        if line.strip()
    })
    allowed = frozenset(allowed_committed_files)
    unexpected = [
        path for path in committed_files
        if path not in allowed
    ]
    return is_ancestor, committed_files, unexpected


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
    profile_path: Optional[Union[Path, str]] = None,
    repo_root: Optional[Union[Path, str]] = None,
    pdf_backends: Optional[Sequence[PdfBackend]] = None,
    command_runner: CommandRunner = run_command,
) -> dict[str, Any]:
    """Create a best-effort bundle and return its manifest/result status."""
    profile = load_validation_profile(profile_path)
    selected = resolve_settings(phi, epsilon, profile)
    kinematic = _safe_token(kinematic, "kinematic")
    source_root = Path(outdir).expanduser()
    output_path = Path(output).expanduser()
    if not output_path.is_absolute():
        output_path = Path.cwd() / output_path
    output_path = output_path.resolve()
    source_paths = []
    for selected_phi, selected_epsilon in selected:
        source_paths.extend((
            source_root / checkpoint_basename(
                selected_phi, kinematic, selected_epsilon, profile
            ),
            source_root / phase_d_checkpoint_basename(
                selected_phi, kinematic, selected_epsilon, profile
            ),
            source_root / parent_preserving_correction_basename(
                selected_phi, kinematic, selected_epsilon, profile
            ),
            source_root / full_background_subtraction_basename(
                selected_phi, kinematic, selected_epsilon, profile
            ),
            source_root / full_background_page_manifest_basename(
                selected_phi, kinematic, selected_epsilon, profile
            ),
        ))
    _validate_output_path(output_path, source_paths)
    repository = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[1]
    repository = repository.resolve()
    backends = list(discover_pdf_backends() if pdf_backends is None else pdf_backends)
    source_state_text, git_head = collect_source_state(repository, command_runner)
    source_identity = profile["source_identity"]
    required_analysis_commit = source_identity["required_analysis_commit"]
    allowed_committed_files = source_identity["allowed_committed_files"]
    source_checks_text, source_checks = collect_source_checks(
        repository,
        command_runner,
        required_analysis_commit=required_analysis_commit,
    )
    issues: list[dict[str, Any]] = []
    for check in source_checks:
        if check["returncode"] != 0:
            _issue(issues, "source_check_failed", artifact=check["name"], detail="returncode={}".format(check["returncode"]))
    required_analysis_commit_is_ancestor, committed_files, unexpected_committed_files = (
        _committed_identity(source_checks, allowed_committed_files)
    )
    if not required_analysis_commit_is_ancestor:
        _issue(
            issues,
            "required_analysis_commit_not_present",
            artifact="required_analysis_commit_ancestor",
            detail=required_analysis_commit,
        )
    if unexpected_committed_files:
        _issue(
            issues,
            "unexpected_committed_files_after_required_analysis_commit",
            artifact="committed_files_after_required_analysis_commit",
            detail=json.dumps(unexpected_committed_files, sort_keys=True),
        )

    manifest: dict[str, Any] = {
        "schema_version": COLLECTOR_SCHEMA_VERSION,
        "validation_profile": profile["validation_profile"],
        "validation_profile_source": os.fspath(
            (DEFAULT_PROFILE_PATH if profile_path is None else Path(profile_path)).resolve()
        ),
        "generated_at_utc": _datetime.datetime.now(_datetime.timezone.utc).isoformat().replace("+00:00", "Z"),
        "git_head": git_head,
        "required_analysis_commit": required_analysis_commit,
        "required_analysis_commit_is_ancestor": required_analysis_commit_is_ancestor,
        "committed_files_after_required_analysis_commit": committed_files,
        "unexpected_committed_files_after_required_analysis_commit": unexpected_committed_files,
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
        prefix=".kaonlt-e7-review-", suffix=".zip", dir=output_path.parent
    )
    os.close(staging_descriptor)
    staging_path = Path(staging_name)
    with ExitStack() as cleanup:
        cleanup.callback(staging_path.unlink, missing_ok=True)
        archive = zipfile.ZipFile(staging_path, "w", compression=zipfile.ZIP_DEFLATED)
        cleanup.callback(archive.close)
        def record_json_artifact(
            setting_manifest: dict[str, Any], setting: Mapping[str, str],
            key: str, path: Path, validator: Callable[[Path, Mapping[str, str]], tuple[dict[str, Any], list[dict[str, str]]]],
            directory: str,
        ) -> bool:
            archive_path = "{}/{}".format(directory, path.name)
            record = _source_artifact(path, archive_path, key)
            setting_manifest["artifacts"][key] = record
            if record["status"] != "exists":
                code = "missing_source_artifact" if record["status"] == "missing" else "source_artifact_unreadable"
                _issue(issues, code, setting=setting, artifact=key, detail=record.get("error"))
                setting_manifest["errors"].append(code)
                return False
            archive.writestr(archive_path, path.read_bytes())
            metadata, metadata_errors = validator(path, setting)
            record["metadata"] = metadata
            for error in metadata_errors:
                _issue(issues, error["code"], setting=setting, artifact=key, detail=error["detail"])
                setting_manifest["errors"].append(error["code"])
            return not metadata_errors

        for selected_phi, selected_epsilon in selected:
            setting = {"phi": selected_phi, "epsilon": selected_epsilon, "kinematic": kinematic}
            directory = "{}_{}".format(selected_phi, selected_epsilon)
            archive.writestr(directory + "/", b"")
            checkpoint_path = source_root / checkpoint_basename(
                selected_phi, kinematic, selected_epsilon, profile
            )
            phase_d_path = source_root / phase_d_checkpoint_basename(
                selected_phi, kinematic, selected_epsilon, profile
            )
            correction_path = source_root / parent_preserving_correction_basename(
                selected_phi, kinematic, selected_epsilon, profile
            )
            full_background_path = source_root / full_background_subtraction_basename(
                selected_phi, kinematic, selected_epsilon, profile
            )
            page_manifest_path = source_root / full_background_page_manifest_basename(
                selected_phi, kinematic, selected_epsilon, profile
            )
            slim_archive = "{}/{}".format(
                directory, meeting_summary_basename(
                    selected_phi, kinematic, selected_epsilon, profile
                )
            )
            setting_manifest: dict[str, Any] = {
                "phi": selected_phi,
                "epsilon": selected_epsilon,
                "kinematic": kinematic,
                "artifacts": {},
                "errors": [],
            }

            record_json_artifact(
                setting_manifest, setting, "phase_c_checkpoint", checkpoint_path,
                _checkpoint_metadata, directory,
            )
            record_json_artifact(
                setting_manifest, setting, "phase_d_checkpoint", phase_d_path,
                _phase_d_metadata, directory,
            )
            record_json_artifact(
                setting_manifest, setting, "parent_preserving_correction", correction_path,
                _parent_preserving_metadata, directory,
            )
            page_manifest_valid = record_json_artifact(
                setting_manifest, setting, "full_background_page_manifest", page_manifest_path,
                lambda path, selected_setting: _page_manifest_metadata(
                    path, selected_setting, full_background_path.name,
                ), directory,
            )

            full_background = _source_artifact(
                full_background_path, None, "full_background_subtraction_pdf"
            )
            setting_manifest["artifacts"]["full_background_subtraction_pdf"] = full_background
            meeting_summary: dict[str, Any] = {
                "artifact": "meeting_summary_pdf",
                "archive_path": slim_archive,
                "status": "unavailable",
                "byte_size": None,
                "sha256": None,
            }
            setting_manifest["artifacts"]["meeting_summary_pdf"] = meeting_summary
            if full_background["status"] != "exists":
                code = (
                    "missing_source_artifact"
                    if full_background["status"] == "missing"
                    else "source_artifact_unreadable"
                )
                _issue(
                    issues, code, setting=setting,
                    artifact="full_background_subtraction_pdf",
                    detail=full_background.get("error"),
                )
                setting_manifest["errors"].append(code)
            elif not page_manifest_valid:
                meeting_summary["error"] = "full_background_page_manifest_invalid"
            else:
                extraction = extract_pdf_pages(
                    full_background_path,
                    backends=backends,
                    command_runner=command_runner,
                    page_selection=profile["artifacts"]["procedure_pdf"]["page_selection"],
                )
                full_background["extraction_attempts"] = extraction.get("attempts", [])
                if extraction.get("page_count") is not None:
                    full_background["original_page_count"] = extraction["page_count"]
                if "pages" in extraction:
                    full_background["extracted_pages"] = extraction["pages"]
                page_manifest_record = setting_manifest["artifacts"]["full_background_page_manifest"]
                expected_page_count = page_manifest_record["metadata"]["page_count"]
                if extraction["available"] and extraction["page_count"] != expected_page_count:
                    extraction["available"] = False
                    extraction["error_code"] = "checkpoint_metadata_mismatch"
                    extraction.pop("payload", None)
                    _issue(
                        issues, "checkpoint_metadata_mismatch", setting=setting,
                        artifact="full_background_subtraction_pdf",
                        detail="source PDF page count does not match renderer page manifest",
                    )
                    setting_manifest["errors"].append("checkpoint_metadata_mismatch")
                if extraction["available"]:
                    slim_payload = extraction["payload"]
                    archive.writestr(slim_archive, slim_payload)
                    full_background.update({
                        "original_page_count": extraction["page_count"],
                        "extracted_pages": extraction["pages"],
                    })
                    meeting_summary.update({
                        "status": "exists",
                        "byte_size": len(slim_payload),
                        "sha256": sha256_bytes(slim_payload),
                        "backend": extraction["backend"],
                        "backend_identity": extraction["backend_identity"],
                    })
                else:
                    full_background["extraction_error"] = extraction["error_code"]
                    if extraction["error_code"] != "checkpoint_metadata_mismatch":
                        _issue(
                            issues, extraction["error_code"], setting=setting,
                            artifact="full_background_subtraction_pdf",
                        )
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
    parser.add_argument("--phi", help="optional profile setting override")
    parser.add_argument("--epsilon", help="optional profile epsilon override")
    parser.add_argument(
        "--profile",
        help="optional JSON bundle profile; defaults to the shipped profile",
    )
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
            profile_path=arguments.profile,
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
