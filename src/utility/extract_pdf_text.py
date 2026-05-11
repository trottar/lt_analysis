#! /usr/bin/python

"""
Best-effort PDF text extractor for local analysis notes and papers.

This script is intentionally dependency-light so it can run in the current
lt_analysis environment even when common PDF libraries are unavailable.

Extraction strategy:
1. Try `pypdf` / `PyPDF2` if installed.
2. Try a lightweight built-in parser for common text-based PDF streams.
3. Fall back to printable-string scraping from the raw PDF bytes.

Notes
-----
- This is NOT OCR. Image-only / scanned PDFs will still need OCR or manual
  screenshots.
- The built-in parser is best-effort. It works best on text-based PDFs with
  standard text operators and Flate-compressed streams.

Examples
--------
python src/utility/extract_pdf_text.py docs/2009_Williams_Qfit.pdf
python src/utility/extract_pdf_text.py docs --output-dir docs/pdf_text
python src/utility/extract_pdf_text.py docs/2009_Williams_Qfit.pdf --stdout
python src/utility/extract_pdf_text.py docs --cache-tar docs/pdf_text.tar.gz
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import tarfile
import time
import zlib
from io import BytesIO
from pathlib import Path
from typing import Iterable


PDF_SUFFIX = ".pdf"
TEXT_SUFFIX = ".txt"

STREAM_RE = re.compile(rb"stream\r?\n(.*?)\r?\nendstream", re.S)
BT_ET_RE = re.compile(rb"BT(.*?)ET", re.S)
TEXT_OP_RE = re.compile(
    rb"\(([^()]*(?:\\.[^()]*)*)\)\s*Tj"
    rb"|"
    rb"\[(.*?)\]\s*TJ"
    rb"|"
    rb"\(([^()]*(?:\\.[^()]*)*)\)\s*'"
    rb"|"
    rb"\(([^()]*(?:\\.[^()]*)*)\)\s*\"",
    re.S,
)
PAREN_RE = re.compile(rb"\(([^()]*(?:\\.[^()]*)*)\)")
PRINTABLE_STR_RE = re.compile(rb"[\x20-\x7e]{6,}")


def _configure_stdout() -> None:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass


def _iter_pdf_paths(paths: Iterable[str]) -> list[Path]:
    pdfs: list[Path] = []
    seen: set[Path] = set()
    for raw_path in paths:
        path = Path(raw_path)
        if path.is_dir():
            matches = sorted(path.rglob(f"*{PDF_SUFFIX}"))
        else:
            matches = [path]
        for match in matches:
            if match.suffix.lower() != PDF_SUFFIX:
                continue
            resolved = match.resolve()
            if resolved not in seen:
                seen.add(resolved)
                pdfs.append(match)
    return pdfs


def _unescape_pdf_string(payload: bytes) -> str:
    out = bytearray()
    i = 0
    while i < len(payload):
        byte = payload[i]
        if byte == 0x5C and i + 1 < len(payload):  # backslash
            i += 1
            esc = payload[i]
            mapping = {
                ord("n"): b"\n",
                ord("r"): b"\r",
                ord("t"): b"\t",
                ord("b"): b"\b",
                ord("f"): b"\f",
                ord("("): b"(",
                ord(")"): b")",
                ord("\\"): b"\\",
            }
            if esc in mapping:
                out.extend(mapping[esc])
            elif 48 <= esc <= 55:
                oct_digits = bytes([esc])
                for _ in range(2):
                    if i + 1 < len(payload) and 48 <= payload[i + 1] <= 55:
                        i += 1
                        oct_digits += bytes([payload[i]])
                    else:
                        break
                out.append(int(oct_digits, 8))
            else:
                out.append(esc)
        else:
            out.append(byte)
        i += 1
    return out.decode("latin1", errors="ignore")


def _normalize_text(text: str) -> str:
    lines = []
    for line in text.splitlines():
        compact = re.sub(r"\s+", " ", line).strip()
        if compact:
            lines.append(compact)
    return "\n".join(lines).strip()


def _extract_with_pypdf(pdf_path: Path) -> tuple[str, str] | None:
    for module_name in ("pypdf", "PyPDF2"):
        try:
            module = __import__(module_name)
        except Exception:
            continue
        try:
            reader = module.PdfReader(str(pdf_path))
            page_text = []
            for page in reader.pages:
                try:
                    text = page.extract_text() or ""
                except Exception:
                    text = ""
                if text.strip():
                    page_text.append(text)
            joined = _normalize_text("\n\n".join(page_text))
            if joined:
                return joined, module_name
        except Exception:
            continue
    return None


def _extract_tj_text(block: bytes) -> list[str]:
    fragments: list[str] = []
    for match in TEXT_OP_RE.finditer(block):
        if match.group(1):
            fragments.append(_unescape_pdf_string(match.group(1)))
        elif match.group(2):
            pieces = [_unescape_pdf_string(pm.group(1)) for pm in PAREN_RE.finditer(match.group(2))]
            if pieces:
                fragments.append("".join(pieces))
        elif match.group(3):
            fragments.append(_unescape_pdf_string(match.group(3)))
        elif match.group(4):
            fragments.append(_unescape_pdf_string(match.group(4)))
    return fragments


def _decode_candidate_stream(raw_stream: bytes) -> list[bytes]:
    candidates = [raw_stream]
    try:
        candidates.append(zlib.decompress(raw_stream))
    except Exception:
        pass
    return candidates


def _extract_with_builtin_parser(pdf_path: Path) -> tuple[str, str] | None:
    data = pdf_path.read_bytes()
    blocks: list[str] = []

    for stream_match in STREAM_RE.finditer(data):
        raw_stream = stream_match.group(1)
        for candidate in _decode_candidate_stream(raw_stream):
            if b"BT" not in candidate and b"Tj" not in candidate and b"TJ" not in candidate:
                continue

            bt_blocks = BT_ET_RE.findall(candidate)
            if not bt_blocks:
                bt_blocks = [candidate]

            for block in bt_blocks:
                parts = _extract_tj_text(block)
                if parts:
                    joined = _normalize_text(" ".join(parts))
                    if joined:
                        blocks.append(joined)

    if not blocks:
        return None

    deduped: list[str] = []
    seen: set[str] = set()
    for block in blocks:
        if block not in seen:
            seen.add(block)
            deduped.append(block)
    return "\n".join(deduped), "builtin-stream-parser"


def _extract_with_strings(pdf_path: Path, min_len: int = 6) -> tuple[str, str] | None:
    data = pdf_path.read_bytes()
    pattern = re.compile(rb"[\x20-\x7e]{" + str(max(1, int(min_len))).encode("ascii") + rb",}")

    chunks = [match.group(0).decode("latin1", errors="ignore") for match in pattern.finditer(data)]
    cleaned = []
    for chunk in chunks:
        compact = _normalize_text(chunk)
        if compact and compact not in cleaned:
            cleaned.append(compact)

    if not cleaned:
        return None
    return "\n".join(cleaned), "raw-string-scrape"


def extract_pdf_text(pdf_path: Path, min_string_len: int = 6) -> tuple[str, str]:
    if not pdf_path.exists():
        raise FileNotFoundError(f"PDF not found: {pdf_path}")

    for extractor in (
        _extract_with_pypdf,
        _extract_with_builtin_parser,
        lambda p: _extract_with_strings(p, min_string_len),
    ):
        result = extractor(pdf_path)
        if result and result[0].strip():
            return result
    return "", "no-text-found"


def build_report(pdf_path: Path, text: str, method: str) -> str:
    header = [
        f"Source PDF: {pdf_path.resolve()}",
        f"Extraction method: {method}",
        "Warning: best-effort extraction only; scanned/image PDFs may require OCR.",
        "",
    ]
    if not text.strip():
        header.append("No text could be extracted from this PDF.")
        return "\n".join(header)
    return "\n".join(header) + text + "\n"


def _common_parent(paths: list[Path]) -> Path:
    if not paths:
        raise ValueError("At least one path is required")
    common = os.path.commonpath([str(path.resolve().parent) for path in paths])
    return Path(common)


def _default_output_dir(pdf_paths: list[Path]) -> Path:
    return _common_parent(pdf_paths) / "pdf_text"


def _default_cache_tar(output_dir: Path) -> Path:
    return output_dir.with_suffix(".tar.gz")


def _member_name(pdf_path: Path, common_root: Path) -> str:
    relative = pdf_path.resolve().relative_to(common_root.resolve())
    return relative.with_suffix(TEXT_SUFFIX).as_posix()


def _output_path_for_pdf(pdf_path: Path, output_dir: Path, common_root: Path) -> Path:
    return output_dir / _member_name(pdf_path, common_root)


def _load_tar_cache(cache_tar: Path) -> dict[str, tuple[str, float]]:
    cached_reports: dict[str, tuple[str, float]] = {}
    if not cache_tar.exists():
        return cached_reports

    with tarfile.open(cache_tar, "r:*") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            extracted = archive.extractfile(member)
            if extracted is None:
                continue
            payload = extracted.read().decode("utf-8", errors="replace")
            cached_reports[member.name] = (payload, float(member.mtime))
    return cached_reports


def _write_tar_cache(cache_tar: Path, cached_reports: dict[str, tuple[str, float]]) -> None:
    cache_tar.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(cache_tar, "w:gz") as archive:
        for member_name in sorted(cached_reports):
            report, member_mtime = cached_reports[member_name]
            payload = report.encode("utf-8")
            info = tarfile.TarInfo(name=member_name)
            info.size = len(payload)
            info.mtime = member_mtime
            archive.addfile(info, BytesIO(payload))


def _is_cache_valid(member_mtime: float, pdf_path: Path) -> bool:
    try:
        pdf_mtime = pdf_path.stat().st_mtime
    except FileNotFoundError:
        return False
    return member_mtime >= pdf_mtime


def _write_report(out_path: Path, report: str) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report, encoding="utf-8")


def main() -> int:
    _configure_stdout()

    parser = argparse.ArgumentParser(description="Best-effort PDF-to-text extractor for lt_analysis docs.")
    parser.add_argument("paths", nargs="+", help="PDF files or directories to process")
    parser.add_argument("-o", "--output", help="Exact output .txt path (single PDF only)")
    parser.add_argument("--output-dir", help="Directory for generated text files")
    parser.add_argument("--cache-tar", help="Tarball cache for extracted text reports")
    parser.add_argument("--no-cache", action="store_true", help="Disable tar cache reads/writes")
    parser.add_argument("--stdout", action="store_true", help="Write extracted text to stdout instead of files")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing text files")
    parser.add_argument("--min-string-len", type=int, default=6, help="Minimum printable-string length for raw fallback mode")
    args = parser.parse_args()

    pdf_paths = _iter_pdf_paths(args.paths)
    if not pdf_paths:
        print("No PDF files found.")
        return 1

    if args.output and len(pdf_paths) != 1:
        parser.error("--output can only be used with a single PDF input")

    if args.output:
        output_dir = Path(args.output).parent
        common_root = pdf_paths[0].resolve().parent
    else:
        output_dir = Path(args.output_dir) if args.output_dir else _default_output_dir(pdf_paths)
        common_root = _common_parent(pdf_paths)
        output_dir.mkdir(parents=True, exist_ok=True)

    cache_tar = None if args.no_cache else Path(args.cache_tar) if args.cache_tar else _default_cache_tar(output_dir)
    cached_reports = {} if cache_tar is None else _load_tar_cache(cache_tar)
    cache_dirty = False
    failures = 0
    for pdf_path in pdf_paths:
        out_path = Path(args.output) if args.output else _output_path_for_pdf(pdf_path, output_dir, common_root)
        member_name = out_path.name if args.output else out_path.relative_to(output_dir).as_posix()
        body = ""

        cached_payload = cached_reports.get(member_name)
        if not args.overwrite and cached_payload and _is_cache_valid(cached_payload[1], pdf_path):
            report, _ = cached_payload
            method = "cache-tar"
            parts = report.split("\n\n", 1)
            body = parts[1] if len(parts) == 2 else report
        else:
            text, method = extract_pdf_text(pdf_path, min_string_len=args.min_string_len)
            report = build_report(pdf_path, text, method)
            body = text
            if cache_tar is not None:
                cached_reports[member_name] = (
                    report,
                    float(max(int(time.time()), int(pdf_path.stat().st_mtime))),
                )
                cache_dirty = True

        if args.stdout:
            print("=" * 100)
            print(report)
            continue

        if out_path.exists() and not args.overwrite:
            print(f"Skipping existing file: {out_path}")
            continue

        _write_report(out_path, report)
        status = "ok" if body.strip() else "empty"
        print(f"[{status}] {pdf_path} -> {out_path} ({method})")
        if not body.strip():
            failures += 1

    if cache_tar is not None and cache_dirty:
        _write_tar_cache(cache_tar, cached_reports)
        print(f"[cache] wrote tarball cache: {cache_tar}")

    return 0 if failures == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
