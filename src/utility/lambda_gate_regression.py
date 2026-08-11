"""Safe shared Lambda-gate regression-table updates.

The table lives in the common analysis output directory, so independent
epsilon/phi jobs can update it concurrently.  Keep the locking and atomic
replacement mechanics outside ROOT-dependent analysis modules for direct
non-ROOT regression coverage.
"""

from __future__ import annotations

import csv
import os
import uuid
from contextlib import contextmanager

try:
    import fcntl as _fcntl
except ImportError:  # Windows test environments have no POSIX advisory lock.
    _fcntl = None


def write_csv_atomically(path, fieldnames, rows):
    """Write a complete CSV to a unique sibling and atomically publish it."""
    temporary_path = "{}.{}.{}.tmp".format(path, os.getpid(), uuid.uuid4().hex)
    try:
        with open(temporary_path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({field: (row or {}).get(field) for field in fieldnames})
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)


@contextmanager
def advisory_linux_file_lock(path):
    """Take a cooperative lock for a shared Linux-output CSV update."""
    lock_path = "{}.lock".format(path)
    with open(lock_path, "a+", encoding="utf-8") as handle:
        if _fcntl is not None:
            _fcntl.flock(handle.fileno(), _fcntl.LOCK_EX)
        try:
            yield
        finally:
            if _fcntl is not None:
                _fcntl.flock(handle.fileno(), _fcntl.LOCK_UN)


def upsert_sorted_csv(path, row, fieldnames, key_fields):
    """Read, replace one setting row, then atomically publish sorted rows."""
    with advisory_linux_file_lock(path):
        rows_by_key = {}
        if os.path.exists(path):
            with open(path, newline="", encoding="utf-8") as handle:
                for existing in csv.DictReader(handle):
                    key = tuple(str(existing.get(field) or "") for field in key_fields)
                    rows_by_key[key] = {
                        field: existing.get(field) for field in fieldnames
                    }
        key = tuple(str((row or {}).get(field) or "") for field in key_fields)
        rows_by_key[key] = {
            field: (row or {}).get(field) for field in fieldnames
        }
        write_csv_atomically(
            path, fieldnames, [rows_by_key[key] for key in sorted(rows_by_key)]
        )
    return path
