#!/usr/bin/env python3
from __future__ import annotations

import sys

try:
    from ltsep_paths import emit_path_field_csv
except ModuleNotFoundError:
    from farm_env.ltsep_paths import emit_path_field_csv


def main() -> int:
    caller_path = sys.argv[1] if len(sys.argv) > 1 else None
    print(emit_path_field_csv(caller_path), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
