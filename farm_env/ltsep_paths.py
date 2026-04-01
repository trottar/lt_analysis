#!/usr/bin/env python3
from __future__ import annotations

import getpass
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

try:
    import pwd
except ImportError:  # pragma: no cover - unavailable on Windows
    pwd = None  # type: ignore[assignment]


PATH_FIELD_NAMES = (
    "VOLATILEPATH",
    "ANALYSISPATH",
    "HCANAPATH",
    "REPLAYPATH",
    "UTILPATH",
    "PACKAGEPATH",
    "OUTPATH",
    "ROOTPATH",
    "SKIMPATH",
    "REPORTPATH",
    "CUTPATH",
    "PARAMPATH",
    "SCRIPTPATH",
    "ANATYPE",
    "USER",
    "HOST",
    "SIMCPATH",
    "LTANAPATH",
)


@dataclass(frozen=True)
class LtsepPaths:
    rootpath: Path
    skimpath: Path
    anatype: str
    user: str
    host: str

    @property
    def replay_source_dir(self) -> Path:
        return _analysis_dir(self.rootpath, self.anatype)

    @property
    def skim_source_dir(self) -> Path:
        return _analysis_dir(self.skimpath, self.anatype)


def _normalize_base(path_value: object) -> Path:
    text = os.path.expandvars(str(path_value)).strip()
    return Path(text).expanduser()


def _analysis_dir(base_path: Path, anatype: str) -> Path:
    expected_leaf = f"{anatype}LT"
    path_text = str(base_path)
    if "None" in path_text:
        return Path(path_text.replace("None", expected_leaf)).expanduser()
    if base_path.name == expected_leaf:
        return base_path
    return base_path / expected_leaf


def _safe_login_name() -> str:
    for env_name in ("LOGNAME", "USER", "LNAME", "USERNAME"):
        value = os.environ.get(env_name)
        if value:
            return value
    if pwd is not None:
        try:
            return pwd.getpwuid(os.getuid()).pw_name
        except Exception:
            pass
    return getpass.getuser()


def _ensure_safe_getlogin() -> None:
    try:
        os.getlogin()
    except Exception:
        os.getlogin = _safe_login_name  # type: ignore[assignment]


def resolve_ltsep_root(caller_path: Optional[str] = None):
    _ensure_safe_getlogin()
    from ltsep import Root

    probe = os.path.realpath(caller_path or __file__)
    return Root(probe, "Plot_LTSep")


def resolve_ltsep_paths(caller_path: Optional[str] = None) -> LtsepPaths:
    lt = resolve_ltsep_root(caller_path)
    return LtsepPaths(
        rootpath=_normalize_base(lt.ROOTPATH),
        skimpath=_normalize_base(lt.SKIMPATH),
        anatype=str(lt.ANATYPE),
        user=str(lt.USER),
        host=str(lt.HOST),
    )


def emit_path_field_csv(caller_path: Optional[str] = None) -> str:
    lt = resolve_ltsep_root(caller_path)
    return ",".join(str(getattr(lt, field_name, "")) for field_name in PATH_FIELD_NAMES)
