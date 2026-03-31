#!/usr/bin/env python3
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


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


def resolve_ltsep_paths(caller_path: Optional[str] = None) -> LtsepPaths:
    from ltsep import Root

    probe = os.path.realpath(caller_path or __file__)
    lt = Root(probe, "Plot_LTSep")
    return LtsepPaths(
        rootpath=_normalize_base(lt.ROOTPATH),
        skimpath=_normalize_base(lt.SKIMPATH),
        anatype=str(lt.ANATYPE),
        user=str(lt.USER),
        host=str(lt.HOST),
    )
