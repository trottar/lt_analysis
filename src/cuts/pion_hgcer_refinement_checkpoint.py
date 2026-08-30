"""Detached, persistent Phase-A/B/C HGCer checkpoint serialization."""

from __future__ import annotations

import json
import math
import os
from collections.abc import Mapping
from types import MappingProxyType


CHECKPOINT_SCHEMA_VERSION = "pion_hgcer_refinement_checkpoint/v1"

_FORBIDDEN_FIELD_NAMES = {
    "C_B",
    "C_final",
    "refined_pion_weight",
    "applied_refinement_weight",
}


def _json_safe(value, context="payload"):
    """Copy only JSON primitives/mappings/sequences and reject opaque objects."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("{}_contains_nonfinite_float".format(context))
        return value
    if isinstance(value, MappingProxyType):
        value = dict(value)
    if isinstance(value, Mapping):
        result = {}
        for key, child in value.items():
            string_key = str(key)
            if string_key in _FORBIDDEN_FIELD_NAMES:
                raise ValueError("{}_contains_forbidden_field".format(context))
            result[string_key] = _json_safe(child, "{}.{}".format(context, string_key))
        return result
    if isinstance(value, (list, tuple)):
        return [
            _json_safe(child, "{}[{}]".format(context, index))
            for index, child in enumerate(value)
        ]
    try:
        if hasattr(value, "item"):
            return _json_safe(value.item(), context)
        if hasattr(value, "tolist"):
            return _json_safe(value.tolist(), context)
    except Exception as exc:
        raise ValueError("{}_is_not_json_safe".format(context)) from exc
    raise ValueError("{}_is_not_json_safe".format(context))


def _summary(value):
    return _json_safe(value if isinstance(value, Mapping) else {})


def _method_payload(value, summary):
    method = value if isinstance(value, Mapping) else {}
    return {
        "summary": _summary(summary),
        "fingerprint": method.get("fingerprint"),
        "cells": method.get("cells") or [],
        "status": method.get("status", "unavailable"),
        "available": bool(method.get("available", False)),
        "reason": method.get("reason"),
    }


def pion_hgcer_refinement_checkpoint_filename(
    phi_setting, kinematic_token, epsilon_setting
):
    """Return the deterministic setting-level Phase-A/B/C checkpoint basename."""
    pieces = (phi_setting, kinematic_token, epsilon_setting)
    normalized = []
    for value in pieces:
        token = str(value or "").strip()
        if not token or any(character in token for character in "\\\\/:"):
            raise ValueError("checkpoint_filename_token_invalid")
        normalized.append(token)
    return "{}_pion_hgcer_refinement_checkpoint_{}_{}.json".format(
        *normalized
    )


def build_pion_hgcer_refinement_checkpoint(
    *,
    setting,
    phase_a,
    method_a,
    method_b,
    phase_a_summary=None,
    method_a_summary=None,
    method_b_summary=None,
):
    """Build a pure JSON-safe A+B+C checkpoint without comparing A to B."""
    setting = setting if isinstance(setting, Mapping) else {}
    phase = phase_a if isinstance(phase_a, Mapping) else {}
    method_b_payload = _method_payload(method_b, method_b_summary)
    method_b_payload["parent_region_references"] = (
        (method_b or {}).get("parent_region_references") or []
        if isinstance(method_b, Mapping)
        else []
    )
    payload = {
        "schema_version": CHECKPOINT_SCHEMA_VERSION,
        "setting": {
            "kinematic_token": setting.get("kinematic_token"),
            "Q2": setting.get("Q2"),
            "W": setting.get("W"),
            "epsilon_setting": setting.get("epsilon_setting"),
            "phi_setting": setting.get("phi_setting"),
        },
        "phase_a": {
            "summary": _summary(phase_a_summary),
            "contract_fingerprint": phase.get("contract_fingerprint"),
            "coordinate_fingerprint": phase.get("coordinate_fingerprint"),
            "host_state": phase.get("host_state"),
            "source_target_state": phase.get("source_target_state"),
        },
        "method_a": _method_payload(method_a, method_a_summary),
        "method_b": method_b_payload,
        "canonical_t_edges": phase.get("canonical_t_edges") or phase.get("t_edges") or [],
        "delta_edges": phase.get("delta_edges") or [],
        "host_state_summary": {
            "phase_a_host_state": phase.get("host_state"),
            "method_b_host_state": (
                method_b.get("host_state") if isinstance(method_b, Mapping) else None
            ),
            "source_target_state": phase.get("source_target_state"),
        },
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }
    result = _json_safe(payload)
    json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return result


def write_pion_hgcer_refinement_checkpoint_json(path, payload):
    """Serialize a supplied detached checkpoint and return its path unchanged."""
    target = os.fspath(path)
    serialized = _json_safe(payload)
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(
            serialized,
            handle,
            sort_keys=True,
            indent=2,
            allow_nan=False,
        )
        handle.write("\n")
    return target


__all__ = [
    "CHECKPOINT_SCHEMA_VERSION",
    "build_pion_hgcer_refinement_checkpoint",
    "pion_hgcer_refinement_checkpoint_filename",
    "write_pion_hgcer_refinement_checkpoint_json",
]
