"""Authoritative kaon-derived experimental coordinate contracts.

The analysis must not infer its coordinate system from optional derived ROOT
branches.  This module keeps the resolved kaon MM and |t| corrections as a
small, serializable contract which every experimental source can share.
"""

from __future__ import annotations

import hashlib
import json
import math


COORDINATE_SCHEMA_VERSION = "kaon_data_coordinate/v1"


def _json_ready(value):
    """Return only deterministic JSON primitives for a shift-provenance hash."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise RuntimeError("kaon_data_coordinate_nonfinite_provenance")
        return float(value)
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise RuntimeError(
            "kaon_data_coordinate_unserializable_provenance:{}".format(type(value).__name__)
        ) from exc
    if not math.isfinite(numeric):
        raise RuntimeError("kaon_data_coordinate_nonfinite_provenance")
    return numeric


def _resolved_shift(entry, coordinate_name, *, required=True):
    if not isinstance(entry, dict):
        if required:
            raise RuntimeError("kaon_data_coordinate_missing_{}_shift".format(coordinate_name))
        return 0.0
    try:
        # Shift-prep entries use ``{"shift": ...}``, while the immutable
        # contract stores the resolved values as ``mm_shift`` / ``t_shift``.
        contract_key = "{}_shift".format(coordinate_name)
        shift = float(
            entry[contract_key] if contract_key in entry else entry["shift"]
        )
    except (KeyError, TypeError, ValueError) as exc:
        if required:
            raise RuntimeError("kaon_data_coordinate_missing_{}_shift".format(coordinate_name)) from exc
        return 0.0
    if not math.isfinite(shift):
        raise RuntimeError("kaon_data_coordinate_nonfinite_{}_shift".format(coordinate_name))
    return shift


def coordinate_fingerprint(contract):
    """Fingerprint the physical transform and its shift-fit provenance."""
    if not isinstance(contract, dict):
        raise RuntimeError("kaon_data_coordinate_contract_missing")
    payload = {
        key: value
        for key, value in contract.items()
        if key not in ("coordinate_fingerprint", "fingerprint")
    }
    encoded = json.dumps(
        _json_ready(payload), sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def build_kaon_data_coordinate_contract(
    phi_setting,
    mm_shift_entry,
    t_shift_entry,
    *,
    require_t_shift=True,
):
    """Build the immutable-by-convention experimental coordinate contract."""
    mm_shift = _resolved_shift(mm_shift_entry, "mm", required=True)
    t_shift = _resolved_shift(t_shift_entry, "t", required=require_t_shift)
    contract = {
        "schema_version": COORDINATE_SCHEMA_VERSION,
        "coordinate_system": "kaon_derived_experimental_analysis",
        "phi_setting": str(phi_setting),
        "source_particle": "kaon",
        "raw_mm_field": "MM",
        "raw_t_expression": "-MandelT",
        "mm_shift": float(mm_shift),
        "t_shift": float(t_shift),
        "mm_shift_provenance": _json_ready(dict(mm_shift_entry or {})),
        "t_shift_provenance": _json_ready(dict(t_shift_entry or {})),
    }
    contract["coordinate_fingerprint"] = coordinate_fingerprint(contract)
    return contract


def validate_kaon_data_coordinate_contract(contract, *, phi_setting=None, require_t_shift=True):
    """Validate a contract received across an analysis ownership boundary."""
    if not isinstance(contract, dict):
        raise RuntimeError("kaon_data_coordinate_contract_missing")
    if contract.get("schema_version") != COORDINATE_SCHEMA_VERSION:
        raise RuntimeError("kaon_data_coordinate_schema_mismatch")
    if contract.get("coordinate_system") != "kaon_derived_experimental_analysis":
        raise RuntimeError("kaon_data_coordinate_system_mismatch")
    if contract.get("source_particle") != "kaon":
        raise RuntimeError("kaon_data_coordinate_source_particle_mismatch")
    if phi_setting is not None and str(contract.get("phi_setting")) != str(phi_setting):
        raise RuntimeError("kaon_data_coordinate_phi_setting_mismatch")
    _resolved_shift(contract, "mm", required=True)
    _resolved_shift(contract, "t", required=require_t_shift)
    declared = str(contract.get("coordinate_fingerprint") or "")
    if not declared or declared != coordinate_fingerprint(contract):
        raise RuntimeError("kaon_data_coordinate_fingerprint_mismatch")
    return contract


def raw_event_coordinates(evt):
    """Read the raw experimental convention without consulting derived branches."""
    return float(getattr(evt, "MM")), -float(getattr(evt, "MandelT"))


def analysis_event_coordinates(evt, contract):
    """Return raw and kaon-analysis coordinates for one experimental record."""
    # The contract is validated once at each ownership boundary.  This helper
    # runs inside ROOT event loops, so do not rehash provenance per record.
    if not isinstance(contract, dict):
        raise RuntimeError("kaon_data_coordinate_contract_missing")
    raw_mm, raw_t = raw_event_coordinates(evt)
    return {
        "raw_mm": raw_mm,
        "raw_t": raw_t,
        "analysis_mm": raw_mm + float(contract["mm_shift"]),
        "analysis_t": raw_t + float(contract["t_shift"]),
    }
