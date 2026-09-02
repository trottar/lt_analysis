"""JSON-only persistence for detached non-authoritative Phase-D diagnostics."""

from __future__ import annotations

from collections.abc import Mapping
import json
import math
import os
from types import MappingProxyType


PHASE_D_CHECKPOINT_SCHEMA_VERSION = "pion_hgcer_phase_d_checkpoint/v1"

_FORBIDDEN_FIELD_NAMES = {
    "C_A",
    "C_B",
    "C",
    "C_final",
    "combined_correction",
    "final_correction",
    "refined_pion_weight",
    "applied_refinement_weight",
    "w_pi_refined",
}


def _json_safe(value, context="payload"):
    """Detach JSON primitives only and reject nonfinite or corrective content."""
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


def _filename_token(value):
    if isinstance(value, (bool, list, tuple, dict, set, MappingProxyType)) or isinstance(
        value, Mapping
    ):
        raise ValueError("phase_d_checkpoint_filename_token_invalid")
    if not isinstance(value, (str, int, float)):
        raise ValueError("phase_d_checkpoint_filename_token_invalid")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("phase_d_checkpoint_filename_token_invalid")
    raw = str(value)
    token = raw.strip()
    if (
        not token
        or token != raw
        or any(character.isspace() for character in token)
        or any(character in token for character in "\\\\/:")
        or ".." in token
    ):
        raise ValueError("phase_d_checkpoint_filename_token_invalid")
    return token


def _setting(setting):
    if not isinstance(setting, Mapping):
        raise ValueError("phase_d_checkpoint_setting_invalid")
    required = (
        "kinematic_token",
        "Q2",
        "W",
        "epsilon_setting",
        "epsilon_filename_token",
        "phi_setting",
        "particle_type",
    )
    if any(key not in setting for key in required):
        raise ValueError("phase_d_checkpoint_setting_invalid")
    result = {
        "kinematic_token": _filename_token(setting["kinematic_token"]),
        "Q2": _json_safe(setting["Q2"], "setting.Q2"),
        "W": _json_safe(setting["W"], "setting.W"),
        "epsilon_setting": _filename_token(setting["epsilon_setting"]).lower(),
        "epsilon_filename_token": _filename_token(
            setting["epsilon_filename_token"]
        ).lower(),
        "phi_setting": _filename_token(setting["phi_setting"]),
        "particle_type": _filename_token(setting["particle_type"]).lower(),
    }
    if result["particle_type"] != "kaon":
        raise ValueError("phase_d_checkpoint_particle_type_invalid")
    if result["epsilon_setting"] not in {"high", "low"}:
        raise ValueError("phase_d_checkpoint_epsilon_setting_invalid")
    if result["epsilon_filename_token"] != "{}e".format(
        result["epsilon_setting"]
    ):
        raise ValueError("phase_d_checkpoint_epsilon_filename_token_invalid")
    return result


def pion_hgcer_phase_d_checkpoint_filename(
    phi_setting, particle_type, kinematic_token, epsilon_filename_token
):
    """Return the deterministic setting-level Phase-D diagnostic basename."""
    phi = _filename_token(phi_setting)
    particle = _filename_token(particle_type).lower()
    kinematic = _filename_token(kinematic_token)
    epsilon = _filename_token(epsilon_filename_token).lower()
    if particle != "kaon":
        raise ValueError("phase_d_checkpoint_particle_type_invalid")
    if epsilon not in {"highe", "lowe"}:
        raise ValueError("phase_d_checkpoint_epsilon_filename_token_invalid")
    return "{}_{}_pion-background_hgcer_ab_comparison_{}_{}.json".format(
        phi, particle, kinematic, epsilon
    )


def _source_fingerprint(method_a_comparison, method_b_comparison, ab_comparison):
    candidates = [
        value.get("source_checkpoint_payload_fingerprint")
        for value in (ab_comparison, method_a_comparison, method_b_comparison)
        if isinstance(value, Mapping)
        and isinstance(value.get("source_checkpoint_payload_fingerprint"), str)
        and value.get("source_checkpoint_payload_fingerprint").strip()
    ]
    if not candidates:
        return None
    if any(value != candidates[0] for value in candidates):
        raise ValueError("phase_d_checkpoint_source_fingerprint_invalid")
    return candidates[0]


def build_pion_hgcer_phase_d_checkpoint(
    *, setting, method_a_comparison, method_b_comparison, ab_comparison
):
    """Build a self-contained detached D.2/D.3/D.4 review artifact."""
    if not isinstance(method_a_comparison, Mapping) or not isinstance(
        method_b_comparison, Mapping
    ) or not isinstance(ab_comparison, Mapping):
        raise ValueError("phase_d_checkpoint_comparison_invalid")
    resolved_setting = _setting(setting)
    method_a = _json_safe(method_a_comparison, "method_a_comparison")
    method_b = _json_safe(method_b_comparison, "method_b_comparison")
    comparison = _json_safe(ab_comparison, "ab_comparison")
    status = comparison.get("status")
    available = comparison.get("available")
    if not isinstance(status, str) or not isinstance(available, bool):
        raise ValueError("phase_d_checkpoint_comparison_invalid")
    required_flags = (
        "non_authoritative",
        "comparison_performed",
        "classification_performed",
        "classification_scope",
        "decision_performed",
        "statistical_compatibility_claimed",
        "production_objects_mutated",
        "refinement_applied",
    )
    if any(key not in comparison for key in required_flags):
        raise ValueError("phase_d_checkpoint_comparison_invalid")
    if (
        comparison["non_authoritative"] is not True
        or comparison["decision_performed"] is not False
        or comparison["statistical_compatibility_claimed"] is not False
        or comparison["production_objects_mutated"] is not False
        or comparison["refinement_applied"] is not False
    ):
        raise ValueError("phase_d_checkpoint_authority_invalid")
    payload = {
        "schema_version": PHASE_D_CHECKPOINT_SCHEMA_VERSION,
        "setting": resolved_setting,
        "source_checkpoint_payload_fingerprint": _source_fingerprint(
            method_a, method_b, comparison
        ),
        "method_a_comparison": method_a,
        "method_b_comparison": method_b,
        "ab_comparison": comparison,
        "status": status,
        "available": available,
        "reason": comparison.get("reason"),
        "non_authoritative": True,
        "comparison_performed": comparison["comparison_performed"],
        "classification_performed": comparison["classification_performed"],
        "classification_scope": comparison["classification_scope"],
        "decision_performed": False,
        "statistical_compatibility_claimed": False,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }
    return _json_safe(payload)


def write_pion_hgcer_phase_d_checkpoint_json(path, payload):
    """Write a deterministic, detached Phase-D JSON artifact."""
    target = os.fspath(path)
    serialized = _json_safe(payload)
    with open(target, "w", encoding="utf-8") as handle:
        json.dump(serialized, handle, sort_keys=True, indent=2, allow_nan=False)
        handle.write("\n")
    return target


__all__ = (
    "PHASE_D_CHECKPOINT_SCHEMA_VERSION",
    "build_pion_hgcer_phase_d_checkpoint",
    "pion_hgcer_phase_d_checkpoint_filename",
    "write_pion_hgcer_phase_d_checkpoint_json",
)
