"""Pure-Python coverage for detached Phase-D checkpoint persistence."""

from __future__ import annotations

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
path = str(REPO_ROOT / "src" / "cuts")
if path not in sys.path:
    sys.path.insert(0, path)


import pion_hgcer_phase_d_checkpoint as phase_d


SOURCE_FINGERPRINT = "phase-c-source-fingerprint"


def _setting(Q2="4p4", W="2p74"):
    return {
        "kinematic_token": "Q4p4W2p74",
        "Q2": Q2,
        "W": W,
        "epsilon_setting": "low",
        "epsilon_filename_token": "lowe",
        "phi_setting": "Left",
        "particle_type": "kaon",
    }


def _representation(name):
    return {
        "schema_version": "synthetic-{}".format(name),
        "status": "available",
        "available": True,
        "source_checkpoint_payload_fingerprint": SOURCE_FINGERPRINT,
        "cells": [{"t_index": 0, "diagnostic": name}],
    }


def _comparison(status="available"):
    available = status == "available"
    result = {
        "schema_version": "pion_hgcer_ab_comparison/v1",
        "status": status,
        "available": available,
        "reason": None if available else "method_b_comparison_unavailable",
        "non_authoritative": True,
        "comparison_performed": available,
        "classification_performed": available,
        "classification_scope": "availability_only_non_prescriptive" if available else "none",
        "decision_performed": False,
        "statistical_compatibility_claimed": False,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }
    if available:
        result["source_checkpoint_payload_fingerprint"] = SOURCE_FINGERPRINT
        result["cells"] = [{"comparison": {"availability": "both_comparable"}}]
    return result


class PionHGCerPhaseDCheckpointTests(unittest.TestCase):
    def _payload(self, *, status="available", Q2="4p4", W="2p74"):
        return phase_d.build_pion_hgcer_phase_d_checkpoint(
            setting=_setting(Q2, W),
            method_a_comparison=_representation("a"),
            method_b_comparison=_representation("b"),
            ab_comparison=_comparison(status),
        )

    def test_available_payload_is_complete_detached_and_preserves_setting_scalars(self):
        payload = self._payload()
        self.assertEqual(payload["schema_version"], phase_d.PHASE_D_CHECKPOINT_SCHEMA_VERSION)
        self.assertTrue(payload["available"])
        self.assertEqual(payload["source_checkpoint_payload_fingerprint"], SOURCE_FINGERPRINT)
        self.assertEqual(payload["setting"]["Q2"], "4p4")
        self.assertEqual(payload["setting"]["W"], "2p74")
        self.assertIn("method_a_comparison", payload)
        self.assertIn("method_b_comparison", payload)
        self.assertIn("ab_comparison", payload)
        payload["method_a_comparison"]["cells"][0]["diagnostic"] = "changed"
        self.assertEqual(_representation("a")["cells"][0]["diagnostic"], "a")
        numeric = self._payload(Q2=4.4, W=2.74)
        self.assertEqual(numeric["setting"]["Q2"], 4.4)
        self.assertEqual(numeric["setting"]["W"], 2.74)

    def test_unavailable_comparison_still_builds_a_review_artifact(self):
        payload = self._payload(status="unavailable")
        self.assertEqual(payload["status"], "unavailable")
        self.assertFalse(payload["available"])
        self.assertEqual(payload["reason"], "method_b_comparison_unavailable")
        self.assertFalse(payload["comparison_performed"])
        self.assertFalse(payload["classification_performed"])
        self.assertEqual(payload["classification_scope"], "none")
        self.assertFalse(payload["decision_performed"])

    def test_filename_validation_is_deterministic(self):
        self.assertEqual(
            phase_d.pion_hgcer_phase_d_checkpoint_filename(
                "Left", "kaon", "Q4p4W2p74", "lowe"
            ),
            "Left_kaon_pion-background_hgcer_ab_comparison_Q4p4W2p74_lowe.json",
        )
        for args in (
            ("Left", "pion", "Q4p4W2p74", "lowe"),
            ("Left", "kaon", "Q4p4W2p74", "low"),
            ("Left/unsafe", "kaon", "Q4p4W2p74", "lowe"),
        ):
            with self.subTest(args=args):
                with self.assertRaises(ValueError):
                    phase_d.pion_hgcer_phase_d_checkpoint_filename(*args)

    def test_json_safety_forbidden_fields_and_writer_are_deterministic(self):
        payload = self._payload()
        before = copy.deepcopy(payload)
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary) / "first.json"
            second = Path(temporary) / "second.json"
            phase_d.write_pion_hgcer_phase_d_checkpoint_json(first, payload)
            phase_d.write_pion_hgcer_phase_d_checkpoint_json(second, payload)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertTrue(first.read_bytes().endswith(b"\n"))
            self.assertEqual(json.loads(first.read_text(encoding="utf-8")), payload)
        self.assertEqual(payload, before)
        bad = _comparison()
        bad["C_final"] = 1.0
        with self.assertRaisesRegex(ValueError, "forbidden"):
            phase_d.build_pion_hgcer_phase_d_checkpoint(
                setting=_setting(),
                method_a_comparison=_representation("a"),
                method_b_comparison=_representation("b"),
                ab_comparison=bad,
            )
        bad = _representation("a")
        bad["value"] = float("nan")
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            phase_d.build_pion_hgcer_phase_d_checkpoint(
                setting=_setting(),
                method_a_comparison=bad,
                method_b_comparison=_representation("b"),
                ab_comparison=_comparison(),
            )
        with self.assertRaisesRegex(ValueError, "json_safe"):
            phase_d.write_pion_hgcer_phase_d_checkpoint_json("ignored.json", {"bad": object()})


if __name__ == "__main__":
    unittest.main()
