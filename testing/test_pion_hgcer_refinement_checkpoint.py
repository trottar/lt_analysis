"""Pure-Python contracts for detached Phase-A/B/C checkpoint serialization."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[1]
for relative_path in ("src/cuts", "src/utility"):
    path = str(REPO_ROOT / relative_path)
    if path not in sys.path:
        sys.path.insert(0, path)

import pion_hgcer_refinement_checkpoint as checkpoint


def _phase(status="available"):
    return {
        "status": status,
        "available": status == "available",
        "contract_fingerprint": "phase-a-fingerprint",
        "coordinate_fingerprint": "coordinates-fingerprint",
        "host_state": "proton_cleaned",
        "source_target_state": "post_proton_noRF",
        "canonical_t_edges": [0.0, 1.0],
        "delta_edges": [0.0, 10.0],
    }


def _method(name, status="available"):
    return {
        "status": status,
        "available": status == "available",
        "reason": None if status == "available" else "synthetic_unavailable",
        "fingerprint": "{}-fingerprint".format(name),
        "host_state": "proton_cleaned",
        "cells": [{"t_index": 0, "delta_index": 0, "method": name}],
        "parent_region_references": ([{"t_index": 0, "region_name": "pi_n"}] if name == "b" else []),
    }


def _payload(method_a_status="available", method_b_status="available"):
    return checkpoint.build_pion_hgcer_refinement_checkpoint(
        setting={
            "kinematic_token": "Q4p4W2p74",
            "Q2": 4.4,
            "W": 2.74,
            "epsilon_setting": "lowe",
            "phi_setting": "Left",
        },
        phase_a=_phase(),
        phase_a_summary={"status": "available", "record_count": 12},
        method_a=_method("a", method_a_status),
        method_a_summary={"status": method_a_status},
        method_b=_method("b", method_b_status),
        method_b_summary={"status": method_b_status},
    )


class PionHGCerRefinementCheckpointTests(unittest.TestCase):
    def test_payload_retains_full_detached_a_b_c_tables_and_provenance(self):
        payload = _payload()
        self.assertEqual(payload["phase_a"]["contract_fingerprint"], "phase-a-fingerprint")
        self.assertEqual(payload["method_a"]["cells"], _method("a")["cells"])
        self.assertEqual(payload["method_b"]["cells"], _method("b")["cells"])
        self.assertEqual(payload["method_b"]["parent_region_references"], _method("b")["parent_region_references"])
        self.assertEqual(payload["host_state_summary"]["phase_a_host_state"], "proton_cleaned")
        self.assertTrue(payload["non_authoritative"])
        self.assertFalse(payload["production_objects_mutated"])
        self.assertFalse(payload["refinement_applied"])

    def test_writer_produces_deterministic_valid_json(self):
        payload = _payload()
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary) / "first.json"
            second = Path(temporary) / "second.json"
            checkpoint.write_pion_hgcer_refinement_checkpoint_json(first, payload)
            checkpoint.write_pion_hgcer_refinement_checkpoint_json(second, payload)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(json.loads(first.read_text(encoding="utf-8")), payload)

    def test_unavailable_method_is_retained_not_dropped(self):
        for method_a_status, method_b_status in (("unavailable", "available"), ("available", "unavailable")):
            with self.subTest(method_a_status=method_a_status, method_b_status=method_b_status):
                payload = _payload(method_a_status, method_b_status)
                self.assertEqual(payload["method_a"]["status"], method_a_status)
                self.assertEqual(payload["method_b"]["status"], method_b_status)

    def test_checkpoint_rejects_opaque_objects_nonfinite_values_and_corrections(self):
        with self.assertRaisesRegex(ValueError, "not_json_safe"):
            checkpoint.write_pion_hgcer_refinement_checkpoint_json("ignored.json", {"object": object()})
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            checkpoint.write_pion_hgcer_refinement_checkpoint_json("ignored.json", {"value": float("nan")})
        method = _method("b")
        method["cells"][0]["C_final"] = 1.0
        with self.assertRaisesRegex(ValueError, "forbidden_field"):
            checkpoint.build_pion_hgcer_refinement_checkpoint(
                setting={}, phase_a=_phase(), method_a=_method("a"), method_b=method
            )

    def test_filename_rule_and_a_b_neutral_schema(self):
        self.assertEqual(
            checkpoint.pion_hgcer_refinement_checkpoint_filename("Left", "Q4p4W2p74", "lowe"),
            "Left_pion_hgcer_refinement_checkpoint_Q4p4W2p74_lowe.json",
        )
        with self.assertRaises(ValueError):
            checkpoint.pion_hgcer_refinement_checkpoint_filename("Left", "../bad", "lowe")
        payload = _payload()
        forbidden = {"method_comparison", "method_agreement", "phase_d", "classification", "C_B", "C_final"}
        def walk(value):
            if isinstance(value, dict):
                self.assertTrue(forbidden.isdisjoint(value))
                for child in value.values():
                    walk(child)
            elif isinstance(value, list):
                for child in value:
                    walk(child)
        walk(payload)

    def test_runtime_writes_checkpoint_after_independent_method_b_handoff(self):
        source = (REPO_ROOT / "src/cuts/rand_sub.py").read_text(encoding="utf-8")
        method_b_position = source.index("build_pion_hgcer_method_b(")
        checkpoint_position = source.index(
            "build_pion_hgcer_refinement_checkpoint("
        )
        writer_position = source.index(
            "write_pion_hgcer_refinement_checkpoint_json("
        )
        self.assertLess(method_b_position, checkpoint_position)
        self.assertLess(checkpoint_position, writer_position)
        self.assertIn('histDict["pion_hgcer_refinement_checkpoint_artifacts"]', source)
        self.assertNotIn("method_comparison", source[method_b_position:writer_position])


if __name__ == "__main__":
    unittest.main()
