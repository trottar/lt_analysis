"""Synthetic contract-shaped forensic fixtures; no patched scientific gates."""

from copy import deepcopy
import contextlib
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from testing import diagnose_e8_4_f1_baseline_reproducibility as diagnostic


def seal(artifact):
    c = artifact["contract"]
    values = {"method_a_training_population_fingerprint": diagnostic.f1._hash(c["method_a_training_records"]),
              "application_population_fingerprint": diagnostic.f1._hash(c["application_records"]),
              "acceptance_feature_metadata_fingerprint": diagnostic.f1._hash(c["feature_metadata"]),
              "application_child_assignment_projection_fingerprint": diagnostic.f1._hash([{k: r[k] for k in diagnostic.CHILD} for r in c["application_records"]])}
    c.update(values)
    c["fingerprint_inputs"] = {"schema_version": c["schema_version"], "fingerprint_schema_version": c["fingerprint_schema_version"],
                               **{k: c[k] for k in diagnostic.PROVENANCE}, "host_state": c["host_state"], "source_target_state": c["source_target_state"],
                               **{k: c[k] for k in ("t_edges", "delta_edges", "phi_edges")}, **values,
                               "method_a_closure": c["method_a_training_summary"]["by_t_delta"], "feature_metadata": c["feature_metadata"]}
    c["fingerprint"] = diagnostic.f1._hash(c["fingerprint_inputs"])
    return artifact


def make_f1():
    rows, training = [], []
    for t in range(3):
        for index, (source, coefficient, weight, phi) in enumerate([
                ("prompt", 1.0, 0.5, -90.0), ("rand", -0.5, 1.0, 90.0),
                ("dummy", -0.25, 0.0, 0.0), ("dummy_rand", 0.125, 2.0, 181.0)]):
            p = None if phi > 180 else 0 if phi < 0 else 1
            row = {"source_label": source, "entry_index": t * 10 + index, "t_index": t,
                   "t_low": float(t), "t_high": float(t + 1), "analysis_t": t + 0.5,
                   "analysis_MM": 1.1 + t * 0.01, "phi_index": p,
                   "phi_low": None if p is None else [-180.0, 0.0][p],
                   "phi_high": None if p is None else [0.0, 180.0][p], "phi_degrees": phi,
                   "phi_status": "outside_phi" if p is None else "inside_phi",
                   "signed_source_coefficient": coefficient, "baseline_pion_weight_w0": weight,
                   "signed_baseline_event_contribution": coefficient * weight,
                   "P_hgcer_npeSum": 3.0, "nommcuts": True, "allcuts": True,
                   **{f: 0.1 for f in diagnostic.FEATURES}}
            rows.append(row)
            if source == "prompt":
                training.append({**row, "response_class": "control"})
                training.append({**row, "entry_index": 1000 + t, "response_class": "low", "P_hgcer_npeSum": 1.0})
    metadata = {"primary_acceptance_features": list(diagnostic.FEATURES),
                "training_population": "prompt_noRF_nommcuts_P_hgcer_npeSum_gt_0",
                "training_low_definition": "0_lt_P_hgcer_npeSum_le_2", "training_control_definition": "P_hgcer_npeSum_gt_2",
                "application_population": "authoritative_physical_pion_control_P_hgcer_npeSum_gt_2",
                "parent_coordinate": "canonical_t", "downstream_yield_coordinates": ["canonical_t", "canonical_phi"],
                "phi_is_training_feature": False, "method_b_numerical_dependency": False,
                "probability_map_constructed": False, "weight_adjustment_constructed": False,
                "future_normalization_policy": "future_parent_t_only_no_tphi_child_renormalization", "absolute_leakage_probability_claimed": False}
    flags = {"non_authoritative": True, "production_objects_mutated": False, "refinement_applied": False,
             "production_application_performed": False, "event_application_performed": False}
    c = {**flags, "schema_version": diagnostic.f1.METHOD_A_ACCEPTANCE_EVENT_CONTRACT_SCHEMA_VERSION,
         "fingerprint_schema_version": diagnostic.f1.METHOD_A_ACCEPTANCE_EVENT_CONTRACT_FINGERPRINT_SCHEMA_VERSION,
         "status": "available", "available": True, "reason": None, "diagnostic_stage": "complete",
         "method_b_numerical_dependency": False, "future_weight_adjustment_constructed": False,
         "host_state": "post_proton", "source_target_state": "post_proton_noRF", "coordinate_fingerprint": "coordinate-test",
         **{k: "synthetic-" + k for k in diagnostic.PROVENANCE}, "feature_metadata": metadata,
         "t_edges": [0.0, 1.0, 2.0, 3.0], "delta_edges": [-10.0, 10.0], "phi_edges": [-180.0, 0.0, 180.0],
         "application_records": rows, "method_a_training_records": training, "method_a_training_summary": {"by_t_delta": []}}
    return seal({**flags, "schema_version": diagnostic.f1.METHOD_A_ACCEPTANCE_EVENT_CONTRACT_ARTIFACT_SCHEMA_VERSION,
                 "setting": dict(diagnostic.SETTING), "contract": c})


def make_f4(f1):
    inv = diagnostic.inventory(diagnostic.validate_f1(f1))["f4_consumed_inside_phi"]["by_t"]
    parents = [{"setting_id": s, "setting": {**diagnostic.SETTING, "phi_setting": s.split("-")[0],
                                           "epsilon_filename_token": s.split("-")[1], "epsilon_setting": "low" if s.endswith("lowe") else "high"},
                "canonical_t_index": t, "canonical_t_low": float(t), "canonical_t_high": float(t + 1), "application_event_count": inv[t]["count"],
                "baseline_parent_sum": inv[t]["signed_baseline_sum"], "absolute_baseline_parent_sum": inv[t]["absolute_baseline_sum"]}
               for s in ("Left-lowe", "Left-highe", "Center-lowe", "Center-highe", "Right-highe") for t in range(3)]
    flags = {"non_authoritative": True, "production_objects_mutated": False, "production_application_performed": False,
             "method_b_numerical_dependency": False, "event_correction_persisted": False, "child_renormalization_performed": False,
             "correction_applied_to_production": False}
    body = {**flags, "schema_version": "pion_hgcer_method_a_parent_preserving_correction/v1",
            "fingerprint_schema_version": "pion_hgcer_method_a_parent_preserving_correction_fingerprint/v1",
            "status": "available", "available": True, "accepted_basis": "hgcer3", "parents": parents}
    artifact = {**flags, "schema_version": "pion_hgcer_method_a_parent_preserving_correction_artifact/v1",
                "correction": body, "provenance": {"input_paths": {}}}
    return seal_f4(artifact)


def seal_f4(artifact):
    c = artifact["correction"]
    c["fingerprint_inputs"] = {"parents": deepcopy(c["parents"])}
    c["fingerprint"] = diagnostic.f1._hash(c["fingerprint_inputs"])
    artifact["artifact_fingerprint"] = diagnostic.f1._hash({"schema_version": artifact["schema_version"], "correction_fingerprint": c["fingerprint"], "input_paths": artifact["provenance"]["input_paths"]})
    return artifact


def encoded(obj):
    return (json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def digest(obj):
    return hashlib.sha256(encoded(obj)).hexdigest()


class DiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.f1 = make_f1()
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.current = self.root / "current.json"
        self.current.write_bytes(encoded(self.f1))
        self.sha = hashlib.sha256(self.current.read_bytes()).hexdigest()
        self.output = self.root / "diagnostic.json"

    def cli(self, extra=(), output=None):
        args = ["--current-f1", str(self.current), "--expected-current-sha256", self.sha, "--output", str(output or self.output), *extra]
        before = self.current.read_bytes()
        with contextlib.redirect_stderr(io.StringIO()):
            code = diagnostic.main(args)
        self.assertEqual(before, self.current.read_bytes())
        return code

    def test_inventory_cancellation_zero_and_outside(self):
        original = deepcopy(self.f1)
        result = diagnostic.build_comparison(self.f1, self.sha)
        self.assertEqual(self.f1, original)
        self.assertEqual(result["inventory"]["all_application"]["count"], 12)
        self.assertEqual(result["inventory"]["f4_consumed_inside_phi"]["count"], 9)
        for parent in result["recomputed_baseline_parents"]:
            t = parent["canonical_t_index"]
            rows = [r for r in self.f1["contract"]["application_records"] if r["t_index"] == t and r["phi_status"] == "inside_phi"]
            self.assertEqual(parent["baseline_parent_sum"], math.fsum(r["signed_baseline_event_contribution"] for r in rows))
            self.assertEqual(parent["absolute_baseline_parent_sum"], 1.0)
        self.assertEqual(result["inventory"]["all_application"]["baseline_weight_distribution"]["zero_count"], 3)
        self.assertEqual(result["event_level_source_attribution"], "NOT VERIFIED")

    def test_reordering_preserves_science_and_keeps_raw_roles(self):
        reordered = deepcopy(self.f1)
        reordered["contract"]["application_records"].reverse()
        reordered["contract"]["method_a_training_records"].reverse()
        seal(reordered)
        result = diagnostic.build_comparison(reordered, digest(reordered), self.f1, self.sha)
        self.assertNotEqual(result["input_sha256"]["current_f1"], result["input_sha256"]["reference_f1"])
        self.assertEqual(result["inventory"], result["reference_comparison"]["reference_inventory"])
        self.assertTrue(all(v["exact_equal"] for v in result["reference_comparison"]["variables"].values()))
        self.assertEqual(result["reference_role"], "historical_comparison_not_reviewed_authority")

    def test_coefficient_only_never_reports_weight_change(self):
        changed = deepcopy(self.f1)
        for row in changed["contract"]["application_records"]:
            row["signed_source_coefficient"] *= 1.00001
            row["signed_baseline_event_contribution"] = row["signed_source_coefficient"] * row["baseline_pion_weight_w0"]
        seal(changed)
        result = diagnostic.build_comparison(changed, digest(changed), self.f1, self.sha)["reference_comparison"]["variables"]
        self.assertEqual(result["signed_source_coefficient"]["changed_count"], 12)
        self.assertEqual(result["baseline_pion_weight_w0"]["changed_count"], 0)
        self.assertEqual(result["signed_source_coefficient"]["first_changed_ids"][0], ["dummy", 2])
        self.assertEqual(len(result["signed_source_coefficient"]["first_changed_ids"]), 10)

    def test_weight_only(self):
        changed = deepcopy(self.f1)
        row = changed["contract"]["application_records"][1]
        row["baseline_pion_weight_w0"] = 0.75
        row["signed_baseline_event_contribution"] = -0.375
        seal(changed)
        variables = diagnostic.build_comparison(changed, digest(changed), self.f1, self.sha)["reference_comparison"]["variables"]
        self.assertEqual(variables["baseline_pion_weight_w0"]["changed_count"], 1)
        self.assertEqual(variables["signed_source_coefficient"]["changed_count"], 0)
        self.assertEqual(variables["baseline_pion_weight_w0"]["signed_delta_sum"], -0.25)

    def test_assignment_missing_added_and_mm_changes(self):
        changed = deepcopy(self.f1)
        rows = changed["contract"]["application_records"]
        rows[1].update(t_index=1, t_low=1.0, t_high=2.0, analysis_t=1.5, analysis_MM=1.2, phi_index=0, phi_low=-180.0, phi_high=0.0, phi_degrees=-10.0)
        rows[2]["entry_index"] = 999
        seal(changed)
        comparison = diagnostic.build_comparison(changed, digest(changed), self.f1, self.sha)["reference_comparison"]
        self.assertEqual(comparison["missing_count"], 1)
        self.assertEqual(comparison["added_count"], 1)
        self.assertEqual(comparison["first_missing_ids"], [["dummy", 2]])
        for name in ("t_index", "analysis_MM", "phi_index"):
            self.assertEqual(comparison["variables"][name]["changed_count"], 1)

    def test_reviewed_f4_exact_mismatch_and_missing_reviewed_f1(self):
        f4 = make_f4(self.f1)
        # Same observed difference as the supplied t0 failure, no tolerance PASS.
        current_value, reviewed_value = 0.12072578661296547, 0.12072578658633845
        changed = deepcopy(self.f1)
        row = changed["contract"]["application_records"][0]
        row["baseline_pion_weight_w0"] = current_value / 2
        row["signed_baseline_event_contribution"] = row["baseline_pion_weight_w0"]
        row = changed["contract"]["application_records"][1]
        row["baseline_pion_weight_w0"] = current_value
        row["signed_baseline_event_contribution"] = -current_value / 2
        seal(changed)
        f4["correction"]["parents"][0]["absolute_baseline_parent_sum"] = reviewed_value
        seal_f4(f4)
        # Isolate synthetic fixture pins; never declare synthetic bytes to have
        # the production authority hash or patch an acceptance function.
        with patch.object(diagnostic, "REVIEWED_F4_SHA", digest(f4)):
            result = diagnostic.build_comparison(changed, digest(changed), reviewed_f4=f4, reviewed_f4_sha256=digest(f4))
        metric = result["reviewed_f4_parent_comparison"][0]["metrics"]["absolute_baseline_parent_sum"]
        self.assertFalse(metric["exact_equal"])
        self.assertAlmostEqual(metric["absolute_difference"], abs(current_value - reviewed_value), delta=1e-16)
        self.assertEqual(result["event_level_source_attribution"], "NOT VERIFIED")
        self.assertFalse(result["runtime_acceptance_granted"])

    def test_producer_filename_token_setting(self):
        a = deepcopy(self.f1)
        a["setting"].update(Q2="4p4", W="2p74")
        self.assertEqual(diagnostic.build_comparison(a, digest(a))["setting"], diagnostic.SETTING)
        f4 = make_f4(a)
        for row in f4["correction"]["parents"]:
            row["setting"].update(Q2="4p4", W="2p74")
        seal_f4(f4)
        with patch.object(diagnostic, "REVIEWED_F4_SHA", digest(f4)):
            diagnostic.build_comparison(a, digest(a), reviewed_f4=f4, reviewed_f4_sha256=digest(f4))

    def test_valid_cli_byte_determinism_and_no_mutation(self):
        self.assertEqual(self.cli(), 0)
        second = self.root / "second.json"
        self.assertEqual(self.cli(output=second), 0)
        self.assertEqual(self.output.read_bytes(), second.read_bytes())
        self.assertEqual(hashlib.sha256(self.current.read_bytes()).hexdigest(), self.sha)
        self.assertEqual(json.loads(self.output.read_bytes())["schema_version"], diagnostic.SCHEMA)

    def test_reference_cli_role_and_immutability(self):
        reference = deepcopy(self.f1)
        reference["contract"]["application_records"].reverse()
        seal(reference)
        path = self.root / "reference.json"
        path.write_bytes(encoded(reference))
        self.assertEqual(self.cli(["--reference-f1", str(path), "--expected-reference-sha256", digest(reference)]), 0)
        self.assertEqual(path.read_bytes(), encoded(reference))

    def test_builder_requires_fixed_reviewed_f4_authority(self):
        expected = "79e7ceda7221cbeeead4ed5bc306b0e0e670741a27beaa980e22349c555e96d7"
        self.assertEqual(diagnostic.REVIEWED_F4_SHA, expected)
        f4 = make_f4(self.f1)
        before = deepcopy(f4)
        fixture_sha = digest(f4)
        with patch.object(diagnostic, "REVIEWED_F4_SHA", fixture_sha):
            result = diagnostic.build_comparison(self.f1, self.sha, reviewed_f4=f4, reviewed_f4_sha256=fixture_sha)
        self.assertEqual(diagnostic.REVIEWED_F4_SHA, expected)
        self.assertEqual(result["input_sha256"]["reviewed_f4"], fixture_sha)
        self.assertTrue(all(metric["exact_equal"] for parent in result["reviewed_f4_parent_comparison"] for metric in parent["metrics"].values()))
        self.assertFalse(result["runtime_acceptance_granted"])
        self.assertEqual(f4, before)
        # Both a different raw serialization and changed, resealed parent data
        # remain self-consistent but cannot acquire reviewed-authority status.
        changed = deepcopy(f4)
        changed["correction"]["parents"][0]["baseline_parent_sum"] = 1.0
        seal_f4(changed)
        with patch.object(diagnostic, "REVIEWED_F4_SHA", fixture_sha):
            diagnostic.validate_f4(changed, self.f1["contract"]["t_edges"])
            with self.assertRaisesRegex(ValueError, "reviewed F.4 authority SHA-256 mismatch"):
                diagnostic.build_comparison(self.f1, self.sha, reviewed_f4=changed, reviewed_f4_sha256=digest(changed))
        for alternative in (f4, changed):
            diagnostic.validate_f4(alternative, self.f1["contract"]["t_edges"])
            with self.assertRaisesRegex(ValueError, "reviewed F.4 authority SHA-256 mismatch"):
                diagnostic.build_comparison(self.f1, self.sha, reviewed_f4=alternative, reviewed_f4_sha256=digest(alternative))

    def test_cli_isolated_fixture_pin_accepts_only_pinned_f4_bytes(self):
        production_pin = diagnostic.REVIEWED_F4_SHA
        f4 = make_f4(self.f1)
        fixture_sha = digest(f4)
        path = self.root / "reviewed-fixture.json"
        path.write_bytes(encoded(f4))
        changed = deepcopy(f4)
        changed["correction"]["parents"][0]["baseline_parent_sum"] = 1.0
        seal_f4(changed)
        diagnostic.validate_f4(changed, self.f1["contract"]["t_edges"])
        alternative = self.root / "alternative.json"
        alternative.write_bytes(encoded(changed))
        before = {p: p.read_bytes() for p in (self.current, path, alternative)}
        with patch.object(diagnostic, "REVIEWED_F4_SHA", fixture_sha):
            self.assertEqual(self.cli(["--reviewed-f4", str(path), "--expected-reviewed-f4-sha256", fixture_sha]), 0)
            result = json.loads(self.output.read_bytes())
            self.assertEqual(result["input_sha256"]["reviewed_f4"], fixture_sha)
            self.assertTrue(all(metric["exact_equal"] for parent in result["reviewed_f4_parent_comparison"] for metric in parent["metrics"].values()))
            for supplied_sha in (digest(changed), fixture_sha):
                rejected_output = self.root / (supplied_sha + ".json")
                self.assertEqual(self.cli(["--reviewed-f4", str(alternative), "--expected-reviewed-f4-sha256", supplied_sha], output=rejected_output), 2)
                self.assertFalse(rejected_output.exists())
        self.assertEqual(diagnostic.REVIEWED_F4_SHA, production_pin)
        self.assertEqual({p: p.read_bytes() for p in before}, before)

    def test_cli_rejects_self_consistent_alternative_f4_before_reads(self):
        f4 = make_f4(self.f1)
        diagnostic.validate_f4(f4, self.f1["contract"]["t_edges"])
        path = self.root / "f4.json"
        raw = encoded(f4)
        path.write_bytes(raw)
        with patch.object(diagnostic, "read_pinned", wraps=diagnostic.read_pinned) as read:
            self.assertEqual(self.cli(["--reviewed-f4", str(path), "--expected-reviewed-f4-sha256", digest(f4)]), 2)
            read.assert_not_called()
        self.assertEqual(path.read_bytes(), raw)
        self.assertFalse(self.output.exists())

    def test_cli_fixed_pin_reaches_independent_byte_check(self):
        f4 = make_f4(self.f1)
        path = self.root / "f4.json"
        raw = encoded(f4)
        path.write_bytes(raw)
        # The genuine reviewed bytes are not a local fixture. The correct pin
        # must pass argument preflight, then reject these different real bytes.
        with patch.object(diagnostic, "read_pinned", wraps=diagnostic.read_pinned) as read:
            self.assertEqual(self.cli(["--reviewed-f4", str(path), "--expected-reviewed-f4-sha256", diagnostic.REVIEWED_F4_SHA]), 2)
            self.assertEqual(read.call_args_list[-1].args, (path, diagnostic.REVIEWED_F4_SHA))
            self.assertEqual(read.call_count, 2)
        self.assertEqual(path.read_bytes(), raw)
        self.assertFalse(self.output.exists())

    def test_bad_hash_and_pairing(self):
        for extra in (["--expected-current-sha256", "0" * 64], ["--expected-current-sha256", "bad"], ["--reference-f1", str(self.current)], ["--expected-reference-sha256", self.sha], ["--reviewed-f4", str(self.current)], ["--expected-reviewed-f4-sha256", self.sha]):
            with self.subTest(extra=extra):
                self.assertEqual(self.cli(extra), 2)
                self.assertFalse(self.output.exists())

    def test_malformed_duplicate_nonfinite_json(self):
        for raw in (b'{', b'[]', b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":Infinity}', b'{"x":1e400}'):
            with self.subTest(raw=raw):
                self.current.write_bytes(raw)
                self.sha = hashlib.sha256(raw).hexdigest()
                self.assertEqual(self.cli(), 2)
                self.assertFalse(self.output.exists())

    def test_invalid_f1_semantics_no_output(self):
        mutations = [
            lambda a: a["setting"].update(phi_setting="Right"),
            lambda a: a["setting"].update(Q2=3.0),
            lambda a: a["setting"].pop("W"),
            lambda a: a["contract"].update(schema_version="bad"),
            lambda a: a["contract"].update(production_objects_mutated=True),
            lambda a: a["contract"]["application_records"].append(deepcopy(a["contract"]["application_records"][0])),
            lambda a: a["contract"]["application_records"][0].pop("entry_index"),
            lambda a: a["contract"]["application_records"][0].update(t_index=3),
            lambda a: a["contract"]["application_records"][0].update(t_low=-1),
            lambda a: a["contract"]["application_records"][0].update(analysis_t=1.9),
            lambda a: a["contract"]["application_records"][0].update(phi_index=1),
            lambda a: a["contract"]["application_records"][0].update(phi_degrees=90),
            lambda a: a["contract"]["application_records"][3].update(phi_degrees=180),
            lambda a: a["contract"]["application_records"][0].pop("baseline_pion_weight_w0"),
            lambda a: a["contract"]["application_records"][0].update(signed_baseline_event_contribution=100),
            lambda a: a["contract"]["application_records"][0].update(baseline_pion_weight_w0=-1),
            lambda a: a["contract"]["application_records"][0].update(entry_index=True),
            lambda a: a["contract"]["application_records"][0].update(P_hgcer_npeSum=1),
            lambda a: a["contract"].update(fingerprint="wrong"),
        ]
        for index, mutate in enumerate(mutations):
            with self.subTest(index=index):
                artifact = deepcopy(self.f1)
                mutate(artifact)
                # Fingerprint-preserving edit is deliberately NOT used to hide invalid input.
                self.current.write_bytes(encoded(artifact))
                self.sha = digest(artifact)
                self.assertEqual(self.cli(), 2)
                self.assertFalse(self.output.exists())

    def test_wrong_f4_schema_and_parent_inventory(self):
        for mode in ("schema", "inventory", "fingerprint", "flag", "setting", "geometry"):
            f4 = make_f4(self.f1)
            if mode == "schema":
                f4["schema_version"] = "bad"
            elif mode == "inventory":
                f4["correction"]["parents"].pop()
                seal_f4(f4)
            elif mode == "fingerprint":
                f4["correction"]["parents"][0]["baseline_parent_sum"] = 10
            elif mode == "flag":
                f4["production_objects_mutated"] = True
            elif mode == "setting":
                f4["correction"]["parents"][0]["setting"]["Q2"] = "3p0"
                seal_f4(f4)
            else:
                f4["correction"]["parents"][0]["canonical_t_low"] = -1
                seal_f4(f4)
            path = self.root / "f4.json"
            path.write_bytes(encoded(f4))
            # Fixture pin isolation allows schema/geometry checks to execute
            # beyond the unchanged production authority gate.
            with patch.object(diagnostic, "REVIEWED_F4_SHA", digest(f4)):
                with self.assertRaises(ValueError):
                    diagnostic.build_comparison(self.f1, self.sha, reviewed_f4=f4, reviewed_f4_sha256=digest(f4))
                self.assertEqual(self.cli(["--reviewed-f4", str(path), "--expected-reviewed-f4-sha256", digest(f4)]), 2)
            self.assertFalse(self.output.exists())

    def test_existing_output_alias_symlinks_and_missing_parent(self):
        self.output.write_bytes(b"preserve")
        self.assertEqual(self.cli(), 2)
        self.assertEqual(self.output.read_bytes(), b"preserve")
        self.assertEqual(self.cli(output=self.current), 2)
        self.assertEqual(self.cli(output=self.root / "absent" / "out.json"), 2)
        alias = self.root / "alias.json"
        alias.hardlink_to(self.current)
        self.assertEqual(self.cli(output=alias), 2)
        alias.unlink()
        alias.symlink_to(self.current)
        self.assertEqual(self.cli(output=alias), 2)
        alias.unlink()
        alias.symlink_to(self.root / "missing.json")
        self.assertEqual(self.cli(output=alias), 2)
        self.assertFalse((self.root / "missing.json").exists())

    def test_reference_identical_bytes_path_or_inode_rejected(self):
        ref = self.root / "ref.json"
        ref.write_bytes(self.current.read_bytes())
        self.assertEqual(self.cli(["--reference-f1", str(ref), "--expected-reference-sha256", self.sha]), 2)
        ref.unlink()
        ref.hardlink_to(self.current)
        self.assertEqual(self.cli(["--reference-f1", str(ref), "--expected-reference-sha256", self.sha]), 2)
        self.assertEqual(self.cli(["--reference-f1", str(self.current), "--expected-reference-sha256", self.sha]), 2)
        self.assertFalse(self.output.exists())

    def test_input_symlinks_are_refused(self):
        link = self.root / "link.json"
        link.symlink_to(self.current)
        self.assertEqual(self.cli(["--current-f1", str(link)]), 2)
        self.assertFalse(self.output.exists())

    def test_phi_last_edge_and_resealed_invalid_identity(self):
        a = deepcopy(self.f1)
        a["contract"]["application_records"][1]["phi_degrees"] = 180.0
        seal(a)
        diagnostic.validate_f1(a)
        a["contract"]["application_records"][1]["signed_baseline_event_contribution"] += 1e-8
        seal(a)
        with self.assertRaisesRegex(ValueError, "signed baseline identity"):
            diagnostic.validate_f1(a)

    def test_output_race_never_overwrites_and_temp_is_removed(self):
        original = diagnostic.os.link
        def race(src, dst):
            Path(dst).write_bytes(b"concurrent")
            return original(src, dst)
        with patch.object(diagnostic.os, "link", side_effect=race):
            self.assertEqual(self.cli(), 2)
        self.assertEqual(self.output.read_bytes(), b"concurrent")
        self.assertEqual(list(self.root.glob(".f1-diagnostic-*")), [])

    def test_no_scientific_runtime_imports_in_fresh_process(self):
        code = "from testing import diagnose_e8_4_f1_baseline_reproducibility; import sys; assert not any(x in sys.modules for x in ('ROOT','numpy','scipy'))"
        subprocess.run([sys.executable, "-B", "-c", code], cwd=Path(__file__).resolve().parents[1], check=True)


if __name__ == "__main__":
    unittest.main()
