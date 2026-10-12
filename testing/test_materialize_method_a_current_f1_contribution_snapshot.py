"""Actual frozen builders with a genuinely hgcer3-supported synthetic population."""

from contextlib import redirect_stderr, redirect_stdout
import copy
import io
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import zipfile

from testing import materialize_method_a_current_f1_contribution_snapshot as owner
import test_pion_hgcer_method_a_tphi_propagation as fixtures


def population():
    artifacts = fixtures._artifacts()
    for artifact in artifacts:
        contract = artifact["contract"]
        for row in contract["method_a_training_records"] + contract["application_records"]:
            low = row.get("response_class") == "low"
            index = row["entry_index"] % 100
            # Independent nuisance coordinates: delta/track cannot classify response.
            # HGC x alone carries the signal, yielding unique hgcer3 without editing
            # any recommendation, optimizer, fingerprint or scientific gate.
            row["SHMS_delta"] = math.sin(index * 1.7)
            row["SHMS_xptar"] = math.cos(index * 2.3)
            row["SHMS_yptar"] = math.sin(index * 3.1)
            row["P_hgcer_xAtCer"] = (-1 if low else 1) + .05 * math.cos(index * 1.3)
            row["P_hgcer_yAtCer"] = math.sin(index * 2.7)
            if "baseline_pion_weight_w0" in row:
                row["analysis_MM"] = 1.12
                row["analysis_t"] = (row["t_low"] + row["t_high"]) / 2
        fixtures.f1_fixtures._seal_f1_artifact(artifact)
    return artifacts


class MaterializationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = population()

    def prepare(self, root, artifacts=None):
        for name in ("f1", "outputs", "bundles"):
            (root/name).mkdir()
        paths = {}
        for role, artifact in zip(owner.SETTINGS, self.inputs if artifacts is None else artifacts):
            path = root / "f1" / owner.canonical_name(role)
            path.write_bytes(owner.json_bytes(artifact))
            paths[role] = path
        return paths

    def args(self, root):
        return owner.parser().parse_args(["--f1-dir", str(root/"f1"), "--kinematic", owner.KINEMATIC,
            "--outdir", str(root/"outputs"), "--bundle-dir", str(root/"bundles")])

    def assert_empty(self, root):
        self.assertEqual(list((root/"outputs").iterdir()), [])
        self.assertEqual(list((root/"bundles").iterdir()), [])

    def test_actual_builders_complete_archive_locator_and_downstream_preflight(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); paths = self.prepare(root)
            original = {s: p.read_bytes() for s, p in paths.items()}
            pins = copy.deepcopy((owner.f4.ACCEPTED_F3_RUNTIME_AUTHORITY_BY_KINEMATIC,
                owner.debug.lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC))
            locator = owner.run(self.args(root))
            package, = (root/"bundles").iterdir()
            receipt, data = owner.downstream.preflight(locator)
            self.assertEqual(receipt["preflight"], "PASS")
            self.assertIsNotNone(data)
            self.assertEqual(original, {s: p.read_bytes() for s, p in paths.items()})
            self.assertEqual(pins, (owner.f4.ACCEPTED_F3_RUNTIME_AUTHORITY_BY_KINEMATIC,
                owner.debug.lineage.F6_3_CANDIDATE_F4_VALIDATION_AUTHORITY_BY_KINEMATIC))
            with zipfile.ZipFile(package) as archive:
                self.assertEqual(archive.namelist(), list(owner.MEMBERS))
                self.assertIsNone(archive.testzip())
                members = {n: archive.read(n) for n in owner.MEMBERS}
            manifest = owner.strict_json(members["manifest.json"])
            self.assertEqual(set(manifest["f1_inputs"]), set(owner.SETTINGS))
            self.assertEqual(manifest["canonical_owner_equivalence"], "BLOCKED")
            for key in ("canonical_owner_equivalence_evaluated", "method_a_promoted",
                        "production_application_performed", "production_objects_mutated",
                        "method_b_numerical_dependency"):
                self.assertIs(manifest[key], False)
            for raw in members.values():
                self.assertNotIn(str(root).encode(), raw)
            f2 = owner.strict_json(members[owner.NAMES["f2"]])
            self.assertEqual(f2["representation"]["recommendation"]["recommended_basis"], "hgcer3")
            self.assertEqual([s["overall_candidate_passed"] for s in f2["representation"]["candidate_summaries"]],
                             [False, False, True, False])
            f4 = owner.strict_json(members[owner.NAMES["f4"]])
            self.assertEqual(len(f4["correction"]["parents"]), 15)
            self.assertTrue(all(p["closure_passed"] for p in f4["correction"]["parents"]))
            a, h, _, _ = owner.load_f1(paths)
            source = manifest["source"]
            self.assertEqual(owner.verify_zip(package, a, h, source, members), owner.sha256(package.read_bytes()))
            # The earlier seven-input combination must fail through the UNCHANGED
            # downstream preflight, after a current F.1 raw byte change.
            paths["Center-highe"].write_bytes(original["Center-highe"] + b" ")
            failed, data = owner.downstream.preflight(locator)
            self.assertEqual(failed["preflight"], "FAIL")
            self.assertIn("f3_fingerprint_input_content_mismatch", failed["first_failure"])
            self.assertIsNone(data)
            self.assertFalse(failed["outputs_written"])

    def test_input_roles_names_schemas_fingerprints_and_geometry(self):
        for mode in ("missing", "duplicate", "misassigned", "filename", "schema", "fingerprint", "geometry"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); paths = self.prepare(root)
                if mode == "missing": paths.pop("Left-lowe")
                elif mode == "duplicate": paths["Left-highe"] = paths["Left-lowe"]
                elif mode == "filename":
                    p = paths["Left-lowe"]; q = p.with_name("wrong.json"); p.rename(q); paths["Left-lowe"] = q
                else:
                    p = paths["Left-lowe"]; value = owner.strict_json(p.read_bytes())
                    if mode == "misassigned": value["setting"]["phi_setting"] = "Right"
                    elif mode == "schema": value["schema_version"] = "wrong"
                    elif mode == "fingerprint": value["contract"]["fingerprint"] = "f"*64
                    else:
                        value["contract"]["phi_edges"][0] = -179.0
                        fixtures.f1_fixtures._seal_f1_artifact(value)
                    p.write_bytes(owner.json_bytes(value))
                with self.assertRaises((ValueError, KeyError)):
                    owner.load_f1(paths)
                self.assert_empty(root)

    def test_duplicate_roles_json_and_nonfinite_json(self):
        with self.assertRaises(ValueError):
            owner.shared._f1_inputs(["Left-lowe=x"]*5)
        for raw in (b'{"x":1,"x":2}', b'{"x":NaN}', b'{"x":Infinity}', b'{"x":1e999}'):
            with self.subTest(raw=raw), self.assertRaises(ValueError): owner.strict_json(raw)

    def test_file_identity_alias_and_replaced_original(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); paths = self.prepare(root)
            _, hashes, resolved, identities = owner.load_f1(paths)
            path = paths["Left-lowe"]
            path.write_bytes(path.read_bytes() + b" ")
            with self.assertRaisesRegex(ValueError, "f1_changed_during_materialization"):
                owner.recheck_f1(resolved, hashes, identities)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); paths = self.prepare(root)
            target = paths["Left-highe"]; target.unlink(); target.hardlink_to(paths["Left-lowe"])
            with self.assertRaisesRegex(ValueError, "f1_duplicate_file"):
                owner.load_f1(paths)

    def test_roots_collision_and_no_overwrite(self):
        with self.assertRaisesRegex(ValueError, "inside_checkout"): owner.external_root(owner.REPO_ROOT)
        with self.assertRaises(ValueError): owner.external_root(owner.REPO_ROOT/"nonexistent-output")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); self.prepare(root)
            final = root/"bundles"/"kaonlt-current-f1-snapshot-fixed.zip"; final.write_bytes(b"keep")
            with mock.patch.object(owner.uuid, "uuid4", return_value=mock.Mock(hex="fixed")):
                with self.assertRaisesRegex(ValueError, "output_already_exists"): owner.run(self.args(root))
            self.assertEqual(final.read_bytes(), b"keep")
            self.assertEqual(list((root/"outputs").iterdir()), [])

    def test_scientific_and_publication_failures_leave_no_success_artifact(self):
        actual_f3 = owner.f3.build_pion_hgcer_method_a_acceptance_map_artifact
        actual_verify = owner.verify_zip
        actual_recheck = owner.recheck_f1
        cases = ("f2_support", "f3_convergence", "f4_closure", "zip_corruption", "final_zip", "changed_f1")
        for mode in cases:
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root = Path(directory); paths = self.prepare(root)
                keep = root/"keep"; keep.write_bytes(b"unrelated")
                def build(*args, **kwargs):
                    result = actual_f3(*args, **kwargs)
                    if mode == "f3_convergence": result["acceptance_map"]["available"] = False
                    return result
                def verify(path, *args):
                    if mode == "zip_corruption": path.write_bytes(b"broken ZIP")
                    if mode == "final_zip" and path.parent == root/"bundles":
                        raise ValueError("final_zip_gate_failure")
                    return actual_verify(path, *args)
                def recheck(*args):
                    if mode == "changed_f1":
                        p = paths["Left-lowe"]; p.write_bytes(p.read_bytes()+b" ")
                    return actual_recheck(*args)
                if mode == "f2_support":
                    value = owner.strict_json(paths["Left-lowe"].read_bytes())
                    value["contract"]["method_a_training_records"].pop(0)
                    fixtures.f1_fixtures._seal_f1_artifact(value)
                    paths["Left-lowe"].write_bytes(owner.json_bytes(value))
                closure = owner.debug.reproduce_snapshot
                with mock.patch.object(owner.f3, "build_pion_hgcer_method_a_acceptance_map_artifact", side_effect=build), \
                     mock.patch.object(owner, "verify_zip", side_effect=verify), \
                     mock.patch.object(owner, "recheck_f1", side_effect=recheck), \
                     mock.patch.object(owner.debug, "reproduce_snapshot", side_effect=(ValueError("closure_failed") if mode == "f4_closure" else closure)):
                    with self.assertRaises((ValueError, zipfile.BadZipFile)): owner.run(self.args(root))
                self.assert_empty(root)
                self.assertEqual(keep.read_bytes(), b"unrelated")

    def test_persisted_stage_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); paths = self.prepare(root); locator = owner.run(self.args(root))
            artifacts, hashes, _, _ = owner.load_f1(paths)
            raw = {s: (locator.parent/n).read_bytes() for s, n in owner.NAMES.items()}
            cases = (("f2", "fingerprint"), ("f2", "group"), ("f3", "model"),
                     ("f3", "f2_sha"), ("f4", "sentinel"), ("f4", "reference"),
                     ("f4", "parent"), ("f4", "closure"), ("f4", "nan"))
            for stage, mode in cases:
                with self.subTest(stage=stage, mode=mode):
                    value = owner.strict_json(raw[stage]); body = value[owner.BODY[stage]]
                    if mode == "fingerprint": value["artifact_fingerprint"] = "f"*64
                    elif mode == "group": body["groups"].pop()
                    elif mode == "model": body["models"][0]["logistic_fit"]["coefficients"][0] += .2
                    elif mode == "f2_sha": body["f2_source_file_sha256"] = "f"*64
                    elif mode == "sentinel": body["f3_runtime_authority"]["accepted"]["farm_source_head"] = "1"*40
                    elif mode == "reference": body["f3_runtime_authority"]["accepted"]["source_file_sha256"] = "f"*64
                    elif mode == "parent": body["parents"].pop()
                    elif mode == "closure": body["parents"][0]["adjusted_parent_sum"] += .1
                    else: body["parents"][0]["parent_normalization"] = float("inf")
                    if mode == "nan":
                        data = raw[stage].replace(b'"parent_normalization": ', b'"parent_normalization": 1e999, "unused": ', 1)
                    else: data = owner.json_bytes(value)
                    with self.assertRaises((ValueError, KeyError)):
                        owner.verify_stages(artifacts, hashes, {**raw, stage: data})

    def test_manifest_reseal_is_not_sufficient_and_zip_bytes_are_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); paths = self.prepare(root); locator = owner.run(self.args(root))
            package, = (root/"bundles").iterdir()
            artifacts, hashes, _, _ = owner.load_f1(paths)
            members = {n: (locator.parent/n).read_bytes() for n in owner.MEMBERS}
            repeat = root/"repeat.zip"; owner.create_zip(repeat, members)
            self.assertEqual(package.read_bytes(), repeat.read_bytes())
            manifest = owner.strict_json(members["manifest.json"]); source = manifest["source"]
            manifest["checks"]["exact_f4_reproduction"] = False; manifest.pop("fingerprint")
            changed = {**members, "manifest.json": owner.json_bytes(owner.debug.seal(manifest))}
            bad = root/"bad.zip"; owner.create_zip(bad, changed)
            with self.assertRaisesRegex(ValueError, "manifest_or_fingerprint"):
                owner.verify_zip(bad, artifacts, hashes, source, changed)

    def test_resealed_inventory_closure_and_embedded_reference_still_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); paths = self.prepare(root); locator = owner.run(self.args(root))
            artifacts, hashes, _, _ = owner.load_f1(paths)
            raw = {s: (locator.parent/n).read_bytes() for s, n in owner.NAMES.items()}
            for stage, mode in (("f2", "group"), ("f3", "model"), ("f4", "parent"),
                                ("f4", "closure"), ("f4", "reference"), ("f4", "sentinel")):
                with self.subTest(stage=stage, mode=mode):
                    value = owner.strict_json(raw[stage]); body = value[owner.BODY[stage]]
                    if mode == "group": body["groups"][0]["canonical_t_index"] = 4
                    elif mode == "model": body["models"][0]["canonical_t_index"] = 4
                    elif mode == "parent": body["parents"][0]["canonical_t_index"] = 4
                    elif mode == "closure": body["parents"][0]["adjusted_parent_sum"] += .1
                    else:
                        field = "farm_source_head" if mode == "sentinel" else "map_fingerprint"
                        body["f3_runtime_authority"]["accepted"][field] = "1"*(40 if mode == "sentinel" else 64)
                    if stage == "f2": body["fingerprint_inputs"]["groups"] = body["groups"]
                    elif stage == "f3": body["fingerprint_inputs"]["models"] = body["models"]
                    else:
                        body["fingerprint_inputs"]["parents"] = body["parents"]
                        body["fingerprint_inputs"]["f3_runtime_authority"] = body["f3_runtime_authority"]
                    body["fingerprint"] = owner.debug.fingerprint(body["fingerprint_inputs"])
                    field = {"f2": "representation_fingerprint", "f3": "map_fingerprint", "f4": "correction_fingerprint"}[stage]
                    value["artifact_fingerprint"] = owner.debug.fingerprint({"schema_version": value["schema_version"],
                        field: body["fingerprint"], "input_paths": {}})
                    with self.assertRaises(ValueError):
                        owner.verify_stages(artifacts, hashes, {**raw, stage: owner.json_bytes(value)})

    def test_event_records_and_private_paths_are_rejected_but_empty_public_provenance_is_preserved(self):
        owner.transferable({"provenance": {"input_paths": {}}})
        for value in ({"application_records": []}, {"correction_factors": [1]},
                      {"provenance": {"input_paths": {"f1": "private"}}},
                      {"path": "/private/data"}, {"path": "C:\\private\\data"}):
            with self.subTest(value=value), self.assertRaises(ValueError): owner.transferable(value)

    def test_cli_failure_has_no_success_stdout_or_private_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); self.prepare(root); stdout, stderr = io.StringIO(), io.StringIO()
            (root/"f1"/owner.canonical_name("Left-lowe")).unlink()
            with redirect_stdout(stdout), redirect_stderr(stderr):
                result = owner.main(["--f1-dir", str(root/"f1"), "--kinematic", owner.KINEMATIC,
                    "--outdir", str(root/"outputs"), "--bundle-dir", str(root/"bundles")])
            self.assertEqual(result, 1); self.assertEqual(stdout.getvalue(), "")
            self.assertNotIn(str(root), stderr.getvalue()); self.assert_empty(root)


if __name__ == "__main__":
    unittest.main()
