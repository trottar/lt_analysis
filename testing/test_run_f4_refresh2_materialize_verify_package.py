"""Synthetic, farm-free checks for the tracked F.4.Refresh.2 execution owner."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
import zipfile

from testing import run_f4_refresh2_materialize_verify_package as owner


def encoded(value):
    return json.dumps(value, sort_keys=True).encode("utf-8")


class Fixture:
    def __init__(self, root: Path):
        self.repo = root / "repo"
        self.repo.mkdir()
        (self.repo / ".git").mkdir()
        self.profile_path = self.repo / owner.PROFILE_RELATIVE
        self.profile_path.parent.mkdir()
        self.profile_path.write_bytes((owner.REPO_ROOT / owner.PROFILE_RELATIVE).read_bytes())
        self.artifacts = root / "artifacts"
        self.artifacts.mkdir()
        self.globus = root / "globus"
        self.globus.mkdir()
        self.inputs = root / "inputs"
        self.inputs.mkdir()
        self.f1 = {}
        for alias in owner.ALIASES:
            self.f1[alias] = self.inputs / f"{alias}.json"
            self.f1[alias].write_bytes(encoded({"alias": alias}))
        self.accepted = {}
        for stage in owner.STAGES:
            self.accepted[stage] = self.inputs / f"accepted-{stage}.json"
            self.accepted[stage].write_bytes(encoded({"accepted": stage}))
        self.comparison = self.inputs / "reviewed-comparison.json"
        self.comparison.write_bytes(encoded({"reviewed": True}))
        self.args = argparse.Namespace(
            python=sys.executable, bundle_commit="a" * 40,
            f1=[f"{alias}={path}" for alias, path in self.f1.items()],
            accepted_f2=self.accepted["f2"], accepted_f3=self.accepted["f3"],
            accepted_f4=self.accepted["f4"], comparison=self.comparison,
            artifact_dir=self.artifacts, kinematic=owner.KINEMATIC,
            output="fresh-review.zip",
        )
        self.calls = []
        self.head = self.args.bundle_commit
        self.dirty = ""
        self.ancestor = True
        self.changed = "\n".join(sorted(owner.ALLOWED_COMMITTED | {"docs/memory/CURRENT.md"}))
        self.materializer_code = 0
        self.package_code = 0
        self.create_outputs = True
        self.create_zip = True
        self.manifest_mutator = None
        self.zip_mutator = None

    @property
    def outputs(self):
        profile = json.loads(self.profile_path.read_text(encoding="utf-8"))
        return owner._profile_outputs(profile, self.artifacts)

    @property
    def zip_path(self):
        return self.globus / self.args.output

    def write_outputs(self):
        outputs = self.outputs
        outputs["f4_refresh1_comparison_input"].write_bytes(self.comparison.read_bytes())
        for stage, key in owner.STAGES.items():
            outputs[key].write_bytes(encoded({"candidate": stage}))
        rows = {stage: {"basename": outputs[key].name,
                        "raw_sha256": owner.sha256_file(outputs[key]),
                        "non_authoritative": True}
                for stage, key in owner.STAGES.items()}
        manifest = {
            "schema_version": owner.MATERIALIZATION_SCHEMA, "non_authoritative": True,
            "accepted_authority_mutated": False, "production_objects_mutated": False,
            "production_application_performed": False, "method_a_promoted": False,
            "source_head": owner.SOURCE_HEAD, "kinematic_token": owner.KINEMATIC,
            "comparison_input": {
                "source_basename": self.comparison.name,
                "raw_sha256": owner.COMPARISON_SHA256,
                "copied_output_basename": outputs["f4_refresh1_comparison_input"].name,
                "copied_output_raw_sha256": owner.COMPARISON_SHA256,
            },
            "scientific_gate": {
                "f2_scientific_payload_match": True,
                "f3_scientific_payload_match": True,
                "f4_scientific_payload_match": False,
                "first_changed_stage": "F4",
            },
            "candidate_outputs": rows, "complete": True, "errors": [],
        }
        if self.manifest_mutator:
            self.manifest_mutator(manifest)
        outputs["f4_refresh2_materialization_manifest"].write_bytes(encoded(manifest))

    def write_zip(self):
        outputs = self.outputs
        records = {key: {"status": "exists", "json_status": "valid",
                         "archive_path": "global/" + path.name,
                         "sha256": owner.sha256_file(path)}
                   for key, path in outputs.items()}
        manifest = {
            "complete": True, "validation_profile": owner.PROFILE_ID,
            "required_analysis_commit": owner.SOURCE_HEAD,
            "git_head": self.args.bundle_commit,
            "requested_kinematic": owner.KINEMATIC,
            "global_artifacts": records, "errors": [],
        }
        if self.zip_mutator:
            self.zip_mutator(manifest)
        with zipfile.ZipFile(self.zip_path, "w") as archive:
            archive.writestr("manifest.json", encoded(manifest))
            for path in outputs.values():
                archive.write(path, "global/" + path.name)

    def runner(self, command, cwd):
        self.calls.append(list(command))
        stdout, returncode = "", 0
        if command[:3] == ["git", "rev-parse", "--show-toplevel"]:
            stdout = str(self.repo)
        elif command == ["git", "rev-parse", "HEAD"]:
            stdout = self.head
        elif command[:2] == ["git", "status"]:
            stdout = self.dirty
        elif command[:2] == ["git", "merge-base"]:
            returncode = 0 if self.ancestor else 1
        elif command[:3] == ["git", "diff", "--name-only"]:
            stdout = self.changed
        elif command[:2] == [self.args.python, owner.MATERIALIZER]:
            returncode = self.materializer_code
            if returncode == 0 and self.create_outputs:
                self.write_outputs()
        elif command[:2] == ["tcsh", owner.WRAPPER]:
            returncode = self.package_code
            if returncode == 0 and self.create_zip:
                self.write_zip()
        else:
            raise AssertionError(f"unexpected command: {command}")
        return subprocess.CompletedProcess(command, returncode, stdout, "synthetic failure" if returncode else "")

    def invoke(self):
        return owner.run_operation(self.args, repo_root=self.repo,
                                   artifact_root=self.artifacts, globus_root=self.globus,
                                   runner=self.runner, profile_path=self.profile_path)

    def calls_for(self, name):
        return [row for row in self.calls if len(row) > 1 and row[1] == name]


class F4Refresh2ExecutionOwnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.fixture = Fixture(Path(self.temporary.name))
        patcher = mock.patch.object(owner, "COMPARISON_SHA256",
                                    owner.sha256_file(self.fixture.comparison))
        patcher.start()
        self.addCleanup(patcher.stop)

    def test_preflight_and_order_and_exact_immutable_package(self):
        fixture = self.fixture
        events = []
        actual_verify_materialization = owner.verify_materialization
        actual_verify_zip = owner.verify_zip

        def verify_materialization(*args):
            events.append("verify_materialization")
            return actual_verify_materialization(*args)

        def verify_zip(*args):
            events.append("verify_zip")
            return actual_verify_zip(*args)

        actual_runner = fixture.runner

        def runner(command, cwd):
            if command[:2] == [fixture.args.python, owner.MATERIALIZER]:
                events.append("materializer")
            if command[:2] == ["tcsh", owner.WRAPPER]:
                events.append("package")
            return actual_runner(command, cwd)

        with mock.patch.object(fixture, "runner", runner), mock.patch.object(
            owner, "verify_materialization", verify_materialization
        ), mock.patch.object(owner, "verify_zip", verify_zip):
            path, sha = fixture.invoke()
        self.assertEqual(events, ["materializer", "verify_materialization", "package", "verify_zip"])
        self.assertEqual(path, fixture.zip_path)
        self.assertEqual(sha, owner.sha256_file(path))
        materializer = fixture.calls_for(owner.MATERIALIZER)[0]
        self.assertEqual(materializer.count("--f1"), 5)
        self.assertNotIn("--overwrite", materializer)
        self.assertEqual(materializer[materializer.index("--source-head") + 1], owner.SOURCE_HEAD)
        package = fixture.calls_for(owner.WRAPPER)[0]
        self.assertEqual(package.count("--immutable"), 5)
        for key, output in fixture.outputs.items():
            index = package.index(str(output))
            self.assertEqual(package[index + 1], owner.sha256_file(output), key)
        self.assertEqual(len(set(fixture.outputs.values())), 5)

    def test_main_prints_exact_requested_zip_and_hash(self):
        path = self.fixture.zip_path
        sha = "f" * 64
        output = io.StringIO()
        argv = ["--python", sys.executable, "--bundle-commit", "a" * 40,
                "--f1", "Left-lowe=x", "--accepted-f2", "x", "--accepted-f3", "x",
                "--accepted-f4", "x", "--comparison", "x", "--artifact-dir", "x",
                "--kinematic", owner.KINEMATIC, "--output", "fresh-review.zip"]
        with mock.patch.object(owner, "run_operation", return_value=(path, sha)), contextlib.redirect_stdout(output):
            self.assertEqual(owner.main(argv), 0)
        self.assertEqual(output.getvalue(), f"review ZIP: {path}\nreview ZIP SHA-256: {sha}\n")

    def test_preflight_failures_write_no_outputs(self):
        cases = (
            ("short_commit", lambda f: setattr(f.args, "bundle_commit", "abc"), "bundle_commit"),
            ("wrong_head", lambda f: setattr(f, "head", "b" * 40), "does_not_match_head"),
            ("dirty", lambda f: setattr(f, "dirty", "?? stray.txt"), "dirty_repository"),
            ("ancestor", lambda f: setattr(f, "ancestor", False), "not_ancestor"),
            ("changed_science", lambda f: setattr(f, "changed", "src/cuts/unreviewed.py"), "unexpected_committed_paths"),
            ("changed_materializer", lambda f: setattr(f, "changed", owner.MATERIALIZER), "unexpected_committed_paths"),
            ("changed_testing_helper", lambda f: setattr(f, "changed", "testing/unreviewed_helper.py"), "unexpected_committed_paths"),
            ("changed_tools_helper", lambda f: setattr(f, "changed", "tools/unreviewed_helper.py"), "unexpected_committed_paths"),
            ("kinematic", lambda f: setattr(f.args, "kinematic", "Q3p0W2p32"), "unsupported_kinematic"),
            ("comparison_sha", lambda f: f.comparison.write_bytes(b"wrong"), "reviewed_comparison_sha_mismatch"),
            ("existing_zip", lambda f: f.zip_path.write_bytes(b"stale"), "zip_destination_unavailable_or_exists"),
            ("wrong_root", lambda f: setattr(f.args, "artifact_dir", f.inputs), "canonical_root"),
            ("unknown_alias", lambda f: f.args.f1.__setitem__(0, "Unknown=x"), "invalid_or_duplicate_f1_alias"),
            ("duplicate_alias", lambda f: f.args.f1.__setitem__(1, f.args.f1[0]), "invalid_or_duplicate_f1_alias"),
            ("missing_alias", lambda f: f.args.f1.pop(), "canonical_five_f1_inputs_required"),
            ("missing_f1", lambda f: next(iter(f.f1.values())).unlink(), "required_input_missing"),
            ("missing_accepted", lambda f: f.accepted["f3"].unlink(), "required_input_missing"),
            ("missing_comparison", lambda f: f.comparison.unlink(), "required_input_missing"),
            ("relative_accepted", lambda f: setattr(f.args, "accepted_f2", Path("relative.json")), "input_paths_must_be_absolute"),
        )
        for label, mutate, message in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                fixture = Fixture(Path(temporary))
                mutate(fixture)
                with self.assertRaisesRegex(ValueError, message):
                    fixture.invoke()
                self.assertFalse(fixture.calls_for(owner.MATERIALIZER))
                self.assertFalse(any(path.exists() for path in fixture.outputs.values()))

    def test_profile_and_stale_output_preflight_failures(self):
        for label, mutate, message in (
            ("profile_id", lambda p: p.__setitem__("validation_profile", "wrong"), "wrong_f4_refresh2_profile_identity"),
            ("required_commit", lambda p: p["source_identity"].__setitem__("required_analysis_commit", "b" * 40), "wrong_required_analysis_commit"),
            ("broadened_allowlist", lambda p: p["source_identity"]["allowed_committed_files"].append("src/unreviewed.py"), "profile_source_identity_not_narrow"),
        ):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                fixture = Fixture(Path(temporary))
                profile = json.loads(fixture.profile_path.read_text(encoding="utf-8"))
                mutate(profile)
                fixture.profile_path.write_bytes(encoded(profile))
                with self.assertRaisesRegex(ValueError, message):
                    fixture.invoke()
                self.assertFalse(fixture.calls_for(owner.MATERIALIZER))
        for key in owner.OUTPUT_SUFFIXES:
            with self.subTest(stale=key), tempfile.TemporaryDirectory() as temporary:
                fixture = Fixture(Path(temporary))
                fixture.outputs[key].write_bytes(b"old")
                with self.assertRaisesRegex(ValueError, "stale_or_partial_materialization_output_exists"):
                    fixture.invoke()
                self.assertEqual(fixture.outputs[key].read_bytes(), b"old")
                self.assertFalse(fixture.calls_for(owner.MATERIALIZER))

    def test_materializer_failure_blocks_package_and_stale_retry(self):
        fixture = self.fixture
        fixture.materializer_code = 3
        with self.assertRaisesRegex(ValueError, "materializer_failed"):
            fixture.invoke()
        self.assertFalse(fixture.calls_for(owner.WRAPPER))
        fixture.materializer_code = 0
        fixture.create_outputs = False
        with self.assertRaisesRegex(ValueError, "materialization_output_missing"):
            fixture.invoke()
        self.assertFalse(fixture.calls_for(owner.WRAPPER))
        fixture.outputs["f4_refresh2_candidate_f2"].write_bytes(b"partial")
        with self.assertRaisesRegex(ValueError, "stale_or_partial_materialization_output_exists"):
            fixture.invoke()

    def test_materialization_manifest_failures_block_package(self):
        fixture = self.fixture
        fixture.write_outputs()
        outputs = fixture.outputs
        manifest_path = outputs["f4_refresh2_materialization_manifest"]
        good = manifest_path.read_bytes()
        cases = [
            ("missing_manifest", None, "materialization_output_missing"),
            ("malformed_json", b"{", None),
            ("wrong_schema", lambda p: p.__setitem__("schema_version", "wrong"), "schema_mismatch"),
            ("source_head", lambda p: p.__setitem__("source_head", "b" * 40), "source_head_mismatch"),
            ("kinematic", lambda p: p.__setitem__("kinematic_token", "wrong"), "kinematic_mismatch"),
            ("comparison", lambda p: p["comparison_input"].__setitem__("raw_sha256", "0" * 64), "comparison_provenance_mismatch"),
            ("f2_gate", lambda p: p["scientific_gate"].__setitem__("f2_scientific_payload_match", False), "scientific_gate_mismatch"),
            ("f3_gate", lambda p: p["scientific_gate"].__setitem__("f3_scientific_payload_match", False), "scientific_gate_mismatch"),
            ("f4_gate", lambda p: p["scientific_gate"].__setitem__("f4_scientific_payload_match", True), "scientific_gate_mismatch"),
            ("first_stage", lambda p: p["scientific_gate"].__setitem__("first_changed_stage", "F3"), "scientific_gate_mismatch"),
            ("basename", lambda p: p["candidate_outputs"]["f2"].__setitem__("basename", "wrong.json"), "candidate_hash_or_name_mismatch"),
            ("raw_hash", lambda p: p["candidate_outputs"]["f4"].__setitem__("raw_sha256", "0" * 64), "candidate_hash_or_name_mismatch"),
        ]
        for field in ("non_authoritative", "accepted_authority_mutated", "production_objects_mutated",
                      "production_application_performed", "method_a_promoted", "complete"):
            cases.append((field, lambda p, key=field: p.__setitem__(key, not p[key]), "flag_mismatch"))
        for label, mutation, message in cases:
            with self.subTest(label=label):
                manifest_path.write_bytes(good)
                if mutation is None:
                    manifest_path.unlink()
                elif isinstance(mutation, bytes):
                    manifest_path.write_bytes(mutation)
                else:
                    payload = json.loads(good)
                    mutation(payload)
                    manifest_path.write_bytes(encoded(payload))
                if message:
                    with self.assertRaisesRegex(ValueError, message):
                        owner.verify_materialization(outputs, fixture.comparison)
                else:
                    with self.assertRaises(json.JSONDecodeError):
                        owner.verify_materialization(outputs, fixture.comparison)
        manifest_path.write_bytes(good)
        outputs["f4_refresh2_candidate_f3"].write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError, "candidate_hash_or_name_mismatch"):
            owner.verify_materialization(outputs, fixture.comparison)

    def test_package_and_zip_failures_never_report_success(self):
        fixture = self.fixture
        fixture.package_code = 2
        with self.assertRaisesRegex(ValueError, "package_wrapper_failed"):
            fixture.invoke()
        self.assertFalse(fixture.zip_path.exists())
        with self.assertRaisesRegex(ValueError, "stale_or_partial_materialization_output_exists"):
            fixture.invoke()

    def test_verification_failure_never_invokes_package(self):
        fixture = self.fixture
        fixture.manifest_mutator = lambda payload: payload["scientific_gate"].__setitem__(
            "first_changed_stage", "F3"
        )
        with self.assertRaisesRegex(ValueError, "materialization_scientific_gate_mismatch"):
            fixture.invoke()
        self.assertFalse(fixture.calls_for(owner.WRAPPER))

    def test_missing_returned_zip_never_reports_success(self):
        fixture = self.fixture
        fixture.create_zip = False
        with self.assertRaisesRegex(ValueError, "returned_zip_missing"):
            fixture.invoke()
        self.assertEqual(len(fixture.calls_for(owner.WRAPPER)), 1)

    def test_returned_zip_negative_cases(self):
        fixture = self.fixture
        fixture.write_outputs()
        hashes = owner.verify_materialization(fixture.outputs, fixture.comparison)
        fixture.write_zip()
        good = fixture.zip_path.read_bytes()
        cases = (
            ("missing_zip", lambda: fixture.zip_path.unlink(), ValueError, "returned_zip_missing"),
            ("bad_zip", lambda: fixture.zip_path.write_bytes(b"bad"), zipfile.BadZipFile, None),
        )
        for label, mutation, kind, message in cases:
            with self.subTest(label=label):
                fixture.zip_path.write_bytes(good)
                mutation()
                with self.assertRaises(kind if message is None else ValueError) as caught:
                    owner.verify_zip(fixture.zip_path, fixture.outputs, hashes, fixture.args.bundle_commit)
                if message:
                    self.assertIn(message, str(caught.exception))
        for label, mutate, message in (
            ("incomplete", lambda p: p.__setitem__("complete", False), "collector_bundle_incomplete"),
            ("profile", lambda p: p.__setitem__("validation_profile", "wrong"), "collector_bundle_identity_mismatch"),
            ("kinematic", lambda p: p.__setitem__("requested_kinematic", "wrong"), "collector_bundle_identity_mismatch"),
            ("missing_global", lambda p: p["global_artifacts"].pop("f4_refresh2_candidate_f2"), "inventory_mismatch"),
            ("sha", lambda p: p["global_artifacts"]["f4_refresh2_candidate_f3"].__setitem__("sha256", "0" * 64), "artifact_mismatch"),
        ):
            with self.subTest(label=label):
                fixture.zip_mutator = mutate
                fixture.write_zip()
                with self.assertRaisesRegex(ValueError, message):
                    owner.verify_zip(fixture.zip_path, fixture.outputs, hashes, fixture.args.bundle_commit)
        fixture.zip_mutator = None
        fixture.write_zip()
        with zipfile.ZipFile(fixture.zip_path, "w") as archive:
            archive.writestr("other.json", "{}")
        with self.assertRaisesRegex(ValueError, "returned_zip_manifest_missing"):
            owner.verify_zip(fixture.zip_path, fixture.outputs, hashes, fixture.args.bundle_commit)
        for raw in (b"{", b'{"complete":NaN}', b'{"complete":true,"complete":true}'):
            with self.subTest(malformed_manifest=raw):
                with zipfile.ZipFile(fixture.zip_path, "w") as archive:
                    archive.writestr("manifest.json", raw)
                with self.assertRaises((ValueError, json.JSONDecodeError)):
                    owner.verify_zip(fixture.zip_path, fixture.outputs, hashes, fixture.args.bundle_commit)

    def test_profile_and_owner_accept_exact_post_hardening_committed_range(self):
        allowed = {
            "testing/pion_hgcer_validation_bundle_profile_f4_refresh2.json",
            "testing/test_pion_hgcer_validation_bundle_profile_f4_refresh2.py",
            "testing/run_f4_refresh2_materialize_verify_package.py",
            "testing/test_run_f4_refresh2_materialize_verify_package.py",
            "testing/test_memory_health.py",
            "tools/check_memory_health.py",
        }
        profile = json.loads((owner.REPO_ROOT / owner.PROFILE_RELATIVE).read_text(encoding="utf-8"))
        self.assertEqual(set(profile["source_identity"]["allowed_committed_files"]), allowed)
        self.assertEqual(owner.ALLOWED_COMMITTED, allowed)
        self.assertEqual(profile["source_identity"]["required_analysis_commit"], owner.SOURCE_HEAD)
        self.fixture.changed = "\n".join(sorted(allowed | {"docs/memory/CURRENT.md"}))
        owner.preflight(self.fixture.args, repo_root=self.fixture.repo,
                        artifact_root=self.fixture.artifacts, globus_root=self.fixture.globus,
                        runner=self.fixture.runner, profile_path=self.fixture.profile_path)
        self.assertFalse(self.fixture.calls_for(owner.MATERIALIZER))


if __name__ == "__main__":
    unittest.main()
