"""Focused, dependency-free tests for the Phase-E.7.2 validation collector."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile


REPO_ROOT = Path(__file__).resolve().parents[1]
COLLECTOR_PATH = REPO_ROOT / "testing" / "collect_pion_hgcer_validation_bundle.py"
PROFILE_PATH = REPO_ROOT / "testing" / "pion_hgcer_validation_bundle_profile.json"
SPEC = importlib.util.spec_from_file_location("_validation_bundle_collector", COLLECTOR_PATH)
collector = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = collector
SPEC.loader.exec_module(collector)


class _FakePdfReader:
    page_count = 12
    fail = False

    def __init__(self, _path):
        if self.fail:
            raise RuntimeError("fake reader failure")
        self.pages = list(range(1, self.page_count + 1))


class _FakePdfWriter:
    def __init__(self):
        self.pages = []

    def add_page(self, page):
        self.pages.append(page)

    def write(self, target):
        target.write(("fake-pdf-pages=" + repr(self.pages)).encode("ascii"))


class _FakePdfModule:
    __version__ = "fake-1"
    PdfReader = _FakePdfReader
    PdfWriter = _FakePdfWriter


def _python_backend():
    return collector.PdfBackend("pypdf", "pypdf fake-1", "python", module=_FakePdfModule)


def _setting(phi="Left", epsilon="lowe", kinematic="Q4p4W2p74"):
    return {
        "phi_setting": phi,
        "epsilon_filename_token": epsilon,
        "kinematic_token": kinematic,
        "particle_type": "kaon",
    }


def _phase_c(phi="Left", epsilon="lowe"):
    return {
        "schema_version": "pion_hgcer_refinement_checkpoint/v1",
        "setting": _setting(phi, epsilon),
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }


def _phase_d(phi="Left", epsilon="lowe"):
    return {
        "schema_version": "pion_hgcer_phase_d_checkpoint/v1",
        "setting": _setting(phi, epsilon),
        "non_authoritative": True,
        "decision_performed": False,
        "statistical_compatibility_claimed": False,
        "production_objects_mutated": False,
        "refinement_applied": False,
    }


def _correction(phi="Left", epsilon="lowe"):
    correction = {
        "schema_version": "pion_hgcer_parent_preserving_correction/v1",
        "status": "available",
        "available": True,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "production_application_performed": False,
        "event_application_performed": False,
        "fingerprint": "correction-fingerprint",
        "t_edges": [0.0, 1.0],
        "delta_edges": [-10.0, 10.0],
        "parents": [{
            "t_index": 0,
            "status": "identity_no_refinable_cells",
            "closure_passed": True,
            "refinable_cell_count": 0,
        }],
        "cells": [{"t_index": 0, "delta_index": 0}],
    }
    return {
        "schema_version": "pion_hgcer_parent_preserving_correction_artifact/v1",
        "setting": _setting(phi, epsilon),
        "correction": correction,
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
        "production_application_performed": False,
        "event_application_performed": False,
    }


def _sidecar(phi="Left", epsilon="lowe", page_count=12, failures=None):
    pages = [
        {"page_id": "full_background.fixture.{}".format(index), "scope": "t1", "authoritative": False}
        for index in range(max(0, page_count - 4))
    ]
    pages.extend(
        {"page_id": page_id, "scope": "setting", "authoritative": False}
        for page_id in collector._E72_PAGE_IDS
    )
    return {
        "schema_version": "full_background_subtraction_page_manifest/v1",
        "setting": _setting(phi, epsilon),
        "pdf_basename": collector.full_background_subtraction_basename(phi, "Q4p4W2p74", epsilon),
        "pages": pages,
        "renderer_failures": [] if failures is None else failures,
    }


def _clean_command_runner(command, _cwd):
    command = [str(item) for item in command]
    stdout = ""
    if command == ["git", "rev-parse", "HEAD"]:
        stdout = "test-head\n"
    elif command == ["git", "rev-parse", "HEAD^"]:
        stdout = "test-parent\n"
    elif command[-1:] == ["--version"]:
        stdout = "Python fake\n"
    return {"command": command, "returncode": 0, "stdout": stdout, "stderr": ""}


class PionHGCerValidationBundleCollectorTests(unittest.TestCase):
    def setUp(self):
        _FakePdfReader.page_count = 12
        _FakePdfReader.fail = False

    def _write_setting(self, outdir, *, phi="Left", epsilon="lowe", phase_c=None,
                       phase_d=None, correction=None, sidecar=None, pdf=True):
        outdir = Path(outdir)
        artifacts = {
            collector.checkpoint_basename(phi, "Q4p4W2p74", epsilon): _phase_c(phi, epsilon) if phase_c is None else phase_c,
            collector.phase_d_checkpoint_basename(phi, "Q4p4W2p74", epsilon): _phase_d(phi, epsilon) if phase_d is None else phase_d,
            collector.parent_preserving_correction_basename(phi, "Q4p4W2p74", epsilon): _correction(phi, epsilon) if correction is None else correction,
            collector.full_background_page_manifest_basename(phi, "Q4p4W2p74", epsilon): _sidecar(phi, epsilon, _FakePdfReader.page_count) if sidecar is None else sidecar,
        }
        for basename, payload in artifacts.items():
            (outdir / basename).write_bytes(json.dumps(payload, sort_keys=True, allow_nan=False).encode("utf-8"))
        pdf_path = outdir / collector.full_background_subtraction_basename(phi, "Q4p4W2p74", epsilon)
        if pdf:
            pdf_path.write_bytes(b"original-full-background-pdf")
        return artifacts, pdf_path

    def _collect(self, temporary, *, backends=None, command_runner=_clean_command_runner,
                 write=True, **kwargs):
        source = Path(temporary) / "source"
        source.mkdir()
        if write:
            self._write_setting(source, **kwargs)
        output = Path(temporary) / "bundle.zip"
        result = collector.collect_validation_bundle(
            outdir=source, kinematic="Q4p4W2p74", output=output,
            phi="Left", epsilon="lowe", repo_root=REPO_ROOT,
            pdf_backends=[_python_backend()] if backends is None else backends,
            command_runner=command_runner,
        )
        return result, output, source

    def test_profile_v2_and_five_setting_selection(self):
        profile = collector.load_validation_profile(PROFILE_PATH)
        self.assertEqual(profile["schema_version"], collector.PROFILE_SCHEMA_VERSION)
        self.assertEqual(profile["validation_profile"], "phase_e7_1_batched_farm_meeting_review/v1")
        self.assertEqual(collector.resolve_settings(), (
            ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
            ("Center", "highe"), ("Right", "highe"),
        ))
        with self.assertRaisesRegex(ValueError, "setting_not_authorized"):
            collector.resolve_settings("Right", "lowe")
        with self.assertRaisesRegex(ValueError, "supplied_together"):
            collector.resolve_settings("Left", None)

    def test_deterministic_artifact_names_and_final_four_selection(self):
        self.assertEqual(
            collector.phase_d_checkpoint_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_pion-background_hgcer_ab_comparison_Q4p4W2p74_lowe.json",
        )
        self.assertEqual(
            collector.parent_preserving_correction_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_pion-background_hgcer_parent-preserving-correction_Q4p4W2p74_lowe.json",
        )
        self.assertEqual(
            collector.full_background_page_manifest_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json",
        )
        self.assertEqual(
            collector.meeting_summary_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_rand_sub_Q4p4W2p74_lowe_meeting-summary.pdf",
        )
        self.assertEqual(collector.select_validation_pages(112), [109, 110, 111, 112])
        self.assertEqual(collector.select_validation_pages(4), [1, 2, 3, 4])
        with self.assertRaisesRegex(ValueError, "too_short_for_meeting_summary"):
            collector.select_validation_pages(3)

    def test_complete_bundle_archives_all_json_and_only_slim_pdf(self):
        with tempfile.TemporaryDirectory() as temporary:
            result, output, source = self._collect(temporary)
            self.assertEqual(result["returncode"], 0)
            self.assertEqual(
                collector.sha256_file(source / collector.full_background_subtraction_basename("Left", "Q4p4W2p74", "lowe")),
                hashlib.sha256(b"original-full-background-pdf").hexdigest(),
            )
            with zipfile.ZipFile(output) as archive:
                names = set(archive.namelist())
                self.assertIn("Left_lowe/" + collector.checkpoint_basename("Left", "Q4p4W2p74", "lowe"), names)
                self.assertIn("Left_lowe/" + collector.phase_d_checkpoint_basename("Left", "Q4p4W2p74", "lowe"), names)
                self.assertIn("Left_lowe/" + collector.parent_preserving_correction_basename("Left", "Q4p4W2p74", "lowe"), names)
                self.assertIn("Left_lowe/" + collector.full_background_page_manifest_basename("Left", "Q4p4W2p74", "lowe"), names)
                self.assertIn("Left_lowe/" + collector.meeting_summary_basename("Left", "Q4p4W2p74", "lowe"), names)
                self.assertNotIn("Left_lowe/" + collector.full_background_subtraction_basename("Left", "Q4p4W2p74", "lowe"), names)
                manifest = json.loads(archive.read("manifest.json"))
            artifacts = manifest["settings"][0]["artifacts"]
            self.assertEqual(artifacts["full_background_subtraction_pdf"]["original_page_count"], 12)
            self.assertEqual(artifacts["full_background_subtraction_pdf"]["extracted_pages"], [9, 10, 11, 12])
            self.assertEqual(artifacts["meeting_summary_pdf"]["backend_identity"], "pypdf fake-1")

    def test_strict_setting_and_json_failures_remain_best_effort(self):
        with tempfile.TemporaryDirectory() as temporary:
            result, output, source = self._collect(temporary)
            correction_path = source / collector.parent_preserving_correction_basename("Left", "Q4p4W2p74", "lowe")
            correction_path.write_text("{\"value\": NaN}", encoding="utf-8")
            output.unlink()
            result = collector.collect_validation_bundle(
                outdir=source, kinematic="Q4p4W2p74", output=output,
                phi="Left", epsilon="lowe", repo_root=REPO_ROOT,
                pdf_backends=[_python_backend()], command_runner=_clean_command_runner,
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            self.assertIn("checkpoint_json_invalid", {entry["code"] for entry in result["manifest"]["errors"]})

    def test_page_sidecar_is_the_extraction_gate(self):
        with tempfile.TemporaryDirectory() as temporary:
            sidecar = _sidecar(failures=["E.7.2: renderer unavailable"])
            result, output, _source = self._collect(temporary, sidecar=sidecar)
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            artifacts = result["manifest"]["settings"][0]["artifacts"]
            self.assertEqual(artifacts["meeting_summary_pdf"]["status"], "unavailable")
            self.assertEqual(artifacts["full_background_page_manifest"]["metadata"]["metadata_status"], "mismatch")

    def test_pdf_page_count_mismatch_and_backend_failure_are_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            sidecar = _sidecar(page_count=11)
            result, output, _source = self._collect(temporary, sidecar=sidecar)
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            self.assertIn("checkpoint_metadata_mismatch", {entry["code"] for entry in result["manifest"]["errors"]})
        with tempfile.TemporaryDirectory() as temporary:
            result, output, _source = self._collect(temporary, backends=[])
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            self.assertIn("pdf_extraction_dependency_unavailable", {entry["code"] for entry in result["manifest"]["errors"]})

    def test_phase_d_and_correction_authority_failures_are_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            bad_phase_d = _phase_d()
            bad_phase_d["decision_performed"] = True
            result, output, _source = self._collect(temporary, phase_d=bad_phase_d)
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
        with tempfile.TemporaryDirectory() as temporary:
            bad_correction = _correction()
            bad_correction["correction"]["parents"][0]["closure_passed"] = False
            result, output, _source = self._collect(temporary, correction=bad_correction)
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())

    def test_identity_parent_is_accepted_and_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temporary:
            result, output, _source = self._collect(temporary)
            self.assertEqual(result["returncode"], 0)
            with self.assertRaisesRegex(ValueError, "output_path_already_exists"):
                collector.collect_validation_bundle(
                    outdir=Path(temporary) / "source", kinematic="Q4p4W2p74", output=output,
                    phi="Left", epsilon="lowe", repo_root=REPO_ROOT,
                    pdf_backends=[_python_backend()], command_runner=_clean_command_runner,
                )

    def test_source_identity_profile_and_collector_are_runtime_detached(self):
        profile = collector.load_validation_profile()
        self.assertEqual(profile["source_identity"]["required_analysis_commit"], "5795b73bb97c73379948efe6c03ea8ee2a5cabd4")
        self.assertEqual(set(profile["source_identity"]["allowed_committed_files"]), {
            "src/cuts/full_background_subtraction_plots.py", "src/cuts/rand_sub.py",
            "testing/test_full_background_subtraction_plots.py",
            "testing/collect_pion_hgcer_validation_bundle.py",
            "testing/test_collect_pion_hgcer_validation_bundle.py",
            "testing/pion_hgcer_validation_bundle_profile.json",
        })
        source = COLLECTOR_PATH.read_text(encoding="utf-8")
        self.assertNotIn("import rand_sub", source)
        self.assertNotIn("build_pion_hgcer", source)
        self.assertIn("git_diff_check_required_analysis_commit_range", source)


if __name__ == "__main__":
    unittest.main()
