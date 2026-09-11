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
        "epsilon_setting": epsilon,
        "Q2": 4.4,
        "W": 2.74,
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
            "t_low": 0.0,
            "t_high": 1.0,
            "parent_status": "identity_no_refinable_cells",
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


def _acceptance(phi="Left", epsilon="lowe", *, application_outside_delta=False):
    feature_metadata = json.loads(json.dumps(collector._F1_FEATURE_METADATA))
    training_records = [
        {
            "source_label": "prompt", "entry_index": 1,
            "coordinate_fingerprint": "coordinates", "t_index": 0,
            "t_low": 0.0, "t_high": 1.0, "SHMS_delta": 0.0,
            "delta_index": 0, "delta_low": -10.0, "delta_high": 10.0,
            "P_hgcer_npeSum": 1.0, "response_class": "low",
            "P_hgcer_xAtCer": 0.1, "P_hgcer_yAtCer": -0.1,
            "SHMS_xptar": 0.01, "SHMS_yptar": -0.02,
            "allcuts": True, "nommcuts": True,
            "analysis_t": 0.4, "analysis_MM": 1.1,
            "diagnostic_weight": 1.0, "Q2": None, "W": None,
            "epsilon": None, "phi": None,
        },
        {
            "source_label": "prompt", "entry_index": 2,
            "coordinate_fingerprint": "coordinates", "t_index": 0,
            "t_low": 0.0, "t_high": 1.0, "SHMS_delta": 0.1,
            "delta_index": 0, "delta_low": -10.0, "delta_high": 10.0,
            "P_hgcer_npeSum": 3.0, "response_class": "control",
            "P_hgcer_xAtCer": 0.2, "P_hgcer_yAtCer": -0.2,
            "SHMS_xptar": 0.02, "SHMS_yptar": -0.01,
            "allcuts": True, "nommcuts": True,
            "analysis_t": 0.5, "analysis_MM": 1.2,
            "diagnostic_weight": 1.0, "Q2": None, "W": None,
            "epsilon": None, "phi": None,
        },
    ]
    application_records = [{
        "source_label": "prompt", "entry_index": 4,
        "coordinate_fingerprint": "coordinates", "t_index": 0,
        "t_low": 0.0, "t_high": 1.0,
        "phi_degrees": 0.0, "phi_index": 0,
        "phi_low": -5.0, "phi_high": 5.0, "phi_status": "inside_phi",
        "SHMS_delta": 0.0, "delta_index": 0,
        "delta_low": -10.0, "delta_high": 10.0,
        "P_hgcer_npeSum": 3.0, "P_hgcer_xAtCer": 0.1,
        "P_hgcer_yAtCer": -0.1, "SHMS_xptar": 0.01,
        "SHMS_yptar": -0.02, "HMS_xptar": None, "HMS_yptar": None,
        "signed_source_coefficient": 1.0, "baseline_pion_weight_w0": 2.0,
        "signed_baseline_event_contribution": 2.0,
        "allcuts": True, "nommcuts": True,
        "analysis_t": 0.4, "analysis_MM": 1.1,
        "Q2": None, "W": None, "epsilon": None, "theta_cm_deg": None,
    }]
    if application_outside_delta:
        application_records[0].update(
            delta_index=None, delta_low=None, delta_high=None,
        )
    audit = {
        "observed_nonpositive_response_record_count": 0,
        "observed_prompt_nommcuts_nonpositive_response_count": 0,
        "zero_or_nonpositive_included_in_training": False,
        "absolute_leakage_probability_claimed": False,
    }
    training_summary = collector._f1_training_summary(
        training_records, [0.0, 1.0], [-10.0, 10.0], audit,
    )
    application_summary = collector._f1_application_summary(
        application_records, 1, 1, [], [],
    )
    child_projection = [{
        "source_label": application_records[0]["source_label"],
        "entry_index": application_records[0]["entry_index"],
        "t_index": application_records[0]["t_index"],
        "phi_index": application_records[0]["phi_index"],
        "phi_low": application_records[0]["phi_low"],
        "phi_high": application_records[0]["phi_high"],
        "phi_status": application_records[0]["phi_status"],
    }]
    inputs = {
        "schema_version": "pion_hgcer_method_a_acceptance_event_contract/v2",
        "fingerprint_schema_version": "pion_hgcer_method_a_acceptance_event_contract_fingerprint/v2",
        "phase_a_contract_fingerprint": "phase-a",
        "phase_a_pion_event_population_fingerprint": "phase-a-pions",
        "method_a_fingerprint": "method-a",
        "method_a_event_population_fingerprint": "method-a-events",
        "part1_config_fingerprint": "part1-config",
        "coordinate_fingerprint": "coordinates", "host_state": "approved",
        "source_target_state": "post_proton_noRF",
        "t_edges": [0.0, 1.0], "delta_edges": [-10.0, 10.0],
        "phi_edges": [-5.0, 5.0],
        "method_a_training_population_fingerprint": collector._canonical_sha256(training_records),
        "application_population_fingerprint": collector._canonical_sha256(application_records),
        "acceptance_feature_metadata_fingerprint": collector._canonical_sha256(feature_metadata),
        "application_child_assignment_projection_fingerprint": collector._canonical_sha256(child_projection),
        "method_a_closure": training_summary["by_t_delta"],
        "feature_metadata": feature_metadata,
    }
    contract = {
        "schema_version": "pion_hgcer_method_a_acceptance_event_contract/v2",
        "fingerprint_schema_version": "pion_hgcer_method_a_acceptance_event_contract_fingerprint/v2",
        "status": "available", "available": True, "reason": None,
        "diagnostic_stage": "complete",
        "non_authoritative": True, "production_objects_mutated": False,
        "refinement_applied": False, "production_application_performed": False,
        "event_application_performed": False, "method_b_numerical_dependency": False,
        "future_weight_adjustment_constructed": False,
        "phase_a_contract_fingerprint": "phase-a",
        "phase_a_pion_event_population_fingerprint": "phase-a-pions",
        "method_a_fingerprint": "method-a",
        "method_a_event_population_fingerprint": "method-a-events",
        "part1_config_fingerprint": "part1-config",
        "coordinate_fingerprint": "coordinates", "host_state": "approved",
        "source_target_state": "post_proton_noRF",
        "t_edges": [0.0, 1.0], "delta_edges": [-10.0, 10.0],
        "phi_edges": [-5.0, 5.0], "feature_metadata": feature_metadata,
        "method_a_training_records": training_records,
        "application_records": application_records,
        "method_a_training_summary": training_summary,
        "application_summary": application_summary,
        "unmatched_parent_cache_identities": [],
        "unmatched_child_cache_identities": [],
        "method_a_training_population_fingerprint": inputs["method_a_training_population_fingerprint"],
        "application_population_fingerprint": inputs["application_population_fingerprint"],
        "acceptance_feature_metadata_fingerprint": inputs["acceptance_feature_metadata_fingerprint"],
        "application_child_assignment_projection_fingerprint": inputs["application_child_assignment_projection_fingerprint"],
        "fingerprint_inputs": inputs,
        "fingerprint": collector._canonical_sha256(inputs),
    }
    return {
        "schema_version": "pion_hgcer_method_a_acceptance_event_contract_artifact/v2",
        "setting": _setting(phi, epsilon), "contract": contract,
        "non_authoritative": True, "production_objects_mutated": False,
        "refinement_applied": False, "production_application_performed": False,
        "event_application_performed": False,
    }


def _sidecar(phi="Left", epsilon="lowe", page_count=12, failures=None):
    pages = [
        {"page_id": "full_background.fixture.{}".format(index), "scope": "t1", "authoritative": False}
        for index in range(max(0, page_count - 9))
    ]
    pages.extend(
        {"page_id": page_id, "scope": "setting", "authoritative": False}
        for page_id in collector._F1_PAGE_IDS
    )
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
                       phase_d=None, correction=None, acceptance=None, sidecar=None, pdf=True):
        outdir = Path(outdir)
        artifacts = {
            collector.checkpoint_basename(phi, "Q4p4W2p74", epsilon): _phase_c(phi, epsilon) if phase_c is None else phase_c,
            collector.phase_d_checkpoint_basename(phi, "Q4p4W2p74", epsilon): _phase_d(phi, epsilon) if phase_d is None else phase_d,
            collector.parent_preserving_correction_basename(phi, "Q4p4W2p74", epsilon): _correction(phi, epsilon) if correction is None else correction,
            collector.method_a_acceptance_contract_basename(phi, "Q4p4W2p74", epsilon): _acceptance(phi, epsilon) if acceptance is None else acceptance,
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

    def test_profile_v4_and_five_setting_selection(self):
        profile = collector.load_validation_profile(PROFILE_PATH)
        self.assertEqual(profile["schema_version"], collector.PROFILE_SCHEMA_VERSION)
        self.assertEqual(profile["validation_profile"], "phase_f1_batched_farm_acceptance_review/v2")
        self.assertEqual(collector.resolve_settings(), (
            ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
            ("Center", "highe"), ("Right", "highe"),
        ))
        with self.assertRaisesRegex(ValueError, "setting_not_authorized"):
            collector.resolve_settings("Right", "lowe")
        with self.assertRaisesRegex(ValueError, "supplied_together"):
            collector.resolve_settings("Left", None)

    def test_deterministic_artifact_names_and_manifest_selected_pages(self):
        self.assertEqual(
            collector.phase_d_checkpoint_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_pion-background_hgcer_ab_comparison_Q4p4W2p74_lowe.json",
        )
        self.assertEqual(
            collector.parent_preserving_correction_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_pion-background_hgcer_parent-preserving-correction_Q4p4W2p74_lowe.json",
        )
        self.assertEqual(
            collector.method_a_acceptance_contract_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_pion-background_hgcer_method-a-acceptance-contract_Q4p4W2p74_lowe.json",
        )
        self.assertEqual(
            collector.full_background_page_manifest_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json",
        )
        self.assertEqual(
            collector.meeting_summary_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_rand_sub_Q4p4W2p74_lowe_meeting-summary.pdf",
        )
        self.assertEqual(
            collector.f1_acceptance_summary_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_rand_sub_Q4p4W2p74_lowe_f1-acceptance-summary.pdf",
        )
        self.assertEqual(collector.select_validation_pages(112), [109, 110, 111, 112])
        self.assertEqual(collector.select_validation_pages(4), [1, 2, 3, 4])
        with self.assertRaisesRegex(ValueError, "too_short_for_requested_summary"):
            collector.select_validation_pages(3)
        self.assertEqual(
            collector.select_validation_pages(12, {"kind": "explicit_pages", "pages": [4, 5, 6, 7, 8]}),
            [4, 5, 6, 7, 8],
        )

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
                self.assertIn("Left_lowe/" + collector.method_a_acceptance_contract_basename("Left", "Q4p4W2p74", "lowe"), names)
                self.assertIn("Left_lowe/" + collector.full_background_page_manifest_basename("Left", "Q4p4W2p74", "lowe"), names)
                self.assertIn("Left_lowe/" + collector.meeting_summary_basename("Left", "Q4p4W2p74", "lowe"), names)
                self.assertIn("Left_lowe/" + collector.f1_acceptance_summary_basename("Left", "Q4p4W2p74", "lowe"), names)
                self.assertNotIn("Left_lowe/" + collector.full_background_subtraction_basename("Left", "Q4p4W2p74", "lowe"), names)
                manifest = json.loads(archive.read("manifest.json"))
            artifacts = manifest["settings"][0]["artifacts"]
            self.assertEqual(artifacts["full_background_subtraction_pdf"]["original_page_count"], 12)
            self.assertEqual(artifacts["full_background_subtraction_pdf"]["extracted_pages"], {
                "meeting_summary_pdf": [9, 10, 11, 12],
                "f1_acceptance_summary_pdf": [4, 5, 6, 7, 8],
            })
            self.assertEqual(artifacts["meeting_summary_pdf"]["backend_identity"], "pypdf fake-1")
            self.assertEqual(artifacts["f1_acceptance_summary_pdf"]["backend_identity"], "pypdf fake-1")
            self.assertEqual(
                artifacts["parent_preserving_correction"]["metadata"]["parent_statuses"],
                ["identity_no_refinable_cells"],
            )
            self.assertEqual(
                artifacts["method_a_acceptance_contract"]["metadata"]["training_summary"]["prompt_control_count"],
                1,
            )
            self.assertEqual(
                artifacts["method_a_acceptance_contract"]["metadata"]["training_summary"]["prompt_low_count"],
                1,
            )
            self.assertEqual(
                artifacts["method_a_acceptance_contract"]["metadata"]["application_summary"]["phase_a_record_count"],
                1,
            )

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

    def test_f1_contract_and_manifest_page_identities_are_strict(self):
        with tempfile.TemporaryDirectory() as temporary:
            bad_acceptance = _acceptance()
            bad_acceptance["contract"]["method_a_training_records"][1]["response_class"] = "low"
            result, output, _source = self._collect(temporary, acceptance=bad_acceptance)
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            metadata = result["manifest"]["settings"][0]["artifacts"]["method_a_acceptance_contract"]["metadata"]
            self.assertEqual(metadata["metadata_status"], "mismatch")
        with tempfile.TemporaryDirectory() as temporary:
            bad_sidecar = _sidecar()
            bad_sidecar["pages"] = [
                page for page in bad_sidecar["pages"]
                if page["page_id"] != collector._F1_PAGE_IDS[2]
            ]
            result, output, _source = self._collect(temporary, sidecar=bad_sidecar)
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            artifacts = result["manifest"]["settings"][0]["artifacts"]
            self.assertEqual(artifacts["f1_acceptance_summary_pdf"]["status"], "unavailable")

    def test_f1_dual_population_fingerprints_and_summaries_are_archived(self):
        with tempfile.TemporaryDirectory() as temporary:
            result, output, _source = self._collect(temporary)
            self.assertEqual(result["returncode"], 0)
            with zipfile.ZipFile(output) as archive:
                manifest = json.loads(archive.read("manifest.json"))
            f1 = manifest["settings"][0]["artifacts"]["method_a_acceptance_contract"]["metadata"]
            self.assertEqual(f1["contract_fingerprint"], _acceptance()["contract"]["fingerprint"])
            self.assertEqual(f1["training_summary"]["by_t_delta"][0]["control_count"], 1)
            self.assertEqual(f1["application_summary"]["by_t_phi"][0]["prompt_record_count"], 1)
            self.assertEqual(
                manifest["f1_aggregate_summary"][0]["application_summary"]["phase_a_record_count"], 1
            )

    def test_f1_application_delta_geometry_preserves_phase_a_outside_lattice(self):
        with tempfile.TemporaryDirectory() as temporary:
            acceptance = _acceptance(application_outside_delta=True)
            application = acceptance["contract"]["application_records"][0]
            self.assertEqual(
                (application["delta_index"], application["delta_low"], application["delta_high"]),
                (None, None, None),
            )
            result, output, _source = self._collect(temporary, acceptance=acceptance)
            self.assertEqual(result["returncode"], 0)
            self.assertTrue(output.exists())
            metadata = result["manifest"]["settings"][0]["artifacts"]["method_a_acceptance_contract"]["metadata"]
            self.assertEqual(metadata["metadata_status"], "match")

        mutations = (
            ("application_none_with_edge", _acceptance(application_outside_delta=True),
             lambda payload: payload["contract"]["application_records"][0].update(delta_low=-10.0)),
            ("application_integer_wrong_edges", _acceptance(),
             lambda payload: payload["contract"]["application_records"][0].update(delta_high=9.0)),
            ("training_none_delta", _acceptance(),
             lambda payload: payload["contract"]["method_a_training_records"][0].update(
                 delta_index=None, delta_low=None, delta_high=None,
             )),
        )
        for label, acceptance, mutation in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                mutation(acceptance)
                result, output, _source = self._collect(temporary, acceptance=acceptance)
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                metadata = result["manifest"]["settings"][0]["artifacts"]["method_a_acceptance_contract"]["metadata"]
                self.assertEqual(metadata["metadata_status"], "mismatch")

    def test_f1_v2_rejects_stale_schemas_flags_and_population_contract_breaks(self):
        mutations = (
            ("v1_artifact", lambda payload: payload.update(
                schema_version="pion_hgcer_method_a_acceptance_event_contract_artifact/v1"
            )),
            ("wrong_contract_schema", lambda payload: payload["contract"].update(
                schema_version="pion_hgcer_method_a_acceptance_event_contract/v1"
            )),
            ("wrong_fingerprint_schema", lambda payload: payload["contract"].update(
                fingerprint_schema_version="pion_hgcer_method_a_acceptance_event_contract_fingerprint/v1"
            )),
            ("production_mutation", lambda payload: payload["contract"].update(
                production_objects_mutated=True
            )),
            ("refinement", lambda payload: payload["contract"].update(refinement_applied=True)),
            ("application", lambda payload: payload["contract"].update(
                production_application_performed=True
            )),
            ("event_application", lambda payload: payload["contract"].update(
                event_application_performed=True
            )),
            ("method_b_dependency", lambda payload: payload["contract"].update(
                method_b_numerical_dependency=True
            )),
            ("future_weight", lambda payload: payload["contract"].update(
                future_weight_adjustment_constructed=True
            )),
            ("training_nonpositive", lambda payload: payload["contract"][
                "method_a_training_records"][0].update(P_hgcer_npeSum=0.0)
            ),
            ("training_class", lambda payload: payload["contract"][
                "method_a_training_records"][0].update(response_class="control")
            ),
            ("application_not_control", lambda payload: payload["contract"][
                "application_records"][0].update(P_hgcer_npeSum=2.0)
            ),
            ("training_summary", lambda payload: payload["contract"][
                "method_a_training_summary"].update(training_record_count=7)
            ),
            ("application_summary", lambda payload: payload["contract"][
                "application_summary"].update(phase_a_record_count=7)
            ),
            ("training_fingerprint", lambda payload: payload["contract"].update(
                method_a_training_population_fingerprint="not-a-fingerprint"
            )),
            ("application_fingerprint", lambda payload: payload["contract"].update(
                application_population_fingerprint="not-a-fingerprint"
            )),
            ("feature_fingerprint", lambda payload: payload["contract"].update(
                acceptance_feature_metadata_fingerprint="not-a-fingerprint"
            )),
            ("child_assignment_fingerprint", lambda payload: payload["contract"].update(
                application_child_assignment_projection_fingerprint="not-a-fingerprint"
            )),
            ("contract_fingerprint", lambda payload: payload["contract"].update(
                fingerprint="not-a-fingerprint"
            )),
            ("setting", lambda payload: payload["setting"].update(phi_setting="Right")),
        )
        for label, mutation in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                bad_acceptance = _acceptance()
                mutation(bad_acceptance)
                result, output, _source = self._collect(temporary, acceptance=bad_acceptance)
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                self.assertEqual(
                    result["manifest"]["settings"][0]["artifacts"]
                    ["method_a_acceptance_contract"]["metadata"]["metadata_status"],
                    "mismatch",
                )

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

    def test_parent_preserving_parent_lattice_is_complete_unique_and_aligned(self):
        def two_parent_correction():
            payload = _correction()
            correction = payload["correction"]
            correction["t_edges"] = [0.0, 1.0, 2.0]
            correction["parents"] = [
                {
                    "t_index": 0, "t_low": 0.0, "t_high": 1.0,
                    "parent_status": "identity_no_refinable_cells",
                    "closure_passed": True, "refinable_cell_count": 0,
                },
                {
                    "t_index": 1, "t_low": 1.0, "t_high": 2.0,
                    "parent_status": "identity_single_refinable_cell",
                    "closure_passed": True, "refinable_cell_count": 1,
                },
            ]
            correction["cells"] = [
                {"t_index": 0, "delta_index": 0},
                {"t_index": 1, "delta_index": 0},
            ]
            return payload

        mutations = (
            (lambda correction: correction["correction"]["parents"].pop(), "missing_parent"),
            (
                lambda correction: correction["correction"]["parents"].__setitem__(
                    1, dict(correction["correction"]["parents"][0])
                ),
                "duplicate_parent_index",
            ),
            (
                lambda correction: correction["correction"]["parents"][1].update(t_high=3.0),
                "parent_edge_mismatch",
            ),
        )
        for mutation, label in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                bad_correction = two_parent_correction()
                mutation(bad_correction)
                result, output, _source = self._collect(temporary, correction=bad_correction)
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                self.assertIn(
                    "checkpoint_metadata_mismatch",
                    {entry["code"] for entry in result["manifest"]["errors"]},
                )

    def test_missing_required_artifacts_remain_explicit_and_best_effort(self):
        for artifact_key, basename in (
            (
                "phase_c_checkpoint",
                collector.checkpoint_basename("Left", "Q4p4W2p74", "lowe"),
            ),
            (
                "full_background_page_manifest",
                collector.full_background_page_manifest_basename("Left", "Q4p4W2p74", "lowe"),
            ),
        ):
            with self.subTest(artifact=artifact_key), tempfile.TemporaryDirectory() as temporary:
                source = Path(temporary) / "source"
                source.mkdir()
                self._write_setting(source)
                (source / basename).unlink()
                output = Path(temporary) / "bundle.zip"
                result = collector.collect_validation_bundle(
                    outdir=source, kinematic="Q4p4W2p74", output=output,
                    phi="Left", epsilon="lowe", repo_root=REPO_ROOT,
                    pdf_backends=[_python_backend()], command_runner=_clean_command_runner,
                )
                self.assertEqual(result["returncode"], 1)
                self.assertTrue(output.exists())
                self.assertIn(
                    "missing_source_artifact",
                    {entry["code"] for entry in result["manifest"]["errors"]},
                )
                self.assertEqual(
                    result["manifest"]["settings"][0]["artifacts"][artifact_key]["status"],
                    "missing",
                )

    def test_unexpected_committed_file_after_required_base_is_rejected(self):
        def unexpected_file_runner(command, cwd):
            result = _clean_command_runner(command, cwd)
            if list(command)[:3] == ["git", "diff", "--name-only"]:
                result = dict(result)
                result["stdout"] = "src/cuts/unexpected_science_change.py\n"
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _source = self._collect(
                temporary, command_runner=unexpected_file_runner,
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            self.assertIn(
                "unexpected_committed_files_after_required_analysis_commit",
                {entry["code"] for entry in result["manifest"]["errors"]},
            )
        def unrelated_test_runner(command, cwd):
            result = _clean_command_runner(command, cwd)
            if list(command)[:3] == ["git", "diff", "--name-only"]:
                result = dict(result)
                result["stdout"] = "testing/test_unrelated_analysis_behavior.py\n"
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _source = self._collect(
                temporary, command_runner=unrelated_test_runner,
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            self.assertIn(
                "unexpected_committed_files_after_required_analysis_commit",
                {entry["code"] for entry in result["manifest"]["errors"]},
            )

        def agents_file_runner(command, cwd):
            result = _clean_command_runner(command, cwd)
            if list(command)[:3] == ["git", "diff", "--name-only"]:
                result = dict(result)
                result["stdout"] = "AGENTS.md\n"
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _source = self._collect(
                temporary, command_runner=agents_file_runner,
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.exists())
            self.assertIn(
                "unexpected_committed_files_after_required_analysis_commit",
                {entry["code"] for entry in result["manifest"]["errors"]},
            )
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
        self.assertEqual(profile["source_identity"]["required_analysis_commit"], "dc4fc6283001739a487ec80068f951b0e388cae6")
        self.assertEqual(set(profile["source_identity"]["allowed_committed_files"]), {
            "testing/collect_pion_hgcer_validation_bundle.py",
            "testing/test_collect_pion_hgcer_validation_bundle.py",
            "testing/pion_hgcer_validation_bundle_profile.json",
        })
        self.assertEqual(profile["source_identity"]["allowed_non_analysis_path_prefixes"], ["docs/memory/"])
        source = COLLECTOR_PATH.read_text(encoding="utf-8")
        self.assertNotIn("import rand_sub", source)
        self.assertNotIn("build_pion_hgcer", source)
        self.assertIn("git_diff_check_required_analysis_commit_range", source)

    def test_reviewed_fix5_ancestry_allows_only_validation_and_memory_followups(self):
        def reviewed_followup_runner(command, cwd):
            result = _clean_command_runner(command, cwd)
            if list(command)[:3] == ["git", "diff", "--name-only"]:
                result = dict(result)
                result["stdout"] = "\n".join((
                    "docs/memory/CURRENT.md",
                    "docs/memory/handoffs/CURRENT_HANDOFF.md",
                    "testing/collect_pion_hgcer_validation_bundle.py",
                    "testing/pion_hgcer_validation_bundle_profile.json",
                    "testing/test_collect_pion_hgcer_validation_bundle.py",
                    "",
                ))
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _source = self._collect(
                temporary, command_runner=reviewed_followup_runner,
            )
            self.assertEqual(result["returncode"], 0)
            self.assertTrue(output.exists())
            self.assertEqual(
                result["manifest"]["unexpected_committed_files_after_required_analysis_commit"], []
            )


if __name__ == "__main__":
    unittest.main()
