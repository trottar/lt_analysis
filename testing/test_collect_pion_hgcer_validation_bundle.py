"""Focused, dependency-free tests for the Phase-E.3 validation bundle collector."""

from __future__ import annotations

import ast
from copy import deepcopy
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile


REPO_ROOT = Path(__file__).resolve().parents[1]
COLLECTOR_PATH = REPO_ROOT / "testing" / "collect_pion_hgcer_validation_bundle.py"
SPEC = importlib.util.spec_from_file_location("_phase_e3_validation_bundle_collector", COLLECTOR_PATH)
collector = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = collector
SPEC.loader.exec_module(collector)


class _FakePdfReader:
    page_count = 112
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
    return collector.PdfBackend(
        name="pypdf",
        identity="pypdf fake-1",
        kind="python",
        module=_FakePdfModule,
    )


def _checkpoint(phi="Left", epsilon="lowe", kinematic="Q4p4W2p74"):
    return {
        "schema_version": "pion_hgcer_refinement_checkpoint/v1",
        "setting": {
            "phi_setting": phi,
            "epsilon_filename_token": epsilon,
            "kinematic_token": kinematic,
            "particle_type": "kaon",
        },
        "non_authoritative": True,
        "production_objects_mutated": False,
        "refinement_applied": False,
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
        _FakePdfReader.page_count = 112
        _FakePdfReader.fail = False

    def _write_setting(self, outdir, *, checkpoint=None, pdf=True):
        outdir = Path(outdir)
        checkpoint_path = outdir / collector.checkpoint_basename("Left", "Q4p4W2p74", "lowe")
        pdf_path = outdir / collector.full_background_subtraction_basename(
            "Left", "Q4p4W2p74", "lowe"
        )
        if checkpoint is not None:
            checkpoint_path.write_bytes(checkpoint)
        if pdf:
            pdf_path.write_bytes(b"original-debug-pdf-bytes")
        return checkpoint_path, pdf_path

    def _collect_left_lowe(self, temporary, *, checkpoint=None, pdf=True, backends=None, command_runner=_clean_command_runner):
        source = Path(temporary) / "source"
        source.mkdir()
        checkpoint_path, pdf_path = self._write_setting(
            source,
            checkpoint=(json.dumps(_checkpoint(), sort_keys=True).encode("utf-8") if checkpoint is None else checkpoint),
            pdf=pdf,
        )
        output = Path(temporary) / "bundle.zip"
        result = collector.collect_validation_bundle(
            outdir=source,
            kinematic="Q4p4W2p74",
            output=output,
            phi="Left",
            epsilon="lowe",
            repo_root=REPO_ROOT,
            pdf_backends=[_python_backend()] if backends is None else backends,
            command_runner=command_runner,
        )
        return result, output, checkpoint_path, pdf_path

    def test_e3_page_selection_is_exact_for_final_three_pages(self):
        self.assertEqual(
            collector.select_validation_pages(60),
            [58, 59, 60],
        )
        self.assertEqual(
            collector.select_validation_pages(45),
            [43, 44, 45],
        )
        self.assertEqual(collector.select_validation_pages(3), [1, 2, 3])
        with self.assertRaisesRegex(ValueError, "pdf_page_count_too_short_for_e3"):
            collector.select_validation_pages(2)
        with self.assertRaisesRegex(ValueError, "pdf_page_count_invalid"):
            collector.select_validation_pages(0)

    def test_e3_fix2_profile_requires_only_explicit_left_lowe(self):
        self.assertEqual(collector.resolve_settings("Left", "lowe"), (("Left", "lowe"),))
        with self.assertRaisesRegex(ValueError, "requires_explicit_left_lowe"):
            collector.resolve_settings()
        with self.assertRaisesRegex(ValueError, "phi_and_epsilon"):
            collector.resolve_settings("Left", None)
        for phi, epsilon in (("Left", "highe"), ("Center", "lowe"), ("Right", "highe")):
            with self.assertRaisesRegex(ValueError, "setting_not_authorized"):
                collector.resolve_settings(phi, epsilon)

    def test_deterministic_names_and_sha256(self):
        self.assertEqual(
            collector.checkpoint_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_pion-background_hgcer_refinement_checkpoint_Q4p4W2p74_lowe.json",
        )
        self.assertEqual(
            collector.full_background_subtraction_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction.pdf",
        )
        self.assertEqual(
            collector.e3_validation_basename("Left", "Q4p4W2p74", "lowe"),
            "Left_kaon_rand_sub_Q4p4W2p74_lowe_E3-validation.pdf",
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bytes.bin"
            path.write_bytes(b"KaonLT")
            self.assertEqual(
                collector.sha256_file(path),
                hashlib.sha256(b"KaonLT").hexdigest(),
            )

    def test_complete_bundle_copies_checkpoint_bytes_and_slim_pdf_only(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            checkpoint_path, pdf_path = self._write_setting(
                source, checkpoint=json.dumps(_checkpoint(), sort_keys=True).encode("utf-8")
            )
            checkpoint_before = checkpoint_path.read_bytes()
            pdf_before = pdf_path.read_bytes()
            output = Path(temporary) / "bundle.zip"
            result = collector.collect_validation_bundle(
                outdir=source,
                kinematic="Q4p4W2p74",
                output=output,
                phi="Left",
                epsilon="lowe",
                repo_root=REPO_ROOT,
                pdf_backends=[_python_backend()],
                command_runner=_clean_command_runner,
            )
            self.assertEqual(result["returncode"], 0)
            self.assertTrue(output.is_file())
            self.assertEqual(checkpoint_path.read_bytes(), checkpoint_before)
            self.assertEqual(pdf_path.read_bytes(), pdf_before)
            with zipfile.ZipFile(output) as archive:
                names = set(archive.namelist())
                checkpoint_name = "Left_lowe/" + checkpoint_path.name
                slim_name = "Left_lowe/" + collector.e3_validation_basename(
                    "Left", "Q4p4W2p74", "lowe"
                )
                self.assertEqual(
                    names,
                    {"Left_lowe/", checkpoint_name, slim_name, "manifest.json", "source_state.txt", "source_checks.txt"},
                )
                self.assertEqual(archive.read(checkpoint_name), checkpoint_before)
                self.assertNotIn(pdf_path.name, names)
                self.assertNotIn("Left_kaon_rand_sub_Q4p4W2p74_lowe_hgcer-debug.pdf", names)
                self.assertNotIn("Left_kaon_rand_sub_Q4p4W2p74_lowe_hgcer-validation.pdf", names)
                manifest = json.loads(archive.read("manifest.json"))
                slim_payload = archive.read(slim_name)
            pdf_artifact = manifest["settings"][0]["artifacts"]["full_background_subtraction_pdf"]
            self.assertEqual(pdf_artifact["original_page_count"], 112)
            self.assertEqual(pdf_artifact["extracted_pages"], collector.select_validation_pages(112))
            self.assertEqual(pdf_artifact["slim_pdf"]["backend"], "pypdf")
            self.assertEqual(pdf_artifact["sha256"], hashlib.sha256(pdf_before).hexdigest())
            self.assertEqual(pdf_artifact["slim_pdf"]["sha256"], hashlib.sha256(slim_payload).hexdigest())
            self.assertEqual(slim_payload, b"fake-pdf-pages=[110, 111, 112]")
            self.assertEqual(manifest["validation_profile"], collector.VALIDATION_PROFILE)
            self.assertEqual(manifest["settings"][0]["kinematic"], "Q4p4W2p74")
            self.assertEqual(manifest["git_head"], "test-head")
            self.assertEqual(manifest["required_analysis_commit"], collector.E3_FIX2_REQUIRED_COMMIT)
            self.assertTrue(manifest["required_analysis_commit_is_ancestor"])
            self.assertEqual(manifest["committed_files_after_required_analysis_commit"], [])
            self.assertTrue(manifest["complete"])

    def test_collector_has_no_analysis_runtime_imports(self):
        source = COLLECTOR_PATH.read_text(encoding="utf-8")
        imports = []
        for node in ast.walk(ast.parse(source)):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
        self.assertTrue(set(imports).issubset({
            "__future__", "argparse", "contextlib", "dataclasses", "datetime", "hashlib",
            "importlib", "io", "json", "os", "pathlib", "re", "shutil",
            "subprocess", "sys", "tempfile", "typing", "zipfile",
        }))

    def test_missing_artifacts_still_write_bundle_and_return_nonzero(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "empty"
            source.mkdir()
            output = Path(temporary) / "bundle.zip"
            result = collector.collect_validation_bundle(
                outdir=source,
                kinematic="Q4p4W2p74",
                output=output,
                phi="Left",
                epsilon="lowe",
                repo_root=REPO_ROOT,
                pdf_backends=[],
                command_runner=_clean_command_runner,
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.is_file())
            self.assertEqual(
                [issue["code"] for issue in result["manifest"]["errors"]],
                ["missing_source_artifact", "missing_source_artifact"],
            )
            with zipfile.ZipFile(output) as archive:
                self.assertIn("Left_lowe/", archive.namelist())
                self.assertIn("manifest.json", archive.namelist())

    def test_missing_checkpoint_or_full_background_pdf_still_writes_e3_zip(self):
        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "missing-checkpoint"
            source.mkdir()
            _checkpoint_path, full_background_path = self._write_setting(
                source, checkpoint=None, pdf=True
            )
            output = Path(temporary) / "missing-checkpoint.zip"
            result = collector.collect_validation_bundle(
                outdir=source, kinematic="Q4p4W2p74", output=output,
                phi="Left", epsilon="lowe", repo_root=REPO_ROOT,
                pdf_backends=[_python_backend()], command_runner=_clean_command_runner,
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.is_file())
            self.assertTrue(full_background_path.is_file())
            self.assertIn(
                "missing_source_artifact",
                [issue["code"] for issue in result["manifest"]["errors"]],
            )
            with zipfile.ZipFile(output) as archive:
                self.assertIn(
                    "Left_lowe/" + collector.e3_validation_basename(
                        "Left", "Q4p4W2p74", "lowe"
                    ),
                    archive.namelist(),
                )

        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "missing-full-background"
            source.mkdir()
            checkpoint_path, _full_background_path = self._write_setting(
                source, checkpoint=json.dumps(_checkpoint()).encode("utf-8"), pdf=False
            )
            output = Path(temporary) / "missing-full-background.zip"
            result = collector.collect_validation_bundle(
                outdir=source, kinematic="Q4p4W2p74", output=output,
                phi="Left", epsilon="lowe", repo_root=REPO_ROOT,
                pdf_backends=[_python_backend()], command_runner=_clean_command_runner,
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.is_file())
            with zipfile.ZipFile(output) as archive:
                self.assertIn("Left_lowe/" + checkpoint_path.name, archive.namelist())
                self.assertNotIn(
                    "Left_lowe/" + collector.e3_validation_basename(
                        "Left", "Q4p4W2p74", "lowe"
                    ),
                    archive.namelist(),
                )
            artifact = result["manifest"]["settings"][0]["artifacts"][
                "full_background_subtraction_pdf"
            ]
            self.assertEqual(artifact["status"], "missing")

    def test_short_full_background_pdf_records_narrow_e3_page_error(self):
        _FakePdfReader.page_count = 2
        with tempfile.TemporaryDirectory() as temporary:
            result, output, _checkpoint_path, _pdf_path = self._collect_left_lowe(
                temporary
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.is_file())
            self.assertIn(
                "pdf_page_count_too_short_for_e3",
                [issue["code"] for issue in result["manifest"]["errors"]],
            )
            artifact = result["manifest"]["settings"][0]["artifacts"][
                "full_background_subtraction_pdf"
            ]
            self.assertEqual(artifact["original_page_count"], 2)
            self.assertEqual(artifact["extracted_pages"], [])
            self.assertEqual(artifact["extraction_error"], "pdf_page_count_too_short_for_e3")
            with zipfile.ZipFile(output) as archive:
                self.assertNotIn(
                    "Left_lowe/" + collector.e3_validation_basename(
                        "Left", "Q4p4W2p74", "lowe"
                    ),
                    archive.namelist(),
                )

    def test_checkpoint_json_and_metadata_failures_are_recorded_without_aborting(self):
        with tempfile.TemporaryDirectory() as temporary:
            invalid_result, invalid_output, _checkpoint_path, _pdf_path = self._collect_left_lowe(
                temporary, checkpoint=b"not-json"
            )
            self.assertEqual(invalid_result["returncode"], 1)
            self.assertTrue(invalid_output.is_file())
            self.assertIn(
                "checkpoint_json_invalid",
                [issue["code"] for issue in invalid_result["manifest"]["errors"]],
            )
        with tempfile.TemporaryDirectory() as temporary:
            mismatched = deepcopy(_checkpoint())
            mismatched["setting"]["phi_setting"] = "Center"
            result, output, _checkpoint_path, _pdf_path = self._collect_left_lowe(
                temporary, checkpoint=json.dumps(mismatched).encode("utf-8")
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.is_file())
            metadata = result["manifest"]["settings"][0]["artifacts"]["phase_c_checkpoint"]["checkpoint_metadata"]
            self.assertEqual(metadata["metadata_status"], "mismatch")
            self.assertIn("phi_setting", metadata["mismatches"])
            self.assertIn(
                "checkpoint_metadata_mismatch",
                [issue["code"] for issue in result["manifest"]["errors"]],
            )

    def test_unavailable_extraction_dependency_and_source_check_failure_are_nonfatal(self):
        def failing_check_runner(command, cwd):
            result = _clean_command_runner(command, cwd)
            if "py_compile" in result["command"]:
                result["returncode"] = 7
                result["stderr"] = "synthetic compile failure\n"
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _checkpoint_path, _pdf_path = self._collect_left_lowe(
                temporary, backends=[], command_runner=failing_check_runner
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.is_file())
            codes = [issue["code"] for issue in result["manifest"]["errors"]]
            self.assertIn("pdf_extraction_dependency_unavailable", codes)
            self.assertIn("source_check_failed", codes)
            with zipfile.ZipFile(output) as archive:
                source_checks = archive.read("source_checks.txt").decode("utf-8")
            self.assertIn("synthetic compile failure", source_checks)
            self.assertIn(
                "src/cuts/full_background_subtraction_plots.py", source_checks
            )
            self.assertIn(
                "testing/test_full_background_subtraction_plots.py", source_checks
            )
            self.assertIn(
                "testing/test_collect_pion_hgcer_validation_bundle.py", source_checks
            )

    def test_e3_source_check_profile_and_existing_output_protection(self):
        source_state_text, source_head = collector.collect_source_state(
            REPO_ROOT, _clean_command_runner
        )
        self.assertEqual(source_head, "test-head")
        self.assertIn("git rev-parse HEAD^", source_state_text)
        self.assertIn("git status --short", source_state_text)
        self.assertIn("Python fake", source_state_text)
        source_checks_text, source_checks = collector.collect_source_checks(
            REPO_ROOT, _clean_command_runner
        )
        self.assertEqual(
            [record["name"] for record in source_checks],
            [
                "required_analysis_commit_ancestor",
                "committed_files_after_required_analysis_commit",
                "py_compile",
                "full_background_subtraction_plots_unittest",
                "hgcer_refinement_plots_unittest",
                "validation_bundle_collector_unittest",
                "git_diff_check",
            ],
        )
        self.assertIn("testing.test_full_background_subtraction_plots", source_checks_text)
        self.assertIn("testing.test_pion_hgcer_refinement_plots", source_checks_text)
        self.assertIn("testing.test_collect_pion_hgcer_validation_bundle", source_checks_text)

        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source"
            source.mkdir()
            self._write_setting(
                source, checkpoint=json.dumps(_checkpoint()).encode("utf-8")
            )
            output = Path(temporary) / "existing.zip"
            output.write_bytes(b"do-not-overwrite")
            with self.assertRaisesRegex(ValueError, "output_path_already_exists"):
                collector.collect_validation_bundle(
                    outdir=source, kinematic="Q4p4W2p74", output=output,
                    phi="Left", epsilon="lowe", repo_root=REPO_ROOT,
                    pdf_backends=[_python_backend()], command_runner=_clean_command_runner,
                )
            self.assertEqual(output.read_bytes(), b"do-not-overwrite")

    def test_required_analysis_commit_gate_is_best_effort_for_stale_or_extra_commits(self):
        def stale_runner(command, cwd):
            result = _clean_command_runner(command, cwd)
            command = result["command"]
            if command == ["git", "rev-parse", "HEAD"]:
                result["stdout"] = "85db000000000000000000000000000000000000\n"
            elif command[:4] == ["git", "merge-base", "--is-ancestor", collector.E3_FIX2_REQUIRED_COMMIT]:
                result["returncode"] = 1
                result["stderr"] = "required commit absent\n"
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _checkpoint_path, _pdf_path = self._collect_left_lowe(
                temporary, command_runner=stale_runner
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.is_file())
            self.assertFalse(result["manifest"]["required_analysis_commit_is_ancestor"])
            self.assertIn(
                "required_analysis_commit_not_present",
                [issue["code"] for issue in result["manifest"]["errors"]],
            )

        def allowed_file_runner(command, cwd):
            result = _clean_command_runner(command, cwd)
            if result["command"][:3] == ["git", "diff", "--name-only"]:
                result["stdout"] = "testing/collect_pion_hgcer_validation_bundle.py\n"
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _checkpoint_path, _pdf_path = self._collect_left_lowe(
                temporary, command_runner=allowed_file_runner
            )
            self.assertEqual(result["returncode"], 0)
            self.assertTrue(output.is_file())
            self.assertEqual(
                result["manifest"]["committed_files_after_required_analysis_commit"],
                ["testing/collect_pion_hgcer_validation_bundle.py"],
            )

        def extra_file_runner(command, cwd):
            result = _clean_command_runner(command, cwd)
            if result["command"][:3] == ["git", "diff", "--name-only"]:
                result["stdout"] = "src/cuts/full_background_subtraction_plots.py\n"
            return result

        with tempfile.TemporaryDirectory() as temporary:
            result, output, _checkpoint_path, _pdf_path = self._collect_left_lowe(
                temporary, command_runner=extra_file_runner
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.is_file())
            self.assertEqual(
                result["manifest"]["unexpected_committed_files_after_required_analysis_commit"],
                ["src/cuts/full_background_subtraction_plots.py"],
            )
            self.assertIn(
                "unexpected_committed_files_after_required_analysis_commit",
                [issue["code"] for issue in result["manifest"]["errors"]],
            )

    def test_extraction_failure_is_recorded_without_discarding_the_bundle(self):
        class _AlwaysFailingPdfModule:
            PdfReader = type(
                "AlwaysFailingReader",
                (),
                {
                    "__init__": lambda self, _path: (_ for _ in ()).throw(
                        RuntimeError("synthetic extraction failure")
                    )
                },
            )
            PdfWriter = _FakePdfWriter

        failing_backend = collector.PdfBackend(
            "pypdf", "pypdf failing", "python", module=_AlwaysFailingPdfModule
        )
        with tempfile.TemporaryDirectory() as temporary:
            result, output, _checkpoint_path, _pdf_path = self._collect_left_lowe(
                temporary, backends=[failing_backend]
            )
            self.assertEqual(result["returncode"], 1)
            self.assertTrue(output.is_file())
            self.assertIn(
                "pdf_extraction_failed",
                [issue["code"] for issue in result["manifest"]["errors"]],
            )
            artifact = result["manifest"]["settings"][0]["artifacts"]["full_background_subtraction_pdf"]
            self.assertEqual(artifact["extraction_error"], "pdf_extraction_failed")
            self.assertEqual(artifact["extraction_attempts"][0]["backend"], "pypdf")
            with zipfile.ZipFile(output) as archive:
                self.assertIn("manifest.json", archive.namelist())

    def test_backend_fallback_uses_qpdf_after_python_reader_failure(self):
        class _FailingPdfModule:
            PdfReader = type("FailingReader", (), {"__init__": lambda self, _path: (_ for _ in ()).throw(RuntimeError("python backend failed"))})
            PdfWriter = _FakePdfWriter

        python_backend = collector.PdfBackend("pypdf", "pypdf failing", "python", module=_FailingPdfModule)
        qpdf_backend = collector.PdfBackend("qpdf", "qpdf fake", "qpdf", executable="qpdf")

        def qpdf_runner(command, _cwd):
            command = [str(item) for item in command]
            if "--show-npages" in command:
                return {"command": command, "returncode": 0, "stdout": "110\n", "stderr": ""}
            self.assertIn("108-110", command)
            Path(command[-1]).write_bytes(b"qpdf-slim-pdf")
            return {"command": command, "returncode": 0, "stdout": "", "stderr": ""}

        with tempfile.TemporaryDirectory() as temporary:
            source = Path(temporary) / "source.pdf"
            source.write_bytes(b"source")
            result = collector.extract_pdf_pages(
                source,
                backends=[python_backend, qpdf_backend],
                command_runner=qpdf_runner,
            )
        self.assertTrue(result["available"])
        self.assertEqual(result["backend"], "qpdf")
        self.assertEqual(result["pages"], [108, 109, 110])
        self.assertEqual(result["payload"], b"qpdf-slim-pdf")
        self.assertEqual(result["attempts"][0]["backend"], "pypdf")

    def test_discovery_order_prefers_python_then_non_raster_command_backends(self):
        def importer(name):
            if name == "pypdf":
                return _FakePdfModule
            raise ImportError(name)

        executables = {
            "qpdf": "/usr/bin/qpdf",
            "pdfseparate": "/usr/bin/pdfseparate",
            "pdfunite": "/usr/bin/pdfunite",
            "mutool": "/usr/bin/mutool",
        }
        backends = collector.discover_pdf_backends(
            importer=importer, which=lambda name: executables.get(name)
        )
        self.assertEqual(
            [backend.name for backend in backends],
            ["pypdf", "qpdf", "pdfseparate_pdfunite", "mutool"],
        )


if __name__ == "__main__":
    unittest.main()
