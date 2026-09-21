"""Focused contract tests for the standalone E.8 frozen-payload renderer."""

from __future__ import annotations

import ast
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
RENDERER_PATH = REPO_ROOT / "testing" / "render_pion_hgcer_method_a_acceptance_refinement_figure_library.py"
SPEC = importlib.util.spec_from_file_location("_e8_figure_library_renderer", RENDERER_PATH)
renderer = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = renderer
assert SPEC.loader is not None
SPEC.loader.exec_module(renderer)


def _metric(value: float = 1.0) -> dict[str, object]:
    return {"available": True, "value": value, "reason": None}


def _unavailable(reason: str = "low_population_empty") -> dict[str, object]:
    return {"available": False, "value": None, "reason": reason}


def _interval() -> dict[str, object]:
    return {
        "requested_replica_count": 4,
        "valid_replica_count": 4,
        "invalid_replica_count": 0,
        "interval_available": True,
        "ci_low": 0.1,
        "ci_high": 0.3,
        "reason": None,
    }


def _unavailable_interval() -> dict[str, object]:
    return {
        "requested_replica_count": 4,
        "valid_replica_count": 0,
        "invalid_replica_count": 4,
        "interval_available": False,
        "ci_low": None,
        "ci_high": None,
        "reason": "no_valid_bootstrap_replicas",
    }


def _population(values: object, *, available: bool = True, reason: str | None = None) -> dict[str, object]:
    return {"available": available, "unit_area": values, "reason": reason}


def _one_dimensional(empty: bool) -> dict[str, object]:
    if empty:
        population = _population([0.0, 0.0], available=False, reason="low_population_empty")
        metric = _unavailable()
    else:
        population = _population([0.4, 0.6])
        metric = _metric(0.2)
    return {
        "edges": [-1.0, 0.0, 1.0],
        "L": deepcopy(population),
        "B": deepcopy(population),
        "A": deepcopy(population),
        "metrics": {"DeltaH": deepcopy(metric), "kappa": deepcopy(metric)},
    }


def _joint(x_name: str, y_name: str, empty: bool) -> dict[str, object]:
    if empty:
        population = _population([[0.0, 0.0], [0.0, 0.0]], available=False, reason="control_population_empty")
        metric = _unavailable("control_population_empty")
    else:
        population = _population([[0.2, 0.3], [0.1, 0.4]])
        metric = _metric(0.3)
    return {
        "x_variable": x_name,
        "y_variable": y_name,
        "x_edges": [-1.0, 0.0, 1.0],
        "y_edges": [-2.0, 0.0, 2.0],
        "L": deepcopy(population),
        "B": deepcopy(population),
        "A": deepcopy(population),
        "metrics": {"kappa": metric},
    }


def _child(setting_id: str, t_index: int, phi_index: int, phi_low: float, phi_high: float, empty: bool) -> dict[str, object]:
    metric = _unavailable() if empty else _metric(3.0)
    interval = _unavailable_interval() if empty else _interval()
    return {
        "setting_id": setting_id,
        "canonical_t_index": t_index,
        "phi_index": phi_index,
        "phi_low": phi_low,
        "phi_high": phi_high,
        "population_counts": {
            "N_low": 0 if empty else 4,
            "N_control": 0 if empty else 5,
            "N_full_application": 0 if empty else 6,
        },
        "availability": {
            "has_low_response": not empty,
            "has_prompt_control": not empty,
            "has_full_application": not empty,
            "completely_empty": empty,
        },
        "support": {
            "prompt_control": {"ood_fraction": None if empty else 0.2},
            "full_physical_application": {"ood_fraction": None if empty else 0.3},
        },
        "effective_sample_size": {
            "baseline_w0": deepcopy(metric),
            "method_a_w0_times_C": deepcopy(metric),
        },
        "one_dimensional": {
            "analysis_MM": _one_dimensional(empty),
            "SHMS_xptar": _one_dimensional(empty),
            "SHMS_yptar": _one_dimensional(empty),
        },
        "joint_distributions": {
            "analysis_MM__SHMS_xptar": _joint("analysis_MM", "SHMS_xptar", empty),
            "analysis_MM__SHMS_yptar": _joint("analysis_MM", "SHMS_yptar", empty),
            "SHMS_delta__SHMS_xptar": _joint("SHMS_delta", "SHMS_xptar", empty),
            "SHMS_delta__SHMS_yptar": _joint("SHMS_delta", "SHMS_yptar", empty),
        },
        "signed_background": {
            "kaon_window": {
                "P_B_K": deepcopy(metric),
                "P_A_K": deepcopy(metric),
                "DeltaP_K": deepcopy(metric),
                "f_refine_K": deepcopy(metric),
            },
            "variance_proxy": {"R_V": deepcopy(metric)},
        },
        "bootstrap": {
            "one_dimensional": {"analysis_MM": {"DeltaH": deepcopy(interval), "kappa": deepcopy(interval)}},
            "joint_missing_mass_acceptance": {
                "analysis_MM__SHMS_xptar": {"kappa": deepcopy(interval)},
                "analysis_MM__SHMS_yptar": {"kappa": deepcopy(interval)},
            },
            "kaon_window": {"DeltaP_K": deepcopy(interval)},
        },
    }


def _artifact() -> dict[str, object]:
    parents = []
    position = 0
    for setting_index, (phi, epsilon) in enumerate((
        ("Left", "lowe"), ("Left", "highe"), ("Center", "lowe"),
        ("Center", "highe"), ("Right", "highe"),
    )):
        setting_id = "{}-{}".format(phi, epsilon)
        for t_index in range(3):
            phi_edges = [-180.0 + 40.0 * item for item in range(10)]
            children = []
            for phi_index in range(9):
                children.append(_child(setting_id, t_index, phi_index, phi_edges[phi_index], phi_edges[phi_index + 1], position >= 115))
                position += 1
            parents.append({
                "setting": {"phi_setting": phi, "epsilon_setting": epsilon, "ordinal": setting_index},
                "setting_id": setting_id,
                "canonical_t_index": t_index,
                "canonical_t_low": -0.5 + t_index * 0.1,
                "canonical_t_high": -0.4 + t_index * 0.1,
                "phi_edges": phi_edges,
                "children": children,
            })
    validation = {
        "schema_version": renderer.VALIDATION_SCHEMA,
        "status": "available",
        "available": True,
        "fingerprint": renderer.VALIDATION_FINGERPRINT,
        "non_authoritative": True,
        "validation_only": True,
        "event_correction_persisted": False,
        "production_application_performed": False,
        "production_objects_mutated": False,
        "yield_constructed": False,
        "cross_section_constructed": False,
        "root_object_constructed": False,
        "child_renormalization_performed": False,
        "smoothing_or_interpolation_performed": False,
        "absolute_probability_constructed": False,
        "method_b_numerical_dependency": False,
        "automatic_case_classification": False,
        "case_thresholds_defined": False,
        "final_yield_uncertainty_claimed": False,
        "kaon_window": {"mm_min": -0.1, "mm_max": 0.1},
        "parents": parents,
    }
    return {
        "schema_version": renderer.ARTIFACT_SCHEMA,
        "validation": validation,
        "artifact_fingerprint": renderer.ARTIFACT_FINGERPRINT,
        "non_authoritative": True,
        "validation_only": True,
        "manual_review_required": True,
        "production_application_performed": False,
        "production_objects_mutated": False,
        "yield_constructed": False,
        "cross_section_constructed": False,
        "method_b_numerical_dependency": False,
        "automatic_case_classification": False,
    }


def _outputs(root: Path) -> tuple[Path, Path]:
    return root / renderer.PDF_BASENAME, root / renderer.MANIFEST_BASENAME


class FigureLibraryRendererTests(unittest.TestCase):
    maxDiff = None

    def _render(self, root: Path, artifact: dict[str, object] | None = None) -> tuple[dict[str, object], Path, Path]:
        validation = renderer.validate_artifact(_artifact() if artifact is None else artifact)
        pdf, manifest = _outputs(root)
        result = renderer._write_figure_library(
            validation,
            input_name="Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json",
            renderer_commit="a" * 40,
            profile_commit="b" * 40,
            output_pdf=pdf,
            output_manifest=manifest,
        )
        return result, pdf, manifest

    def test_rendered_inventory_manifest_order_and_input_immutability(self):
        artifact = _artifact()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "synthetic-input.json"
            original = json.dumps(artifact, sort_keys=True, allow_nan=False).encode("utf-8")
            input_path.write_bytes(original)
            manifest, pdf, manifest_path = self._render(root, artifact)
            self.assertEqual(input_path.read_bytes(), original)
            self.assertTrue(pdf.is_file()); self.assertTrue(manifest_path.is_file())
            self.assertEqual(manifest["schema_version"], renderer.MANIFEST_SCHEMA)
            self.assertEqual(manifest["pdf_page_count"], 116)
            self.assertEqual((manifest["parent_count"], manifest["canonical_child_count"], manifest["rendered_child_count"], manifest["empty_child_count"]), (15, 135, 115, 20))
            records = manifest["pages"]
            self.assertEqual(len(records), 136)
            self.assertEqual(records[0], {"pdf_page_number": 1, "page_kind": "cover", "rendered": True, "empty": False})
            rendered = [record for record in records if record.get("page_kind") == "child" and record["rendered"]]
            empty = [record for record in records if record.get("page_kind") == "child" and record["empty"]]
            self.assertEqual([record["pdf_page_number"] for record in rendered], list(range(2, 117)))
            self.assertEqual(len(empty), 20); self.assertTrue(all(record["pdf_page_number"] is None for record in empty))
            self.assertEqual((rendered[0]["setting_id"], rendered[0]["canonical_t_index"], rendered[0]["phi_index"]), ("Left-lowe", 0, 0))
            self.assertEqual(json.loads(manifest_path.read_text(encoding="utf-8")), manifest)
            page_objects = re.findall(rb"/Type /Page(?!s)", pdf.read_bytes())
            self.assertEqual(len(page_objects), 116)

    def test_identical_validated_renders_are_byte_identical(self):
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            self._render(Path(first)); self._render(Path(second))
            first_pdf, first_manifest = _outputs(Path(first)); second_pdf, second_manifest = _outputs(Path(second))
            self.assertEqual(hashlib.sha256(first_pdf.read_bytes()).hexdigest(), hashlib.sha256(second_pdf.read_bytes()).hexdigest())
            self.assertEqual(hashlib.sha256(first_manifest.read_bytes()).hexdigest(), hashlib.sha256(second_manifest.read_bytes()).hexdigest())

    def test_literal_unavailable_reason_and_required_persisted_fields_survive(self):
        artifact = _artifact()
        population = artifact["validation"]["parents"][0]["children"][0]["one_dimensional"]["analysis_MM"]["L"]
        population.update({"available": False, "unit_area": [0.0, 0.0], "reason": "persisted_literal_reason"})
        validation = renderer.validate_artifact(artifact)
        self.assertEqual(renderer._unavailable_text(validation["parents"][0]["children"][0]["one_dimensional"]["analysis_MM"]["L"]), "unavailable:persisted_literal_reason")
        source = RENDERER_PATH.read_text(encoding="utf-8")
        for required in (
            "SHMS_delta__SHMS_xptar", "SHMS_delta__SHMS_yptar",
            "analysis_MM__SHMS_xptar", "analysis_MM__SHMS_yptar",
            "Normalized bin fraction", "descriptive acceptance refinement",
        ):
            self.assertIn(required, source)

    def test_validation_rejects_malformed_or_non_authoritative_payloads(self):
        cases = {
            "schema": lambda value: value.__setitem__("schema_version", "wrong"),
            "available": lambda value: value["validation"].__setitem__("available", False),
            "artifact_fingerprint": lambda value: value.__setitem__("artifact_fingerprint", "0" * 64),
            "validation_fingerprint": lambda value: value["validation"].__setitem__("fingerprint", "0" * 64),
            "production": lambda value: value["validation"].__setitem__("production_application_performed", True),
            "method_b": lambda value: value.__setitem__("method_b_numerical_dependency", True),
            "parents": lambda value: value["validation"].__setitem__("parents", value["validation"]["parents"][:-1]),
            "inventory": lambda value: value["validation"]["parents"][0]["children"][0]["availability"].__setitem__("completely_empty", True),
            "one_d_edges": lambda value: value["validation"]["parents"][0]["children"][0]["one_dimensional"]["analysis_MM"].__setitem__("edges", [0.0, 0.0]),
            "one_d_shape": lambda value: value["validation"]["parents"][0]["children"][0]["one_dimensional"]["analysis_MM"]["L"].__setitem__("unit_area", [1.0]),
            "nonfinite": lambda value: value["validation"]["parents"][0]["children"][0]["one_dimensional"]["analysis_MM"]["L"].__setitem__("unit_area", [float("nan"), 0.0]),
            "two_d_edges": lambda value: value["validation"]["parents"][0]["children"][0]["joint_distributions"]["SHMS_delta__SHMS_xptar"].__setitem__("x_edges", [0.0, 0.0]),
            "two_d_shape": lambda value: value["validation"]["parents"][0]["children"][0]["joint_distributions"]["SHMS_delta__SHMS_xptar"]["L"].__setitem__("unit_area", [[1.0]]),
        }
        for label, mutate in cases.items():
            with self.subTest(label=label):
                value = _artifact(); mutate(value)
                with self.assertRaises(renderer.FigureLibraryError):
                    renderer.validate_artifact(value)

    def test_cli_rejects_sha_paths_commits_and_existing_outputs_without_writing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "input.json"; original = b"not-the-accepted-json"; input_path.write_bytes(original)
            pdf, manifest = _outputs(root)
            arguments = [
                "--input-json", str(input_path), "--output-pdf", str(pdf), "--output-manifest", str(manifest),
                "--renderer-source-commit", "a" * 40, "--bundle-profile-commit", "b" * 40,
            ]
            self.assertEqual(renderer.main(arguments), 1)
            self.assertEqual(input_path.read_bytes(), original); self.assertFalse(pdf.exists()); self.assertFalse(manifest.exists())
            self.assertEqual(renderer.main(arguments[:-4] + ["--renderer-source-commit", "A" * 40, "--bundle-profile-commit", "b" * 40]), 1)
            wrong_pdf = root / "wrong.pdf"
            self.assertEqual(renderer.main(["--input-json", str(input_path), "--output-pdf", str(wrong_pdf), "--output-manifest", str(manifest), "--renderer-source-commit", "a" * 40, "--bundle-profile-commit", "b" * 40]), 1)
            pdf.write_bytes(b"preexisting")
            self.assertEqual(renderer.main(arguments), 1)
            self.assertEqual(pdf.read_bytes(), b"preexisting")
            self.assertEqual(renderer.main(["--input-json", str(input_path), "--output-pdf", str(input_path), "--output-manifest", str(manifest), "--renderer-source-commit", "a" * 40, "--bundle-profile-commit", "b" * 40]), 1)

    def test_parser_hashes_before_parsing_and_rejects_malformed_json(self):
        malformed = b"not-json"
        with self.assertRaisesRegex(renderer.FigureLibraryError, "input_sha256_mismatch"):
            renderer._parse_input_bytes(malformed)
        with self.assertRaisesRegex(renderer.FigureLibraryError, "input_json_invalid"):
            renderer._parse_input_bytes(malformed, expected_sha256=hashlib.sha256(malformed).hexdigest())

    def test_atomic_cleanup_removes_partially_promoted_outputs(self):
        validation = renderer.validate_artifact(_artifact())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); pdf, manifest = _outputs(root)
            original_promote = renderer._promote_without_overwrite

            def short_pdf(path, *_args):
                path.write_bytes(b"pdf")

            def fail_second(temporary_path, destination):
                if destination == manifest:
                    raise OSError("manifest_promotion_failure")
                return original_promote(temporary_path, destination)

            with mock.patch.object(renderer, "_render_pdf", side_effect=short_pdf), mock.patch.object(renderer, "_promote_without_overwrite", side_effect=fail_second):
                with self.assertRaisesRegex(OSError, "manifest_promotion_failure"):
                    renderer._write_figure_library(validation, input_name="input.json", renderer_commit="a" * 40, profile_commit="b" * 40, output_pdf=pdf, output_manifest=manifest)
            self.assertFalse(pdf.exists()); self.assertFalse(manifest.exists())
            self.assertEqual(list(root.glob(".e8_*")), [])

    def test_static_ownership_and_changed_file_boundary(self):
        source = RENDERER_PATH.read_text(encoding="utf-8")
        tree = ast.parse(source)
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module is not None:
                imported.add(node.module)
        forbidden = {
            "ROOT", "PyROOT", "rand_sub", "full_background_subtraction_plots",
            "pion_hgcer_method_a_acceptance_refinement_validation",
            "analyze_pion_hgcer_method_a_acceptance_refinement_validation",
        }
        self.assertFalse(any(name.split(".")[-1] in forbidden for name in imported))
        self.assertFalse(any(name.startswith("pion_hgcer_") for name in imported))
        self.assertNotIn("calculate_yield", source)
        self.assertNotIn("calculate_xsection", source)
        status = subprocess.run(["git", "status", "--short"], cwd=REPO_ROOT, text=True, capture_output=True, check=True).stdout.splitlines()
        changed = {line[3:].replace("\\", "/") for line in status if len(line) >= 4}
        self.assertEqual(changed, {
            "testing/render_pion_hgcer_method_a_acceptance_refinement_figure_library.py",
            "testing/test_render_pion_hgcer_method_a_acceptance_refinement_figure_library.py",
        })


if __name__ == "__main__":
    unittest.main()
