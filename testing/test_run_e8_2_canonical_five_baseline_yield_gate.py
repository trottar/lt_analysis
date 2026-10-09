"""Synthetic owner/filesystem fixtures; no launcher, ROOT or farm execution."""
import ast
from contextlib import contextmanager, ExitStack
from copy import deepcopy
import io
import hashlib
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock, patch
import zipfile
from types import SimpleNamespace

import numpy as np

from testing import run_e8_2_canonical_five_baseline_yield_gate as owner

ROOT = Path(__file__).resolve().parents[1]
SHA = "a" * 40
T = [0.0, 0.1, 0.2, 0.3]
PHI = list(range(10))
MM = np.linspace(1.1, 1.14, 101).tolist()


def write_json(path, value):
    path.write_text(json.dumps(value), encoding="utf-8")


def fixture(worktree, outdir, skip=False, distinct_children=False):
    (worktree / "src/kaon/yields").mkdir(parents=True, exist_ok=True)
    (worktree / "log").mkdir(exist_ok=True)
    for eps, phis in (("lowe", ["Center", "Left"]), ("highe", ["Center", "Left", "Right"])):
        inp = {"ParticleType": "kaon", "EPSSET": eps[:-1], "Q2": "4p4", "W": "2p74", "POL": 1,
               "EPSVAL": 0.3 if eps == "lowe" else 0.7, "OutFilename": f"FullAnalysis_{owner.KINEMATIC}_{eps}",
               "iter_num": 0, "mm_min": MM[0], "mm_max": MM[-1], "bg_active_profile": "no_empirical_residual",
               "bg_stat_scale1": 0, "bg_stat_scale2": 0, "model_provenance": {"source": "fixture-only"}, "_private_cache": {"events": [1]}}
        histlist = []
        support = {"t_bins": T, "phi_bins": PHI, "mm_edges": MM, "settings": phis}
        for phi in phis:
            inp["runNum" + phi] = "1"
            inp["pThetaVal" + phi] = [20 if phi == "Center" else 21 if phi == "Left" else 19]
            (worktree / "log" / f"{phi}_kaon_1_{owner.KINEMATIC}_{eps}.log").write_text("fixture")
            parents = [{"t_bin": j + 1, "t_edges": T[j:j + 2], "pion_parent_id": f"{phi}-{eps}-{j}",
                        "diagnostic_application_status": {"final_status": "skip_bin" if skip and j == 1 else "applied_component",
                        "fallback_mode": "skip_bin" if skip and j == 1 else "error",
                        "production_evaluation": "rejected" if skip and j == 1 else "accepted"}} for j in range(3)]
            hist = {"phi_setting": phi, "t_bins": T, "phi_bins": PHI,
                    "pion_t_parent_collection_frozen": True, "pion_t_amplitude_table": parents}
            histlist.append(hist)
            pdf = owner.isolation.pdf_name(phi, eps)
            (outdir / pdf).write_bytes(b"%PDF-1.4\nfixture")
            pages = [{"page_id": pid, "scope": f"t{n}", "t_index": n - 1, "semantic_stage": semantic,
                      "setting": phi, "epsilon": eps[:-1], "t_edges": T[n - 1:n + 1],
                      "represented_phi_inventory": [{"phi_index": k, "phi_edges": PHI[k:k + 2]} for k in range(9)],
                      "invalid_unavailable_children": [{"phi_index": k, "reason": "fixture frozen skip"} for k in range(9)]
                         if skip and n == 2 else [], "authoritative": False, "presentation_only": True}
                     for pid, (n, semantic) in owner.baseline.PAGE_IDS.items()]
            pages.extend([{"page_id": "full_background.e8_4.unavailable", "reason": "fixture_candidate_missing"},
                          {"page_id": "full_background.e8.handoff"}])
            write_json(outdir / pdf.replace(".pdf", "-manifest.json"), {
                "schema_version": "full_background_subtraction_page_manifest/v1", "pdf_basename": pdf,
                "setting": {"kinematic_token": owner.KINEMATIC, "epsilon_setting": eps[:-1],
                "epsilon_filename_token": eps, "phi_setting": phi, "particle_type": "kaon"},
                "renderer_failures": [], "pages": pages})
            for kind in ("data", "simc"):
                support[f"{kind}_mm_{phi.lower()}_values"] = np.broadcast_to([1.123456, 2.0] + [0.0] * 98, (3, 9, 100)).copy()
                support[f"{kind}_mm_{phi.lower()}_errors"] = np.broadcast_to([0.1, 0.2] + [0.0] * 98, (3, 9, 100)).copy()
                if distinct_children:
                    for j in range(3):
                        for k in range(9):
                            # Distinct data/SIMC and t/phi signatures expose wrong-child lookup.
                            value = 3.123456 + 100 * j + k if kind == "data" else 5.23456 + 200 * j + 2 * k
                            support[f"{kind}_mm_{phi.lower()}_values"][j, k] = [value] + [0.0] * 99
        # Filename producer needs Center angles even for side settings.
        for phi in phis:
            for kind in ("data", "simc"):
                name = owner.table_name(worktree, inp, phi, kind)
                rows = []
                for j in range(3):
                    for k in range(9):
                        if kind == "data" and skip and j == 1:
                            continue
                        value = (3.123456 + 100 * j + k if kind == "data" else 5.23456 + 200 * j + 2 * k) if distinct_children else 3.123456
                        rows.append(f"{value:.4f} {('7.5000' if kind == 'data' else '0.2236')} {k + 1} {j + 1}\n")
                (worktree / "src/kaon/yields" / name).write_text("".join(rows))
        write_json(outdir / f"kaon_FullAnalysis_{owner.KINEMATIC}_{eps}.json", {"inpDict": inp, "histlist": histlist})
        base = outdir / f"kaon_FullAnalysis_{owner.KINEMATIC}_{eps}_correction_ledger_no_empirical_residual"
        write_json(base.with_suffix(".json"), {"active_profile": "no_empirical_residual", "particle_type": "kaon",
            "epsset": eps[:-1], "q2": "4p4", "w": "2p74", "outfilename": inp["OutFilename"],
            "mm_cut_window": [MM[0], MM[-1]], "settings": [{"phi_setting": p} for p in phis]})
        base.with_suffix(".csv").write_text("row_kind,phi_setting\n" + "".join(f"setting_total,{p}\n" for p in phis))
        np.savez_compressed(outdir / f"kaon_xsect_support_{owner.KINEMATIC}_{eps}.npz", **support)


def fake_zip(path, records, profile, commit=SHA):
    manifest = {"complete": True, "errors": [], "git_head": commit, "required_analysis_commit": commit,
                "validation_profile": owner.PROFILE_ID, "requested_settings": owner.SETTINGS,
                "requested_kinematic": owner.KINEMATIC, "global_artifacts": {}, "settings": []}
    rows = {}
    for row in owner.SETTINGS:
        prefix = row["phi"] + "_" + row["epsilon"]
        item = {**row, "kinematic": owner.KINEMATIC, "artifacts": {}}
        rows[prefix] = item
        manifest["settings"].append(item)
    with zipfile.ZipFile(path, "w") as archive:
        for prefix in ["global", *rows]:
            archive.writestr(prefix + "/", "")
        for name, entry in owner.names(profile).items():
            record = records[Path(name).name]
            archived = {"archive_path": name, "sha256": record["sha256"], "byte_size": record["bytes"]}
            prefix = name.split("/")[0]
            target = manifest["global_artifacts"] if prefix == "global" else rows[prefix]["artifacts"]
            target[entry["key"]] = archived
            archive.writestr(name, (path.parent / "volatile/OUTPUT/Analysis/KaonLT" / Path(name).name).read_bytes())
        archive.writestr("manifest.json", json.dumps(manifest))
        archive.writestr("source_state.txt", "fixture source")
        archive.writestr("source_checks.txt", "fixture checks")


class OwnerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.worktree = self.root / "analysis"
        self.worktree.mkdir()
        self.outdir = self.root / "volatile/OUTPUT/Analysis/KaonLT"
        self.outdir.mkdir(parents=True)
        self.output = self.root / "evidence.zip"
        self.profile = owner.resolved_profile(ROOT, SHA, self.output)

    def test_authentic_tables_support_and_missing_method_a(self):
        fixture(self.worktree, self.outdir)
        records, cells, meta = owner.audit_outputs(self.worktree, self.outdir, self.profile, 0, {})
        self.assertEqual(len(cells), 135)
        self.assertEqual(len(records), 28)
        self.assertEqual(cells[0]["data"]["error"], 7.5)
        self.assertEqual(cells[0]["data_values"][:2], [1.123456, 2.0])
        self.assertNotIn("_private_cache", meta["lowe"]["persisted_input"])
        self.assertEqual(meta["lowe"]["settings"]["Left"]["e8_4_unavailable_reasons"], ["fixture_candidate_missing"])
        for record in records.values():
            if "source_relative_basename" in record:
                self.assertEqual((self.worktree / record["source_relative_basename"]).read_bytes(),
                                 (self.outdir / record["copied_basename"]).read_bytes())

    def test_explicit_skip_is_excluded_and_zero_is_valid(self):
        fixture(self.worktree, self.outdir, skip=True)
        _, cells, _ = owner.audit_outputs(self.worktree, self.outdir, self.profile, 0, {})
        skipped = [c for c in cells if c["parent_state"] == "skip_bin"]
        self.assertEqual(len(skipped), 45)
        self.assertTrue(all(c["data"] is None and c["simc_numerical_gate"] == "PASS" for c in skipped))
        self.assertTrue(all(c["t_index"] == 1 and c["simc"]["producer_t_bin"] == 2 and
                            c["simc"]["producer_phi_bin"] == c["phi_index"] + 1 for c in skipped))
        hist = {"pion_t_parent_collection_frozen": True, "pion_t_amplitude_table": [
            {"t_bin": j + 1, "t_edges": T[j:j + 2], "diagnostic_application_status":
             {"final_status": "zero", "fallback_mode": "zero", "production_evaluation": "rejected"}} for j in range(3)]}
        self.assertEqual(owner.parent_states(hist, T), ["zero"] * 3)
        hist["pion_t_amplitude_table"][0]["diagnostic_application_status"]["final_status"] = "unavailable"
        self.assertEqual(owner.parent_states(hist, T), ["zero"] * 3)
        hist["pion_t_amplitude_table"][0]["diagnostic_application_status"]["fallback_mode"] = "error"
        with self.assertRaisesRegex(Exception, "policy_not_verified"):
            owner.parent_states(hist, T)

    def test_quantization_signed_data_and_nonnegative_simc(self):
        data = {"yield_text": "-1.2346", "error": 8.0}
        simc = {"yield_text": "0.0000", "error_text": "0.0000"}
        self.assertEqual(owner.check_cell([-1.23456], [-2], [3], data, simc)["simc_nonnegative_yield"], 0)
        self.assertEqual(data["error"], 8.0)
        with self.assertRaisesRegex(Exception, "data_yield_support_mismatch"):
            owner.check_cell([-1.2344], [-2], [3], data, simc)
        with self.assertRaisesRegex(Exception, "simc_yield_support_mismatch"):
            owner.check_simc([1], [0], simc)

    def test_authentic_one_based_table_converts_and_preserves_raw_coordinates(self):
        path = self.root / "table.dat"
        path.write_text("".join(f"1.0000 1.0000 {k + 1} {j + 1}\n" for j in range(3) for k in range(9)))
        expected = {(j, k) for j in range(3) for k in range(9)}
        rows = owner.read_table(path, expected)
        self.assertEqual(set(rows), expected)
        for (j, k), row in rows.items():
            self.assertEqual((row["producer_phi_bin"], row["producer_t_bin"]), (k + 1, j + 1))

    def test_zero_based_and_out_of_range_tables_fail(self):
        path = self.root / "table.dat"
        zero_based = "".join(f"1.0000 1.0000 {k} {j}\n" for j in range(3) for k in range(9))
        for text in (zero_based, "1.0000 1.0000 0 1\n", "1.0000 1.0000 1 0\n",
                     "1.0000 1.0000 10 1\n", "1.0000 1.0000 1 4\n"):
            with self.subTest(text=text):
                path.write_text(text)
                with self.assertRaisesRegex(Exception, "yield_producer_coordinate_out_of_range"):
                    owner.read_table(path, {(j, k) for j in range(3) for k in range(9)})

    def test_table_swaps_duplicates_omissions_order_and_nonfinite(self):
        path = self.root / "table.dat"
        expected = {(j, k) for j in range(3) for k in range(9)}
        authentic = [f"1.0000 1.0000 {k + 1} {j + 1}\n" for j in range(3) for k in range(9)]
        swapped = [f"1.0000 1.0000 {j + 1} {k + 1}\n" for j in range(3) for k in range(9)]
        cases = [
            ("phi/t swap", swapped, "yield_producer_coordinate_out_of_range"),
            ("duplicate", [authentic[0], *authentic], "duplicate_yield_coordinate"),
            ("missing cell", authentic[:-1], "yield_coordinate_inventory_or_order_invalid"),
            ("order", list(reversed(authentic)), "yield_coordinate_inventory_or_order_invalid"),
            ("nonfinite", ["nan 1.0000 1 1\n"], "yield_table_format_invalid"),
            ("negative error", ["1.0000 -1000.0000 1 1\n"], "yield_table_undefined_error"),
        ]
        for name, lines, reason in cases:
            with self.subTest(name=name):
                path.write_text("".join(lines))
                with self.assertRaisesRegex(Exception, reason):
                    owner.read_table(path, expected)
        # An in-range axis swap must also fail the expected internal inventory.
        path.write_text("1.0000 1.0000 2 1\n")
        with self.assertRaisesRegex(Exception, "yield_coordinate_inventory_or_order_invalid"):
            owner.read_table(path, {(1, 0)})

    def test_converted_coordinates_select_distinct_data_and_simc_npz_children(self):
        fixture(self.worktree, self.outdir, distinct_children=True)
        _, cells, _ = owner.audit_outputs(self.worktree, self.outdir, self.profile, 0, {})
        for cell in cells:
            j, k = cell["t_index"], cell["phi_index"]
            yd, ys = 3.123456 + 100 * j + k, 5.23456 + 200 * j + 2 * k
            self.assertEqual(cell["data_values"], [yd] + [0.0] * 99)
            self.assertEqual(cell["simc_values"], [ys] + [0.0] * 99)
            self.assertEqual(cell["data"]["yield_text"], f"{yd:.4f}")
            self.assertEqual(cell["simc"]["yield_text"], f"{ys:.4f}")
            for kind in ("data", "simc"):
                self.assertEqual((cell[kind]["producer_phi_bin"], cell[kind]["producer_t_bin"]), (k + 1, j + 1))
            self.assertEqual(cell["numerical_gate"], "PASS")
            self.assertEqual(cell["simc_numerical_gate"], "PASS")

    def test_zero_parent_policy_keeps_one_based_populated_yield_inventory(self):
        fixture(self.worktree, self.outdir, distinct_children=True)
        def zero_policy(payload):
            for hist in payload["histlist"]:
                hist["pion_t_amplitude_table"][1]["diagnostic_application_status"] = {
                    "final_status": "zero", "fallback_mode": "zero", "production_evaluation": "rejected"}
        self.change_json("kaon_FullAnalysis_Q4p4W2p74_lowe.json", zero_policy)
        _, cells, _ = owner.audit_outputs(self.worktree, self.outdir, self.profile, 0, {})
        zeros = [c for c in cells if c["parent_state"] == "zero"]
        self.assertEqual(len(zeros), 18)
        self.assertTrue(all(c["t_index"] == 1 and c["data"]["producer_t_bin"] == 2 and
                            c["simc"]["producer_t_bin"] == 2 and c["numerical_gate"] == "PASS" for c in zeros))

    def test_mutated_runtime_artifacts_fail_closed(self):
        mutations = [
            ("pdf", lambda: (self.outdir / owner.isolation.pdf_name("Left", "lowe")).write_bytes(b"bad")),
            ("json", lambda: self.change_json("kaon_FullAnalysis_Q4p4W2p74_lowe.json", lambda p: p["inpDict"].update(EPSSET="high"))),
            ("iteration", lambda: self.change_json("kaon_FullAnalysis_Q4p4W2p74_lowe.json", lambda p: p["inpDict"].update(iter_num=1))),
            ("manifest", lambda: self.change_json("Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json",
                                                   lambda p: p["setting"].update(phi_setting="Right"))),
            ("csv", lambda: (self.outdir / "kaon_FullAnalysis_Q4p4W2p74_lowe_correction_ledger_no_empirical_residual.csv").write_text("row_kind,phi_setting\nsetting_total,Right\n")),
            ("support_axes", lambda: self.change_npz("t_bins", [0, 0.15, 0.2, 0.3])),
            ("support_settings", lambda: self.change_npz("settings", ["Left", "Right"])),
            ("support_shape", lambda: self.change_npz("data_mm_left_values", np.zeros((3, 8, 2)))),
            ("support_nonfinite", lambda: self.change_npz("simc_mm_left_values", np.full((3, 9, 2), np.nan))),
            ("missing_table", lambda: (self.worktree / "src/kaon/yields/yield_data.pl_Q44W274_30_+1000.dat").unlink())]
        for name, mutate in mutations:
            with self.subTest(name=name):
                for path in self.outdir.iterdir():
                    path.unlink()
                fixture(self.worktree, self.outdir)
                mutate()
                with self.assertRaises(Exception):
                    owner.audit_outputs(self.worktree, self.outdir, self.profile, 0, {})

    def change_json(self, name, mutate):
        path = self.outdir / name
        payload = json.loads(path.read_text())
        mutate(payload)
        write_json(path, payload)

    def change_npz(self, key, value):
        path = self.outdir / "kaon_xsect_support_Q4p4W2p74_lowe.npz"
        with np.load(path) as inp:
            payload = {k: inp[k] for k in inp.files}
        payload[key] = value
        np.savez_compressed(path, **payload)

    def test_stale_outputs_duplicate_json_and_method_a_reason(self):
        path = self.root / "file.json"
        path.write_text('{"key":1,"key":2}')
        with self.assertRaisesRegex(Exception, "duplicate_json_key"):
            owner.strict_json(path)
        with self.assertRaisesRegex(Exception, "artifact_stale"):
            owner.fresh(path, path.stat().st_mtime_ns + 1, {})
        fixture(self.worktree, self.outdir)
        self.change_json("Left_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json",
                         lambda p: p["pages"][-2].update(reason=None))
        with self.assertRaisesRegex(Exception, "reason_not_persisted"):
            owner.audit_outputs(self.worktree, self.outdir, self.profile, 0, {})

    def test_exact_launcher_and_completion_markers(self):
        log = self.root / "child.log"
        process = Mock(stdout=io.BytesIO(b"Low Epsilon Completed!\nHigh Epsilon Completed!\n"))
        process.wait.return_value = 0
        with patch.object(owner.isolation.subprocess, "Popen", return_value=process) as launch:
            owner.run_analysis(self.worktree, log, {"fixture": "env"})
        self.assertEqual(launch.call_args.args[0], owner.COMMAND)
        self.assertEqual(launch.call_args.kwargs["cwd"], self.worktree)
        owner.verify_completion(log, 0)
        for text, rc in (("Low Epsilon Completed!", 0), ("High Epsilon Completed!", 0),
                         (log.read_text(), 1), (log.read_text() + "Left / lowe debug analysis completed.", 0)):
            log.write_text(text)
            with self.assertRaises(Exception):
                owner.verify_completion(log, rc)

    def test_preflight_wrong_identity_dirty_and_preservation(self):
        clean = {"branch": "test", "head": SHA, "origin_test": SHA, "porcelain": ""}
        for key, value in (("branch", "main"), ("head", "b" * 40), ("origin_test", "b" * 40),
                           ("porcelain", " M src/main.py")):
            with patch.object(owner, "snapshot", return_value={**clean, key: value}), patch.object(owner, "git"):
                with self.assertRaises(Exception):
                    owner.preflight(ROOT, SHA)
        with patch.object(owner, "snapshot", return_value=clean), patch.object(owner, "git"):
            self.assertEqual(owner.preflight(ROOT, SHA), clean)
        with patch.object(owner, "snapshot", return_value={**clean, "content": "mutated"}):
            with self.assertRaisesRegex(Exception, "bytes_changed"):
                owner.preservation(ROOT, clean)

    def test_operational_success_and_failure_cleanup(self):
        @contextmanager
        def worktree(repo, commit, role):
            entered.append(role)
            try:
                yield self.worktree
            finally:
                cleaned.append(role)
        paths = {"VOLATILEPATH": str(self.root / "volatile")}
        for fail_stage in (None, "source_preflight", "prepare_runtime_overlay", "probe_paths",
                           "external_symlink_preflight", "verify_ltsep_preservation",
                           "analysis", "audit", "collector", "zip"):
            with self.subTest(stage=fail_stage):
                for path in self.outdir.iterdir():
                    path.unlink()
                if self.output.exists():
                    self.output.unlink()
                entered, cleaned = [], []
                overlay = {"baseline_paths": paths}
                def run(*args):
                    args[1].write_text("Low Epsilon Completed!\nHigh Epsilon Completed!\n")
                    if fail_stage == "analysis":
                        raise RuntimeError("synthetic analysis failure")
                    return 0
                def collect(**args):
                    if fail_stage == "collector":
                        return {"returncode": 1}
                    self.output.write_bytes(b"synthetic zip")
                    return {"returncode": 0}
                module = Mock(collect_validation_bundle=collect)
                @contextmanager
                def collection(source):
                    yield module
                with ExitStack() as stack:
                    for name, value in (("preflight", {}), ("source_preflight", []),
                                        ("prepare_runtime_overlay", ({}, overlay)),
                                        ("probe_paths", [{"paths": paths}]), ("external_symlink_preflight", {"passed": True}),
                                        ("preservation", {"passed": True}), ("verify_ltsep_preservation", {"passed": True})):
                        operational = stack.enter_context(patch.object(owner, name, return_value=value))
                        if fail_stage == name:
                            operational.side_effect = RuntimeError("synthetic " + name + " failure")
                    stack.enter_context(patch.object(owner, "owned_worktree", side_effect=worktree))
                    stack.enter_context(patch.object(owner, "collection_module", side_effect=collection))
                    stack.enter_context(patch.object(owner, "run_analysis", side_effect=run))
                    audit = stack.enter_context(patch.object(owner, "audit_outputs", return_value=({}, [], {})))
                    verify = stack.enter_context(patch.object(owner, "verify_zip"))
                    delivery = stack.enter_context(patch.object(owner.isolation, "deliver_evidence"))
                    if fail_stage == "audit":
                        audit.side_effect = RuntimeError("synthetic audit failure")
                    if fail_stage == "zip":
                        verify.side_effect = RuntimeError("synthetic zip failure")
                    if fail_stage:
                        with self.assertRaises(Exception):
                            owner.execute_gate(ROOT, self.outdir, self.output, SHA)
                        delivery.assert_not_called()
                    else:
                        self.assertEqual(owner.execute_gate(ROOT, self.outdir, self.output, SHA), self.output)
                        delivery.assert_called_once()
                    self.assertEqual(entered, cleaned)

    def test_zip_inventory_and_corruption(self):
        fixture(self.worktree, self.outdir)
        records, _, _ = owner.audit_outputs(self.worktree, self.outdir, self.profile, 0, {})
        write_json(self.outdir / owner.artifact_name(self.profile, "run_summary"), {"schema_version": "fixture"})
        records[owner.artifact_name(self.profile, "run_summary")] = {"sha256": owner.sha256(self.outdir / owner.artifact_name(self.profile, "run_summary")),
                                  "bytes": (self.outdir / owner.artifact_name(self.profile, "run_summary")).stat().st_size}
        fake_zip(self.output, records, self.profile)
        owner.verify_zip(self.output, SHA, records, self.profile)
        with zipfile.ZipFile(self.output, "a") as archive:
            archive.writestr("unexpected.dat", "unexpected")
        with self.assertRaisesRegex(Exception, "archive_inventory_invalid"):
            owner.verify_zip(self.output, SHA, records, self.profile)

    def test_no_candidate_gate_calls_or_scientific_imports(self):
        tree = ast.parse((ROOT / "testing/run_e8_2_canonical_five_baseline_yield_gate.py").read_text())
        calls = {node.func.attr for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)}
        self.assertFalse(calls & {"candidate_preflight", "stage_candidates", "restore_candidates", "verify_scientific_equivalence"})

    def attempt_harness(self, stack, *, partial_failure=False):
        """Run real capture/receipt/collector/ZIP checks with a fake farm only."""
        from testing.test_pion_hgcer_validation_bundle_profile_e8_2_canonical_five_yields import source_runner
        primary, external = self.root / "primary-config", self.root / "installed-ltsep-config"
        if not primary.exists():
            primary.write_bytes(b"ordinary checkout identity")
        if not external.exists():
            external.write_bytes(b"external installed configuration")
        frozen = {p: p.read_bytes() for p in (primary, external)}
        producers, cleanups = [], []

        def clock():
            # Explicit synthetic launch clock, independent of coarse/host
            # filesystem timestamps. Successive attempts get distinct mtimes.
            self.fixture_clock = getattr(self, "fixture_clock", 1_600_000_000_000_000_000) + 1_000_000_000
            return self.fixture_clock

        stack.enter_context(patch.object(owner.time, "time_ns", side_effect=clock))

        @contextmanager
        def worktree(repo, commit, role):
            with tempfile.TemporaryDirectory(prefix="owned-" + role + "-", dir=self.root) as temporary:
                try:
                    yield Path(temporary)
                finally:
                    cleanups.append(role)

        def run(worktree, log, env):
            fixture(worktree, self.outdir)
            if partial_failure:
                source = worktree / "src/kaon/yields/yield_data.pl_Q44W274_30_+0000.dat"
                source.write_text(source.read_text().replace("3.1235", "4.1235", 1))
            raw_names = {Path(name).name for name, entry in owner.names(self.profile).items()
                         if entry["key"] != "run_summary" and not entry["key"].startswith("yield_")}
            timestamp = self.fixture_clock + 1_000_000
            for path in [*(self.outdir / name for name in raw_names),
                         *(worktree / "src/kaon/yields").iterdir()]:
                os.utime(path, ns=(timestamp, timestamp))
            producers.append({p.name: p.read_bytes() for p in (worktree / "src/kaon/yields").iterdir()})
            with log.open("xb") as handle:
                handle.write(b"Low Epsilon Completed!\nHigh Epsilon Completed!\n")
            return 0

        def preserve(*args):
            for path, content in frozen.items():
                self.assertEqual(path.read_bytes(), content)
            return {"passed": True}

        @contextmanager
        def collection(source):
            yield SimpleNamespace(collect_validation_bundle=lambda **args:
                owner.collector.collect_validation_bundle(**args, command_runner=source_runner))

        paths = {"VOLATILEPATH": str(self.root / "volatile")}
        for name, value in (("preflight", {}), ("source_preflight", []),
                            ("prepare_runtime_overlay", ({}, {"baseline_paths": paths})),
                            ("probe_paths", [{"paths": paths}]), ("external_symlink_preflight", {"passed": True})):
            stack.enter_context(patch.object(owner, name, return_value=value))
        stack.enter_context(patch.object(owner, "preservation", side_effect=preserve))
        stack.enter_context(patch.object(owner, "verify_ltsep_preservation", side_effect=preserve))
        stack.enter_context(patch.object(owner, "owned_worktree", side_effect=worktree))
        stack.enter_context(patch.object(owner, "collection_module", side_effect=collection))
        stack.enter_context(patch.object(owner, "run_analysis", side_effect=run))
        return producers, cleanups, frozen

    def attempt_files(self, output):
        token = owner.attempt_identity(output)["attempt_id"]
        return {p: p.read_bytes() for p in self.outdir.iterdir() if token in p.name}

    def test_two_successive_attempts_preserve_evidence_and_exact_zip_profiles(self):
        with ExitStack() as stack:
            producers, cleanups, frozen = self.attempt_harness(stack)
            first = self.root / "first.zip"
            second = self.root / "second.zip"
            # Legacy fixed-name evidence is also immutable and cannot block a new attempt.
            legacy = self.outdir / "Q4p4W2p74_e8_2_baseline-five-yield-receipt.json"
            legacy.write_bytes(b"prior failed diagnostic receipt")
            owner.execute_gate(ROOT, self.outdir, first, SHA)
            old = self.attempt_files(first)
            first_zip = first.read_bytes()
            owner.execute_gate(ROOT, self.outdir, second, SHA)
            self.assertTrue(old)
            for path, content in old.items():
                self.assertEqual(path.read_bytes(), content)
            self.assertEqual(first.read_bytes(), first_zip)
            self.assertEqual(legacy.read_bytes(), b"prior failed diagnostic receipt")
            self.assertTrue(set(old).isdisjoint(self.attempt_files(second)))
            self.assertEqual(cleanups, ["analysis", "source", "analysis", "source"])
            for index, output in enumerate((first, second)):
                profile = owner.resolved_profile(ROOT, SHA, output)
                receipt = json.loads((self.outdir / owner.artifact_name(profile, "run_summary")).read_text())
                self.assertEqual(receipt["attempt_id"], owner.attempt_identity(output)["attempt_id"])
                self.assertEqual(receipt["requested_output"], str(output))
                self.assertEqual(len(receipt["cells"]), 135)
                self.assertFalse(receipt["boundaries"]["runtime_acceptance"])
                records = dict(receipt["artifacts"])
                summary = self.outdir / owner.artifact_name(profile, "run_summary")
                records[summary.name] = {"sha256": owner.sha256(summary), "bytes": summary.stat().st_size}
                owner.verify_zip(output, SHA, records, profile)
                copied = [r for r in records.values() if "source_relative_basename" in r]
                self.assertEqual(len(copied), 10)
                for record in copied:
                    source_bytes = producers[index][Path(record["source_relative_basename"]).name]
                    self.assertEqual(hashlib.sha256(source_bytes).hexdigest(), record["sha256"])
                    self.assertEqual((self.outdir / record["copied_basename"]).read_bytes(), source_bytes)
                    self.assertGreater(record["source_mtime_ns"], 0)
                    self.assertGreaterEqual(record["copied_mtime_ns"], record["source_mtime_ns"])
            with self.assertRaises(Exception):
                owner.verify_zip(first, SHA, records, profile)
            for path, content in frozen.items():
                self.assertEqual(path.read_bytes(), content)

    def test_failed_partial_attempt_then_fresh_attempt_preserves_every_owned_file(self):
        failed = self.root / "failed.zip"
        with ExitStack() as stack:
            _, cleanups, _ = self.attempt_harness(stack, partial_failure=True)
            with self.assertRaisesRegex(Exception, "data_yield_support_mismatch"):
                owner.execute_gate(ROOT, self.outdir, failed, SHA)
            self.assertEqual(cleanups, ["analysis"])
        old = self.attempt_files(failed)
        self.assertEqual(len([p for p in old if p.suffix == ".dat"]), 2)
        statuses = [json.loads(raw) for path, raw in old.items() if path.name.endswith("-gate-status.json")]
        self.assertEqual(statuses[0]["status"], "failed")
        with ExitStack() as stack:
            self.attempt_harness(stack)
            fresh = self.root / "fresh.zip"
            owner.execute_gate(ROOT, self.outdir, fresh, SHA)
            self.assertTrue(fresh.is_file())
        for path, content in old.items():
            self.assertEqual(path.read_bytes(), content)
        self.assertFalse(failed.exists())

    def test_same_attempt_collisions_reject_before_mutation(self):
        for key in ("run_summary", "yield_data", "status", "zip", "companion"):
            with self.subTest(key=key):
                output = self.root / ("collision-" + key + ".zip")
                profile = owner.resolved_profile(ROOT, SHA, output)
                if key == "zip":
                    path = output
                elif key == "companion":
                    path = output.with_name(output.stem + "-gate-status.json")
                elif key == "status":
                    path = self.outdir / (owner.attempt_identity(output)["attempt_id"] + "-gate-status.json")
                else:
                    name = owner.artifact_name(profile, key, **({} if key == "run_summary" else {"phi": "Left", "epsilon": "lowe"}))
                    path = self.outdir / name
                path.write_bytes(b"existing failed-attempt evidence")
                before = {p: p.read_bytes() for p in self.outdir.iterdir()}
                with patch.object(owner, "source_preflight") as source, patch.object(owner, "run_analysis") as run:
                    with self.assertRaisesRegex(Exception, "attempt_evidence_already_exists"):
                        owner.execute_gate(ROOT, self.outdir, output, SHA)
                    source.assert_not_called()
                    run.assert_not_called()
                self.assertEqual(path.read_bytes(), b"existing failed-attempt evidence")
                self.assertEqual({p: p.read_bytes() for p in self.outdir.iterdir()}, before)

    def test_attempt_identity_validation_and_same_stem_distinct_roots(self):
        first, second = self.root / "one/run.zip", self.root / "two/run.zip"
        self.assertNotEqual(owner.attempt_identity(first), owner.attempt_identity(second))
        self.assertEqual(owner.attempt_identity(first), owner.attempt_identity(first))
        for bad in (Path("relative.zip"), self.root / "bad name.zip", self.root / "évidence.zip",
                    self.root / ("a" * 81 + ".zip"), self.root / "wrong.json", self.root / "one/../run.zip"):
            with self.subTest(bad=bad), self.assertRaisesRegex(Exception, "identity_invalid"):
                owner.attempt_identity(bad)

    def test_accepted_zero_missing_wide_mm_snapshots_remains_fail_closed(self):
        from testing import test_e8_2_baseline_stage_audit as source_tests
        calculate = source_tests._load_calculate_yield_module()
        cut, wide = source_tests._single_bin_histogram(4.0), source_tests._single_bin_histogram(4.0)
        policy = {"parent_id": "frozen-zero", "fit_accepted": False, "action": "zero", "reason": "accepted frozen zero"}
        with patch.object(calculate, "clone_reset_hist", side_effect=lambda *_args: source_tests._single_bin_histogram(0.0)):
            zero = calculate._apply_zero_parent_pion_subtraction_for_bin(
                {"H_MM_DATA_0_0": cut, "H_MM_nosub_DATA_0_0": wide}, 0, 0, policy)
        self.assertTrue(zero["accepted"] and zero["child_valid"])
        self.assertEqual(zero["application_action"], "zero")
        stage = calculate._build_e8_2_pion_stage_capture(zero)
        self.assertFalse(stage["available"])
        self.assertEqual(stage["pion_application_status"], "accepted")
        self.assertEqual(stage["reason"], "accepted_pion_application_missing_exact_wide_object")
        case = source_tests.E82BaselineStageAuditTests()
        case.calculate_yield = calculate
        source = case._source_from_producer(pion_stage=stage)
        self.assertFalse(source["available"])
        self.assertEqual(source["reason"], "e8_2_authority_failure_accepted_pion_application_missing_exact_wide_object")
        self.assertFalse(source_tests.plots.build_full_background_subtraction_e8_2_payload(source)["available"])
        self.assertEqual((cut.Integral(), wide.Integral()), (4.0, 4.0))
        fixture(self.worktree, self.outdir)
        hist = json.loads((self.outdir / "kaon_FullAnalysis_Q4p4W2p74_lowe.json").read_text())["histlist"][0]
        manifest = json.loads((self.outdir / "Center_kaon_rand_sub_Q4p4W2p74_lowe_full-background-subtraction-manifest.json").read_text())
        manifest["pages"][0]["invalid_unavailable_children"] = [{"phi_index": k, "reason": stage["reason"]} for k in range(9)]
        with self.assertRaisesRegex(Exception, "invalid_populated_children"):
            owner.verify_pages(manifest, "Center", "lowe", hist, ["zero"] * 3)


if __name__ == "__main__":
    unittest.main()
