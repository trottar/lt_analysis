"""Deterministic scalar fixtures only; no private artifacts or ROOT imports."""
from __future__ import annotations

import copy
from contextlib import redirect_stdout, redirect_stderr
import hashlib
import io
import json
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest

from testing import analyze_pion_transfer_validity as validator


def encoded_hash(value):
    return hashlib.sha256((json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()).hexdigest()


def hist(values, edges, variance=None):
    errors = [math.sqrt(v) for v in variance] if variance else [0.] * len(values)
    raw = [float(len(edges)-1), float(edges[0]), float(edges[-1]), *map(float, edges)]
    raw += [x for pair in zip(values, errors) for x in pair]
    return dict(status="available", values=list(values), errors=errors, edges=list(edges),
                stored_sumw2=list(variance) if variance else [0.] * len(values),
                axis_signature=encoded_hash(edges), content_error_fingerprint=hashlib.sha256(
                    struct.pack("!" + "d" * len(raw), *raw)).hexdigest())


def seal(doc):
    audit = doc.get("audit", doc)
    audit["collection_identity"] = encoded_hash([p.get("identity", p) for p in audit["parents"]])
    return doc


def fixture(*, model=(2., 3., 4.), kaon=(1., 2., 3.), edges=(-4., 0., 3., 10.),
            contributions=None, floor=1e-12, clip_min=0., clip_max=None, scale=1.,
            zero_variance=False, flows=(0., 0.), parent_count=1, known_head=False):
    """Four signed source classes; scale event weights, yields and variances.

    Source count is unchanged by weight-unit rescaling. Model and floor have
    the same yield units. Bin labels and setting tokens carry no numerical role.
    """
    n = len(model)
    if contributions is None:
        contributions = ([10., 14., 12.], [-9., -4., -3.], [-4., -4., -5.], [1., 1., 1.])
    p, k = ([0., *[v * scale for v in seq], 0.] for seq in (model, kaon))
    weights = [0.]
    for den, num in zip(p[1:-1], k[1:-1]):
        w = num / den if den > floor * scale else 0.
        if clip_min is not None:
            w = max(w, clip_min)
        if clip_max is not None:
            w = min(w, clip_max)
        weights.append(w)
    weights.append(0.)
    cells = {}
    for s, values in zip(validator.SOURCES, contributions):
        signed = [flows[0] * scale if s == "prompt" else 0., *[v * scale for v in values],
                  flows[1] * scale if s == "prompt" else 0.]
        cells[s] = {"signed": signed, "absolute": list(map(abs, signed)),
                    "sumsq": [0. if zero_variance else v*v for v in signed],
                    "count": [int(v != 0) for v in signed]}
    def spectra(weighted=False):
        src = copy.deepcopy(cells)
        if weighted:
            for cell in src.values():
                for i, w in enumerate(weights):
                    cell["signed"][i] *= w
                    cell["absolute"][i] *= abs(w)
                    cell["sumsq"][i] *= w*w
        combined = {key: [math.fsum(cell[key][i] for cell in src.values()) for i in range(n+2)]
                    for key in validator.MOMENTS}
        return dict(status="available", edges=list(edges), selection="nommcuts", sources=src,
                    combined=combined, invalid_records=[], records_fingerprint="frozen-records")
    source, weighted = spectra(), spectra(True)
    control = hist(source["combined"]["signed"], edges, source["combined"]["sumsq"])
    closure = {key: {"status": "PASS"} for key in ("content", "error_squared", "stored_sumw2")}
    parent = dict(status="available", t_bin_index="native-A", identity_closure="PASS",
        identity=dict(pion_parent_id="p-A", t_bin_index="native-A", t_edges=[.12, .78],
            coordinate_fingerprint="coordinate", pion_control_input_fingerprint=control["content_error_fingerprint"]),
        control=control, source_spectra=source, source_closure=copy.deepcopy(closure),
        cache_control_closure={"status": "PASS"}, fit_input_closure={"status": "PASS"},
        fit=dict(histograms={"H_simc_shape_"+c: hist(v, edges) for c, v in
                 zip(validator.COMPONENTS, (p, k, [0.]*(n+2)))},
            diagnostics=dict(pion=dict(fit_min=edges[0], fit_max=edges[-1],
                anchor_windows={"arbitrary": [edges[0], edges[1]]},
                joint_refinement=dict(success=True, chi2_data=250., ndf=2, chi2_ndf=125.,
                    bound_hit_flags={"pi_n": True, "pi_sidis": False}, uncertainties={"pi_n": .2})),
                kaon={}), fit_validation_pion={"accepted": True}),
        transfer=dict(histograms={name: hist(v, edges) for name, v in
                      zip(("H_pion_control_model", "H_kaon_pion_model", "H_pion_weight_vs_MM"), (p, k, weights))},
            weights=weights, denominator_floor=floor*scale, clip_min=clip_min, clip_max=clip_max,
            selected_amplitudes={side: dict(zip(validator.COMPONENTS, a)) for side, a in
                                 (("pion", (1., 0., 0.)), ("kaon", (0., 1., 0.)))},
            closure={key: {"status": "PASS"} for key in ("P_components", "K_components", "model_algebra", "weight_array_vs_hist")}),
        final=dict(histograms={"H_pion_subtraction_template_MM_nosub":
                   hist(weighted["combined"]["signed"], edges, weighted["combined"]["sumsq"])},
            weight_reference=hist(p, edges), weights=list(weights), nommcuts_weighted_sources=weighted,
            nommcuts_event_reproduction=copy.deepcopy(closure), nommcuts_after_closure={"status": "PASS"}),
        proposal={"never": "use this as applied data"})
    audit = dict(schema_version=validator.SCHEMA, production_objects_mutated=False, diagnostic_only=True,
        scope="Q8p7W9p2/unfamiliar-epsilon/oblique", observed_source_HEAD="a"*40 if known_head else "not_provided",
        cache_identity={"coordinate_fingerprint": "coordinate"}, parents=[],
        capture_completeness={"status": "captured"})
    for i in range(parent_count):
        entry = copy.deepcopy(parent)
        entry["t_bin_index"] = entry["identity"]["t_bin_index"] = "native-"+str(i)
        entry["identity"]["pion_parent_id"] = "parent-"+str(i)
        audit["parents"].append(entry)
    return seal(dict(schema_version=validator.SCHEMA, status="captured_pending_independent_validation", audit=audit))


class TransferTests(unittest.TestCase):
    def parent(self, doc=None):
        report = validator.analyze(doc or fixture())
        self.assertEqual(report["arithmetic_status"], "PASS", report)
        return report["parents"][0]

    def assert_blocked(self, doc, reason):
        report = validator.analyze(doc)
        self.assertEqual(report["arithmetic_status"], "FAIL")
        self.assertIn(reason, report["first_failed_invariant"])
        self.assertTrue(all(not p["regions"] for p in report["parents"]))

    def test_exact_identity_and_read_only(self):
        doc = fixture(contributions=([2., 3., 4.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.]))
        before = copy.deepcopy(doc)
        normal = self.parent(doc)["regions"]["normal"]
        self.assertEqual(normal["direct_all_normal_bins"], {"D_R": 0., "L_R": 0.})
        self.assertEqual(normal["qualified"]["conditional_pull"], 0.)
        self.assertEqual(doc, before)

    def test_model_algebra_pass_does_not_accept_event_model(self):
        p = self.parent()
        r = p["regions"]["normal"]
        self.assertNotEqual(r["direct_all_normal_bins"]["D_R"], 0.)
        self.assertAlmostEqual(r["direct_all_normal_bins"]["D_R"], r["qualified"]["reconstructed"]["D_R"])
        self.assertEqual(p["scientific_agreement"], "NOT VERIFIED")
        self.assertTrue(p["fit"]["producer_validation_flags"]["fit_validation_pion"]["accepted"])
        j = p["fit"]["pion"]["joint_refinement"]
        self.assertEqual(j["chi2_ndf"], 125.)
        self.assertEqual(j["boundary_hit_count"], 1)
        self.assertEqual(j["uncertainties"], {"pi_n": .2})

    def test_zero_and_floor_denominator_positive_numerator(self):
        p = self.parent(fixture(model=(0., 1e-12, 4.)))
        self.assertEqual(p["counts"]["unsupported_positive_numerator"], 2)
        self.assertEqual(p["counts"]["denominator_floor_touch"], 2)
        r = p["regions"]["normal"]
        self.assertEqual(r["qualified_bins"], [3])
        self.assertEqual(r["excluded_bins"], [1, 2])
        self.assertNotEqual(r["direct_all_normal_bins"], r["qualified"]["direct"])

    def test_tiny_supported_denominator_high_leverage(self):
        p = self.parent(fixture(model=(2e-12, 3., 4.)))
        self.assertEqual(p["bins"][1]["weight"], 5e11)
        r = p["regions"]["normal"]
        self.assertGreater(r["qualified"]["conditional_fixed_weight_sigma"], 1e12)

    def test_negative_random_dummy_and_effective_support(self):
        p = self.parent()
        region = p["regions"]["fit.diagnostics.pion.anchor_windows.arbitrary"]
        self.assertEqual(region["bins"], [1])
        self.assertEqual(region["source_support"]["signed"], -2.)
        self.assertEqual(region["source_support"]["absolute"], 24.)
        self.assertEqual(region["source_support"]["sumsq"], 198.)
        self.assertAlmostEqual(region["source_support"]["cancellation_fraction"], 1-2/24)
        self.assertAlmostEqual(region["source_support"]["signed_effective_support"], 4/198)
        self.assertAlmostEqual(region["source_support"]["absolute_effective_support"], 576/198)
        self.assertEqual(region["by_source"]["rand"]["signed"], -9.)
        self.assertEqual(region["by_source"]["dummy"]["signed"], -4.)

    def test_nonzero_weights_with_zero_sumw2_are_rejected(self):
        self.assert_blocked(fixture(zero_variance=True), "nonzero weight with zero Sumw2")
        # Positive combined variance must not hide one impossible source cell.
        d = fixture(); cell = d["audit"]["parents"][0]["source_spectra"]
        cell["sources"]["dummy"]["sumsq"][1] = 0.
        cell["combined"]["sumsq"][1] -= 16.
        self.assert_blocked(d, "nonzero weight with zero Sumw2")
        # Signed cancellation does not make nonzero absolute weights admissible.
        cell = dict(status="available", edges=[0., 1.], sources={s: {
            "signed": [0., 0., 0.], "absolute": [0., 2. if s == "prompt" else 0., 0.],
            "sumsq": [0., 0., 0.], "count": [0., 2 if s == "prompt" else 0, 0]}
            for s in validator.SOURCES})
        cell["combined"] = copy.deepcopy(cell["sources"]["prompt"])
        with self.assertRaisesRegex(validator.EvidenceError, "nonzero weight with zero Sumw2"):
            validator.moments(cell, cell["edges"], "cancellation")

    def test_zero_weight_zero_sumw2_remains_admissible_without_significance(self):
        zeros = [0., 0., 0.]
        d = fixture(contributions=(zeros, zeros, zeros, zeros))
        # Counted events may have exactly zero coefficients.
        source = d["audit"]["parents"][0]["source_spectra"]
        weighted = d["audit"]["parents"][0]["final"]["nommcuts_weighted_sources"]
        for spectra in (source, weighted):
            spectra["sources"]["prompt"]["count"][1] = 1
            spectra["combined"]["count"][1] = 1
        p = self.parent(d)
        r = p["regions"]["normal"]
        self.assertIsNone(r["qualified"]["conditional_pull"])
        self.assertEqual(r["qualified"]["uncertainty_status"], "NOT VERIFIED")
        self.assertIsNone(r["source_support"]["signed_effective_support"])

    def test_fit_input_closure_must_be_explicit_pass(self):
        for closure in (None, {}, {"status": "unavailable"}, {"status": "FAIL"},
                        {"status": "NOT VERIFIED"}, "PASS", True, []):
            with self.subTest(closure=closure):
                d = fixture(parent_count=2)
                d["audit"]["parents"][1]["fit_input_closure"] = closure
                self.assert_blocked(d, "fit-input closure")
        d = fixture(); d["audit"]["parents"][0].pop("fit_input_closure")
        self.assert_blocked(d, "fit-input closure")

    def test_source_record_fingerprints_must_be_present_nonempty_strings(self):
        invalid = (None, "", "  ", 0, 1, True, [], ["records"], {}, {"status": "unavailable"})
        for value in invalid:
            for target in ("source", "final", "both"):
                with self.subTest(value=value, target=target):
                    d = fixture(); p = d["audit"]["parents"][0]
                    if target in ("source", "both"):
                        p["source_spectra"]["records_fingerprint"] = copy.deepcopy(value)
                    if target in ("final", "both"):
                        p["final"]["nommcuts_weighted_sources"]["records_fingerprint"] = copy.deepcopy(value)
                    self.assert_blocked(d, "source record provenance")
        for target in ("source", "final", "both"):
            d = fixture(); p = d["audit"]["parents"][0]
            if target in ("source", "both"):
                p["source_spectra"].pop("records_fingerprint")
            if target in ("final", "both"):
                p["final"]["nommcuts_weighted_sources"].pop("records_fingerprint")
            self.assert_blocked(d, "source record provenance")

    def test_clipping_excludes_identity_bins(self):
        p = self.parent(fixture(clip_max=.6))
        self.assertEqual(p["counts"]["clipped"], 2)
        self.assertEqual(p["regions"]["normal"]["qualified_bins"], [1])
        p = self.parent(fixture(model=(0., 3., 4.), clip_min=.2))
        self.assertEqual(p["counts"]["clipped"], 1)
        self.assertTrue(p["bins"][1]["unsupported_positive_numerator"])

    def test_all_unsupported_and_empty_window(self):
        p = self.parent(fixture(model=(0., 0., 0.)))
        self.assertEqual(p["regions"]["normal"]["qualified"]["status"], "NOT VERIFIED")
        doc = fixture()
        doc["audit"]["parents"][0]["fit"]["diagnostics"]["pion"]["anchor_windows"] = {"empty": [100., 101.]}
        r = self.parent(doc)["regions"]["fit.diagnostics.pion.anchor_windows.empty"]
        self.assertEqual(r["bins"], [])
        self.assertIsNone(r["qualified"]["conditional_pull"])

    def test_serialized_none_bound_stays_unknown_without_blocking_unchanged_ratios(self):
        d = fixture(); p = d["audit"]["parents"][0]
        p["transfer"]["clip_max"] = {"status": "unavailable", "reason": "not retained by producer"}
        result = self.parent(d)
        self.assertEqual(result["numerical_policy"]["configured_clipping_status"], "NOT VERIFIED")
        self.assertEqual(result["counts"]["identity_qualified"], 3)
        p["transfer"]["clip_min"] = 2.
        self.assert_blocked(d, "known clipping bound")
        d = fixture(clip_max=.6); p = d["audit"]["parents"][0]
        p["transfer"]["clip_max"] = {"status": "unavailable"}
        self.assert_blocked(d, "without retained clipping evidence")
        p["transfer"]["bins"] = [{"bin": i, "configured_clipping": True, "stored_weight_policy_match": True} for i in (2, 3)]
        self.assertEqual(self.parent(d)["counts"]["clipped"], 2)

    def test_unsigned_discrepancy_survives_signed_region_cancellation(self):
        r = self.parent(fixture(model=(2., 2., 2.), kaon=(1., 1., 1.),
            contributions=([4., 0., 2.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.])))["regions"]["normal"]
        self.assertEqual(r["direct_all_normal_bins"], {"D_R": 0., "L_R": 2.})
        self.assertEqual(r["qualified"]["reconstructed"], {"D_R": 0., "L_R": 2.})

    def test_flow_support_is_preserved_and_excluded(self):
        p = self.parent(fixture(flows=(17., 23.)))
        self.assertEqual(p["flows"]["underflow"]["P_obs"], 17.)
        self.assertEqual(p["flows"]["overflow"]["P_obs"], 23.)
        self.assertEqual(p["flows"]["underflow"]["C_applied"], 0.)
        self.assertEqual(p["regions"]["normal"]["source_support"]["signed"], 10.)

    def test_small_nonzero_kaon_model_retains_absolute_magnitude(self):
        r = self.parent(fixture(kaon=(1e-30, 2e-30, 3e-30)))["regions"]["normal"]
        self.assertGreater(r["sum_abs_K_model"], 0.)
        self.assertLess(r["direct_all_normal_bins"]["L_R"], 1e-28)
        self.assertGreater(r["normalized_absolute_discrepancy"], 0.)

    def test_zero_model_normalization_does_not_divide(self):
        r = self.parent(fixture(kaon=(0., 0., 0.)))["regions"]["normal"]
        self.assertIsNone(r["normalized_absolute_discrepancy"])

    def test_arbitrary_parent_count_labels_axis_and_metadata(self):
        for count in (1, 2, 4, 5):
            with self.subTest(count=count):
                doc = fixture(parent_count=count, edges=(20., 21., 25., 40.))
                report = validator.analyze(doc, {"native_subset": [1, 3]})
                self.assertEqual(report["arithmetic_status"], "PASS", report)
                self.assertEqual(len(report["parents"]), count)
                p = report["parents"][0]
                self.assertEqual(p["regions"]["native_subset"]["bins"], [1, 3])
                self.assertEqual(p["regions"]["fit.diagnostics.pion.anchor_windows.arbitrary"]["window"], [20., 21.])
        self.parent(fixture(model=(2., 3.), kaon=(1., 2.), edges=(20., 26., 50.),
                            contributions=([3., 4.], [-1., -2.], [-1., -1.], [1., 1.])))

    def test_weight_unit_rescaling_preserves_dimensionless_diagnostics(self):
        baseline = self.parent()["regions"]["normal"]
        for scale in (1e-18, 1e-5, 7., 1e18):
            with self.subTest(scale=scale):
                r = self.parent(fixture(scale=scale))["regions"]["normal"]
                for key in ("normalized_signed_discrepancy", "normalized_absolute_discrepancy"):
                    self.assertAlmostEqual(r[key], baseline[key])
                for key in ("signed_effective_support", "absolute_effective_support", "cancellation_fraction"):
                    self.assertAlmostEqual(r["source_support"][key], baseline["source_support"][key])
                self.assertAlmostEqual(r["qualified"]["conditional_pull"], baseline["qualified"]["conditional_pull"])
                self.assertAlmostEqual(r["direct_all_normal_bins"]["L_R"]/scale, baseline["direct_all_normal_bins"]["L_R"])

    def test_bad_axes_fingerprints_and_provenance_stop_interpretation(self):
        mutations = [
            (lambda d,p: p["control"].pop("edges"), "axis"),
            (lambda d,p: p["transfer"]["histograms"]["H_kaon_pion_model"].update(edges=[-5., 0., 3., 10.]), "axis mismatch"),
            (lambda d,p: p["control"].update(content_error_fingerprint="wrong"), "fingerprint"),
            (lambda d,p: p["identity"].update(coordinate_fingerprint="wrong"), "coordinate identity"),
            (lambda d,p: p["identity"].pop("pion_control_input_fingerprint"), "control input identity"),
            (lambda d,p: p.update(identity_closure="FAIL"), "identity closure"),
            (lambda d,p: p["final"]["nommcuts_weighted_sources"].update(records_fingerprint="other"), "record provenance"),
            (lambda d,p: p["source_closure"]["content"].update(status="FAIL"), "retained invariant"),
            (lambda d,p: p["final"].pop("nommcuts_event_reproduction"), "independent final event reproduction"),
        ]
        for mutate, reason in mutations:
            with self.subTest(reason=reason):
                d = fixture(parent_count=2)
                mutate(d, d["audit"]["parents"][1])
                self.assert_blocked(seal(d), reason)
        d = fixture(); d["audit"]["collection_identity"] = "wrong"
        self.assert_blocked(d, "collection identity")

    def test_missing_or_invalid_core_support_blocks(self):
        for key in ("sumsq", "absolute", "count"):
            d = fixture(); d["audit"]["parents"][0]["source_spectra"]["combined"].pop(key)
            self.assert_blocked(d, key)
        for invalid in (-1., float("nan"), float("inf")):
            d = fixture(); d["audit"]["parents"][0]["source_spectra"]["combined"]["sumsq"][1] = invalid
            self.assert_blocked(d, "sumsq")

    def test_numerical_policy_and_component_algebra_are_independently_checked(self):
        d = fixture(); p = d["audit"]["parents"][0]
        p["transfer"]["selected_amplitudes"]["pion"]["pi_n"] = 2.
        self.assert_blocked(d, "component-model algebra")
        d = fixture(); p = d["audit"]["parents"][0]
        p["transfer"]["clip_max"] = .1
        self.assert_blocked(d, "clipping bound")

    def test_final_fallback_cannot_use_proposal(self):
        d = fixture(); d["audit"]["parents"][0]["final"]["weights"][1] = 0.
        self.assert_blocked(d, "final/transfer weights")

    def test_corrupted_final_even_with_saved_pass_is_rejected(self):
        d = fixture(); p = d["audit"]["parents"][0]
        h = p["final"]["histograms"]["H_pion_subtraction_template_MM_nosub"]
        h["values"][1] += 1.
        p["final"]["histograms"]["H_pion_subtraction_template_MM_nosub"] = hist(h["values"], h["edges"], h["stored_sumw2"])
        self.assert_blocked(d, "final reproduction")

    def test_missing_dimensions_and_source_head_remain_unknown(self):
        d = fixture(); p = d["audit"]["parents"][0]
        p["fit"]["diagnostics"] = {}
        r = validator.analyze(d)
        self.assertEqual(r["arithmetic_status"], "PASS")
        self.assertEqual(r["producer_provenance"]["source_identity_status"], "NOT VERIFIED")
        self.assertEqual(r["parents"][0]["metadata_window_status"], "NOT VERIFIED")
        self.assertEqual(r["parents"][0]["fit"]["pion"]["joint_refinement"]["amplitude_uncertainty_status"], "NOT VERIFIED")
        self.assertEqual(r["parents"][0]["simc_raw_occupancy"]["status"], "NOT VERIFIED")
        self.assertEqual(r["scientific_agreement"], "NOT VERIFIED")

    def test_schema_duplicate_parent_and_mutation_claim_rejected(self):
        d = fixture(); d["schema_version"] = "unknown"
        self.assert_blocked(d, "schema")
        d = fixture(parent_count=2); d["audit"]["parents"][1] = copy.deepcopy(d["audit"]["parents"][0])
        self.assert_blocked(seal(d), "duplicate parent")
        d = fixture(); d["audit"]["production_objects_mutated"] = True
        self.assert_blocked(d, "read-only diagnostic provenance")

    def test_regions_cannot_include_flows_or_invent_coordinates(self):
        for bins in ([0], [4], [1, 1], []):
            r = validator.analyze(fixture(), {"bad": bins})
            self.assertEqual(r["arithmetic_status"], "FAIL")
            self.assertFalse(r["parents"][0]["regions"])

    def test_cli_input_unchanged_reports_and_failure_codes(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_path, output = Path(tmp)/"input.json", Path(tmp)/"output.json"
            raw = json.dumps(fixture()).encode()
            input_path.write_bytes(raw)
            with redirect_stdout(io.StringIO()) as out:
                self.assertEqual(validator.main([str(input_path), "--report", str(output), "--bins", "custom=2"]), 0)
            self.assertEqual(input_path.read_bytes(), raw)
            report = json.loads(out.getvalue())
            self.assertEqual(report, json.loads(output.read_text()))
            self.assertEqual(report["input_sha256"], hashlib.sha256(raw).hexdigest())
            saved = output.read_bytes()
            with redirect_stderr(io.StringIO()):
                self.assertEqual(validator.main([str(input_path), "--report", str(input_path)]), 2)
                self.assertEqual(validator.main([str(input_path), "--report", str(output)]), 2)
            self.assertEqual(output.read_bytes(), saved)
            self.assertEqual(input_path.read_bytes(), raw)
            input_path.write_text('{"schema_version": "wrong"}')
            with redirect_stdout(io.StringIO()):
                self.assertEqual(validator.main([str(input_path)]), 1)
            input_path.write_text('{"x":1,"x":2}')
            with redirect_stderr(io.StringIO()):
                self.assertEqual(validator.main([str(input_path)]), 2)

    def test_standalone_cli_works_without_root_or_repository_imports(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"fixture.json"
            path.write_text(json.dumps(fixture()))
            cli = Path(validator.__file__).resolve()
            proc = subprocess.run([sys.executable, "-I", "-B", str(cli), str(path)], cwd=tmp,
                                  capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(json.loads(proc.stdout)["arithmetic_status"], "PASS")


if __name__ == "__main__":
    unittest.main()
