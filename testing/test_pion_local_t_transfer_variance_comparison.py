"""ROOT-free deterministic C controls; opt-in accepted evidence via environment.

KAONLT_ACCEPTED_AB_ZIP supplies the exact accepted local archive to run the six
real B controls. Without it these evidence-dependent checks explicitly skip.
"""
from copy import deepcopy
import importlib.util
import json
import os
import re
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import numpy as np
from scipy.optimize import least_squares

SPEC = importlib.util.spec_from_file_location("detached_c", Path(__file__).with_name("pion_local_t_transfer_variance_comparison.py"))
c = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(c)


def synthetic():
    n = 12
    e = np.linspace(.7, 1.3, n+1)
    x = np.linspace(0, 1, n)
    t = np.vstack([np.zeros(3), np.column_stack((.1+x, .2+x*x, 1.2-x)), np.zeros(3)])
    init = np.array([.3, .8, .5])
    y = t@np.array([.8, .2, .5])
    y[3] -= .2
    s2 = np.full(n+2, .01)
    q = {"names": c.COMPONENTS, "edges": e, "t": t, "y": y, "s2": s2,
         "ix": np.arange(1, n+1), "initial": init, "tau": np.array([.3, .7, .8]),
         "k": t@np.array([.1, .2, .1])}
    r = least_squares(lambda b: c.base_residual(q, b), init, bounds=(np.zeros(3), np.full(3, np.inf)), method="trf")
    residual = c.base_residual(q, r.x)
    q["joint"] = {"coefficients": dict(zip(q["names"], r.x)), "data_chi2_contribution": float(residual[:-3]@residual[:-3]),
                  "regularization_contribution": float(residual[-3:]@residual[-3:])}
    return q


class Arithmetic(unittest.TestCase):
    def test_analytic_risk_and_signed_cancellation(self):
        # +2 and -2 cancel in content but contribute 4+4 to Sumw2.
        q = synthetic()
        q["y"][1] = 2-2
        q["s2"][1] = 2**2+(-2)**2
        k, p = np.ones(14), np.full(14, 2.)
        self.assertEqual(c.support_risk(k, p, q["s2"], [1])["V"], 2.)
        self.assertEqual(c.transfer_summary(q, p, k)["signed_histogram_applied"][1], 0.)

    def test_variance_scales_squared_source_coefficients(self):
        counts = np.array([12, 8, 3, 5])
        coefficients = np.array([.2, -.1, -.3, .05])
        s2 = np.array([np.sum(counts*coefficients**2)])
        base = c.support_risk(np.array([2.]), np.array([4.]), s2, [0])["V"]
        scaled = c.support_risk(np.array([2.]), np.array([4.]), s2*9, [0])["V"]
        self.assertAlmostEqual(scaled, base*9)

    def test_unsupported_zero_and_nonfinite_denominator(self):
        for value in (0., -1., float("nan"), float("inf")):
            r = c.support_risk(np.array([1.]), np.array([value]), np.array([2.]), [0])
            self.assertEqual(r["status"], "unavailable")
            self.assertIsNone(r["V"])
            self.assertEqual(r["unsupported_bins"], [0])

    def test_zero_numerator_and_zero_variance(self):
        for k, s2 in ((0., 2.), (1., 0.), (0., 0.)):
            r = c.support_risk(np.array([k]), np.array([0.]), np.array([s2]), [0])
            self.assertEqual(r["V"], 0.)
        q = synthetic()
        tr = c.transfer_summary(q, np.zeros(14), np.zeros(14))
        self.assertEqual(tr["weights"], [0.]*14)

    def test_unsupported_bins_are_never_free_zero_weights(self):
        q = synthetic()
        m = q["t"]@q["initial"]
        m[5] = 0
        tr = c.transfer_summary(q, m, q["k"])
        self.assertEqual(tr["unsupported_bins"], [5])
        self.assertIsNone(tr["weights"][5])
        self.assertEqual(tr["unsupported_MM_intervals"], [[q["edges"][4], q["edges"][5]]])

    def test_monotonic_denominator_sensitivity(self):
        q = synthetic()
        m = q["t"]@q["initial"]
        ix = c.domains(q)["normal_full_MM"]
        a = c.support_risk(q["k"], m, q["s2"], ix)["V"]
        b = c.support_risk(q["k"], m*2, q["s2"], ix)["V"]
        self.assertAlmostEqual(a, 4*b)

    def test_zero_regularization_is_B(self):
        q = synthetic()
        r, control = c.reconstruct_b(q)
        v = c.support_risk(q["k"], q["t"]@r.x, q["s2"], q["ix"])["V"]
        result = c.fit_candidate(q, 0., r, v)
        self.assertEqual(control["status"], "PASS")
        np.testing.assert_array_equal(list(result["coefficients"].values()), r.x)

    def test_analytic_regularized_solution(self):
        q = synthetic()
        ix = q["ix"]
        design = q["t"][ix]/np.sqrt(q["s2"][ix])[:, None]
        precision = np.diag(1/q["tau"]**2)
        analytic = np.linalg.solve(design.T@design+precision,
                                   design.T@(q["y"][ix]/np.sqrt(q["s2"][ix]))+precision@q["initial"])
        self.assertTrue(np.all(analytic>0))
        numerical, _ = c.reconstruct_b(q)
        np.testing.assert_allclose(numerical.x, analytic, atol=1e-7, rtol=0)

    def test_bound_hits_are_reported_without_CI_claim(self):
        q = synthetic()
        r = type("Result", (), {"jac": np.vstack([np.eye(3), np.zeros((12, 3))]),
                                 "x": np.array([0., 1., 2.])})()
        d = c.conditional_diagnostics(q, r)
        self.assertEqual(d["bound_hits"], [True, False, False])
        self.assertIn("especially invalid near bounds", d["limits"])

    def test_risk_reduction_retains_shape_tradeoff_and_fixed_K(self):
        q = synthetic()
        frozen = q["k"].copy()
        r, _ = c.reconstruct_b(q)
        v = c.support_risk(q["k"], q["t"]@r.x, q["s2"], q["ix"])["V"]
        sweep = [c.fit_candidate(q, lam, r, v) for lam in c.LAMBDAS]
        risks = [a["transfer"]["domains"]["normal_full_MM"]["V"] for a in sweep]
        self.assertTrue(all(a>=b-1e-8 for a, b in zip(risks, risks[1:])))
        self.assertTrue(sweep[-1]["shape_fit_worsened"])
        self.assertFalse(any(a["selected"] for a in sweep))
        np.testing.assert_array_equal(q["k"], frozen)

    def test_solver_nonconvergence_recorded(self):
        q = synthetic()
        r, _ = c.reconstruct_b(q)
        result = c.fit_candidate(q, 1000., r, 1., max_nfev=1)
        self.assertEqual(result["status"], "unstable")
        self.assertIn("nonconvergence", result["reason"])
        self.assertIn("last_residual_vector", result)
        with patch.object(c, "least_squares", return_value=type("Result", (), {"success": False, "message": "test"})()):
            with self.assertRaisesRegex(c.InvariantError, "B reconstruction nonconvergence"):
                c.reconstruct_b(q)

    def test_B_control_mismatch_stops(self):
        q = synthetic()
        q["joint"]["coefficients"]["pi_n"] += .01
        with self.assertRaisesRegex(c.InvariantError, "coefficient mismatch"):
            c.reconstruct_b(q)

    def test_conditional_correlation_sensitivity_and_limits(self):
        q = synthetic()
        r, _ = c.reconstruct_b(q)
        d = c.conditional_diagnostics(q, r)
        self.assertEqual(d["status"], "conditional_only")
        np.testing.assert_allclose(np.diag(d["correlation"]), 1.)
        self.assertIn("not CI", d["limits"])
        perturb = deepcopy(q)
        perturb["y"][4] += 1e-5
        rr = least_squares(lambda b: c.base_residual(perturb, b), q["initial"], bounds=(np.zeros(3), np.full(3, np.inf)))
        np.testing.assert_allclose((rr.x-r.x)/1e-5, np.array(d["d_beta_d_data_fit_bins"])[:, 3], atol=2e-3)

    def test_rank_deficient_covariance_is_unavailable(self):
        q = synthetic()
        r = type("Result", (), {"jac": np.ones((10, 3)), "x": np.zeros(3)})()
        self.assertEqual(c.conditional_diagnostics(q, r)["status"], "unavailable")

    def test_no_postfit_weight_clipping(self):
        q = synthetic()
        m = np.full(14, 1e-12)
        tr = c.transfer_summary(q, m, q["k"])
        self.assertGreater(tr["weight_bin_quantiles_unweighted_not_event_quantiles"]["max"], 1e10)

    def test_domain_includes_problematic_MM_and_separate_signal(self):
        q = synthetic()
        d = c.domains(q)
        centers = (q["edges"][:-1]+q["edges"][1:])/2
        self.assertTrue(any(.94<=centers[i-1]<=1.03 for i in d["normal_full_MM"]))
        self.assertTrue(all(1.10<=centers[i-1]<=1.16 for i in d["signal_1.10_1.16"]))

    def test_reproducible_json_and_finite_values(self):
        q = synthetic()
        r, _ = c.reconstruct_b(q)
        a = c.fit_candidate(q, .1, r, 1.)
        b = c.fit_candidate(q, .1, r, 1.)
        self.assertEqual(c.digest(a), c.digest(b))
        self.assertEqual(json.loads(c.encoded(a))["status"], "available")
        with self.assertRaises(ValueError):
            c.encoded({"bad": float("nan")})

    def test_out_of_grid_is_rejected(self):
        q = synthetic()
        r, _ = c.reconstruct_b(q)
        with self.assertRaisesRegex(c.InvariantError, "predeclared"):
            c.fit_candidate(q, 7., r, 1.)

    def test_input_scope_and_output_collisions(self):
        with self.assertRaisesRegex(c.InvariantError, "exact Center/Left"):
            c.load_inputs()
        with self.assertRaisesRegex(c.InvariantError, "both Center"):
            c.load_inputs(center_json="one.json")
        with self.assertRaisesRegex(c.InvariantError, "basename"):
            c.load_inputs(center_json="wrong.json", left_json="also-wrong.json")
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/"result.json"
            with self.assertRaisesRegex(c.InvariantError, "collision"):
                c.main(["--archive", str(path), "--output-json", str(path)])
            path.write_text("preserve")
            with self.assertRaisesRegex(c.InvariantError, "new files"):
                c.main(["--archive", "missing.zip", "--output-json", str(path)])
            self.assertEqual(path.read_text(), "preserve")
        with self.assertRaisesRegex(c.InvariantError, "outside source"):
            c.main(["--archive", "missing.zip", "--output-json", str(Path(__file__).with_name("never-create.json"))])

    def test_bad_archive_provenance_rejected_in_isolation(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/"input.zip"
            with zipfile.ZipFile(path, "w") as z:
                z.writestr("manifest.json", json.dumps({"complete": True, "git_head": "wrong"}))
            with self.assertRaisesRegex(c.InvariantError, "archive provenance"):
                c.load_inputs(archive=path)


@unittest.skipUnless(os.environ.get("KAONLT_ACCEPTED_AB_ZIP"), "accepted A/B archive not supplied; six-parent evidence unexecuted")
class AcceptedEvidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.docs, cls.provenance = c.load_inputs(archive=Path(os.environ["KAONLT_ACCEPTED_AB_ZIP"]))
        cls.before = c.digest(cls.docs)
        cls.result = c.compare(cls.docs, cls.provenance)

    def test_all_six_B_solutions_and_objectives(self):
        self.assertEqual(len(self.result["parents"]), 6)
        for row in self.result["parents"]:
            with self.subTest(setting=row["setting"], t=row["t_bin"]):
                self.assertEqual(row["B_reconstruction"]["status"], "PASS")
                self.assertLessEqual(row["B_reconstruction"]["max_absolute_coefficient_difference"], 1e-7)
                self.assertEqual([x["lambda"] for x in row["C"]], list(c.LAMBDAS))
                self.assertTrue(all(x["status"] == "available" for x in row["C"]))
                joint = row["A_B_unchanged"]["candidate_B"]["diagnostics"]["pion"]["joint_refinement"]
                for name, val in joint["coefficients"].items():
                    self.assertLessEqual(abs(row["C"][0]["coefficients"][name]-val), 1e-7)
                for candidate in row["C"]:
                    self.assertEqual(set(candidate["component_statistical_proxies"]), set(c.COMPONENTS))
                    self.assertLess(candidate["transfer"]["domains"]["normal_full_MM"]["model_ratio_closure_max_absolute"], 1e-12)

    def test_A_B_unchanged_and_schema(self):
        self.assertEqual(c.digest(self.docs), self.before)
        self.assertEqual(self.result["schema_version"], c.SCHEMA)
        self.assertFalse(self.result["production_applied"])
        c.encoded(self.result)  # no NaN/Inf serialization
        for row in self.result["parents"]:
            self.assertEqual(row["A_B_unchanged"], self.docs[row["setting"]]["parents"][row["t_bin"]-1])

    def test_input_fingerprint_provenance_alignment_and_Sumw2_rejections(self):
        for mutation, reason in (
            (lambda d: d.update(source_HEAD="wrong"), "source/schema"),
            (lambda d: d["parents"][0]["frozen_inputs"]["pion"]["values"].__setitem__(3, 99.), "fingerprint"),
            (lambda d: d["parents"][0]["candidate_B"].update(pion_component_alignment={}), "alignment"),
            (lambda d: d["parents"][0]["templates"]["pi_n"].update(status="inherited"), "local template"),
            (lambda d: d["parents"][0]["frozen_inputs"]["pion"].pop("stored_sumw2"), "Sumw2 absent"),
        ):
            with self.subTest(reason=reason):
                d = deepcopy(self.docs["Center"])
                mutation(d)
                with self.assertRaisesRegex(c.InvariantError, reason):
                    c.validate_document(d, "Center")

    def test_aligned_columns_not_raw_templates(self):
        p = self.docs["Center"]["parents"][1]
        q = c.problem(p)
        for i, n in enumerate(q["names"]):
            np.testing.assert_array_equal(q["t"][:, i], p["candidate_B"]["H_simc_shape_"+n]["values"])
        self.assertTrue(any(not np.allclose(q["t"][:, i], p["templates"][n]["normalized_values"])
                            for i, n in enumerate(q["names"])))

    def test_all_six_controls_before_positive_sweep(self):
        docs = deepcopy(self.docs)
        docs["Left"]["parents"][2]["candidate_B"]["diagnostics"]["pion"]["joint_refinement"]["coefficients"]["pi_n"] += .01
        with patch.object(c, "fit_candidate", side_effect=AssertionError("positive fit before controls")):
            with self.assertRaisesRegex(c.InvariantError, "coefficient mismatch"):
                c.compare(docs, self.provenance)

    def test_sweep_stops_after_first_nonconvergence(self):
        with patch.object(c, "fit_candidate", return_value={"lambda": 0., "status": "unstable", "reason": "nonconvergence"}) as fit:
            result = c.compare(self.docs, self.provenance)
        self.assertEqual(fit.call_count, 1)
        self.assertEqual(result["status"], "unstable")
        self.assertIn("first_failed_invariant", result)

    def test_pdf_all_four_panels_and_signal_zoom_render(self):
        try:
            import matplotlib  # optional report dependency
        except ImportError:
            self.skipTest("matplotlib unavailable: local PDF rendering unexecuted")
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/"diagnostic.pdf"
            result = {**self.result, "parents": self.result["parents"][:1]}
            c.render_pdf(result, path)
            raw = path.read_bytes()
            self.assertTrue(raw.startswith(b"%PDF-"))
            self.assertEqual(len(re.findall(rb"/Type /Page\b", raw)), 4)


if __name__ == "__main__":
    unittest.main()
