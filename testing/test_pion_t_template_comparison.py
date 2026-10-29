"""Deterministic scalar/fake-hist tests; no farm/ROOT integration acceptance."""
from __future__ import annotations

import ast
from copy import deepcopy
import json
import math
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "src/cuts"), str(REPO / "src/utility")]
import pion_t_template_comparison as comparison


class HighSelectionTests(unittest.TestCase):
    def test_all_five_selections_and_wrong_modes(self):
        for eps,setting in (('low','Center'),('low','Left'),('high','Center'),('high','Left'),('high','Right')):
            inp=dict(Q2='4p4',W='2p74',EPSSET=eps,ParticleType='kaon')
            self.assertTrue(comparison.selected(inp,setting,{comparison.HIGH_ENV:'1'}))
            self.assertEqual(comparison.selected(inp,setting,{}),eps=='low')
            for key,value in (('Q2','3p0'),('W','3p14'),('ParticleType','pion'),('EPSSET','mid')):
                with self.subTest(key=key,setting=setting):
                    self.assertFalse(comparison.selected({**inp,key:value},setting,{comparison.HIGH_ENV:'1'}))
        self.assertFalse(comparison.selected(dict(Q2='4p4',W='2p74',EPSSET='low',ParticleType='kaon'),'Right',
                                            {comparison.HIGH_ENV:'1'}))


class Hist:
    def __init__(self, values=(0., 2., 6., 0.), errors=(0., 1., 2., 0.), edges=(0., 1., 2.)):
        self.values, self.errors, self.edges = list(values), list(errors), list(edges)
    def GetNbinsX(self): return len(self.edges) - 1
    def GetXaxis(self): return self
    def GetBinLowEdge(self, i): return self.edges[i-1]
    def GetXmin(self): return self.edges[0]
    def GetXmax(self): return self.edges[-1]
    def GetBinContent(self, i): return self.values[i]
    def GetBinError(self, i): return self.errors[i]
    def GetSumw2N(self): return len(self.values)
    def GetSumw2(self): return [x*x for x in self.errors]
    def Clone(self, name): return deepcopy(self)
    def SetDirectory(self, value): pass
    def Sumw2(self): pass
    def SetBinContent(self, i, value): self.values[i] = value
    def SetBinError(self, i, value): self.errors[i] = value
    def SetEntries(self, value): self.entries = value


def provenance():
    return dict(status="available", sha256="a"*64, path="fixture-only.root", tree_name="h666",
                coordinate="missmass; -t; phipq*180/pi",
                weight_expression="getattr(evt, 'iter_weight', 1.0)")


def integrated(records=None, **kwargs):
    if records is None:
        records = [(.5, .5, 10., 2.), (1.5, .5, 110., 3.), (1.5, .5, 110., 3.)]
    return comparison.integrate_raw_templates(records, kwargs.get("t", [0., 1., 2.]),
                                             [0., 90., 180.], [0., 1., 2.],
                                             kwargs.get("provenance", provenance()))


def fixture():
    e = integrated()["parents"][0]
    shape = Hist(e["normalized_values"], [math.sqrt(x) for x in e["normalized_variances"]])
    templates = {c: dict(hist=deepcopy(shape), evidence=deepcopy(e)) for c in comparison.COMPONENTS}
    alignment = dict(source="parent", components={c: dict(total_shift_gev=.01, effective_windows=[[0., 2.]])
                                                  for c in comparison.COMPONENTS})
    fit = dict(pion_component_alignment=alignment, H_pion_control_input=Hist(), H_kaon_nosub_input=Hist(),
               template_mm_offset_data=.01, particle_subtraction_phi_setting="Center",
               resolved_subtraction_config=dict(fit="fixture"), A_n=2., B_n=3.,
               staged_amplitudes_scaled=dict(kaon_nosub=dict(pi_n=1.)),
               diagnostics=dict(kaon=dict(covariance=[[1.]], correlation=[[1.]], bound_hits=[])))
    parent = dict(fit_result=fit, t_bin_index=0, t_edges=[0., 1.], analysis_scope="t_bin1",
                  pion_parent_id="fixture-parent", H_pion_control=Hist(), H_proton_cleaned_final_rf=Hist())
    kwargs = dict(selected_shapes={c: shape for c in comparison.COMPONENTS},
                  global_shapes={c: shape for c in comparison.COMPONENTS}, signal_shape=Hist(),
                  sigma0_shape=Hist(), sigma0_diagnostics={}, dynamic_alignment_enabled=True)
    def fitter(pion, kaon, shapes, config, **options):
        candidate = deepcopy(fit)
        candidate.update(H_pion_control_input=pion, H_kaon_nosub_input=kaon,
                         pion_component_alignment=options["pion_component_alignment"])
        return candidate
    return parent, templates, kwargs, fitter


class RawIntegrationTests(unittest.TestCase):
    def test_unequal_phi_populations_are_weighted_not_equal_phi(self):
        e = integrated()["parents"][0]
        self.assertEqual(e["raw"]["values"], [0., 2., 6., 0.])
        self.assertEqual(e["normalized_values"], [0., .25, .75, 0.])
        self.assertNotEqual(e["normalized_values"][1:-1], [.5, .5])
        self.assertEqual(e["normalization_count"], 1)
        self.assertEqual(e["raw_children"][1]["values"], [0., 0., 6., 0.])
        self.assertEqual(e["raw"]["sumw2"], 22.)
        self.assertAlmostEqual(e["normalized_variances"][2], 18./64.)
        self.assertEqual(e["normalized_integral"], 1.)
        self.assertTrue(e["weighted_child_to_parent_closure"])

    def test_canonical_boundaries_and_phi_coverage(self):
        records = [(x, t, phi, 1.) for x, t, phi in
                   [(.5, 0., 0.), (.5, 1., 90.), (1.5, 2., 180.), (.5, 2.1, 0.),
                    (.5, -1., 0.), (.5, .5, 181.), (.5, 1., 90.)]]
        r = integrated(records)
        self.assertEqual([p["raw"]["count"] for p in r["parents"]], [1, 3])
        self.assertEqual(r["excluded"]["count"], 3)
        self.assertEqual(r["parents"][0]["phi_excluded"]["count"], 1)
        self.assertTrue(r["canonical_selection_closure"])
        self.assertEqual(r["parents"][0]["status"], "unavailable")

    def test_signed_weights_flows_and_sumw2(self):
        records = [(-.1, .5, 0., 2.), (2., .5, 0., 3.), (.5, .5, 0., 4.),
                   (.5, .5, 0., -1.), (1.5, .5, 180., 1.), (1.5, .5, 180., 1.)]
        e = integrated(records)["parents"][0]
        self.assertEqual(e["raw"]["signed"], 10.)
        self.assertEqual(e["raw"]["absolute"], 12.)
        self.assertEqual(e["raw"]["sumw2"], 32.)
        self.assertEqual(e["raw_integral"], 5.)
        self.assertEqual(e["raw"]["variances"], [4., 17., 2., 9.])
        self.assertTrue(e["weighted_child_to_parent_closure"])

    def test_empty_sparse_nonpositive_or_negative_bins_rejected(self):
        for records in ([], [(.5, .5, 0., 1.)], [(.5, .5, 0., -1.)]*3,
                        [(.5, .5, 0., -1.), (1.5, .5, 0., 3.)]*3):
            e = integrated(records)["parents"][0]
            self.assertEqual(e["status"], "unavailable")
            self.assertIsNone(e["normalized_values"])
            self.assertEqual(e["normalization_count"], 0)

    def test_invalid_provenance_and_nonfinite_rejected(self):
        for p in ({}, dict(provenance(), weight_expression="raw count"), dict(provenance(), sha256="")):
            self.assertIn("invalid_source_provenance", integrated(provenance=p)["parents"][0]["rejection_reasons"])
        e = integrated([(.5, .5, 0., 1.)]*3 + [(float("nan"), .5, 0., 1.)])["parents"][0]
        self.assertIn("nonfinite_selected_records", e["rejection_reasons"])

    def test_invalid_edges(self):
        for edges in ([0.], [0., 0.], [1., 0.], [0., float("nan")]):
            with self.assertRaises(ValueError):
                integrated(t=edges)

    def test_deterministic_json(self):
        a, b = integrated(), integrated()
        self.assertEqual(comparison.digest(a), comparison.digest(b))
        self.assertEqual(json.loads(json.dumps(a, allow_nan=False)), a)


class FrozenComparisonTests(unittest.TestCase):
    def run_comparison(self, modify=None, fitter_wrap=None):
        parent, templates, kwargs, fitter = fixture()
        if modify:
            modify(parent, templates, kwargs)
        spy = Mock(side_effect=fitter_wrap(fitter) if fitter_wrap else fitter)
        before = comparison.digest(parent)
        report, candidate = comparison.compare_parent(parent, templates, {"frozen": [1]}, fitter=spy, **kwargs)
        self.assertEqual(before, comparison.digest(parent))
        return report, candidate, spy, parent

    def test_fixed_alignment_identical_inputs_and_original_preserved(self):
        r, candidate, spy, parent = self.run_comparison()
        self.assertEqual(r["status"], "available")
        self.assertTrue(r["authoritative_parent_unchanged"])
        pion, kaon, shapes, cfg = spy.call_args.args
        self.assertIsNot(pion, parent["H_pion_control"])
        self.assertIsNot(kaon, parent["H_proton_cleaned_final_rf"])
        self.assertEqual(comparison.histogram(pion), comparison.histogram(parent["H_pion_control"]))
        self.assertEqual(spy.call_args.kwargs["pion_component_alignment"], parent["fit_result"]["pion_component_alignment"])
        self.assertIsNot(spy.call_args.kwargs["pion_component_alignment"], parent["fit_result"]["pion_component_alignment"])
        self.assertNotIn("alignment_outpath", spy.call_args.kwargs)
        self.assertFalse(r["production_applied"])
        self.assertEqual(r["candidate_B"]["diagnostics"]["kaon"]["covariance"], [[1.]])

    def test_disabled_dynamic_alignment_or_missing_map_never_calls_fitter(self):
        for change in (lambda p,t,k: k.update(dynamic_alignment_enabled=False),
                       lambda p,t,k: p["fit_result"].update(pion_component_alignment=None),
                       lambda p,t,k: p["fit_result"]["pion_component_alignment"].update(components={})):
            r, candidate, spy, _ = self.run_comparison(change)
            self.assertEqual(r["status"], "unavailable")
            spy.assert_not_called()

    def test_mismatched_baseline_inputs_reject_before_fit(self):
        def modify(p, t, k): p["H_pion_control"].errors[1] = 99.
        r, _, spy, _ = self.run_comparison(modify)
        self.assertIn("pion_baseline_input_identity_mismatch", r["rejection_reasons"])
        spy.assert_not_called()

    def test_different_global_source_selection_or_variance_rejected(self):
        for modify in (lambda p,t,k: k["global_shapes"]["pi_n"].values.__setitem__(1, .9),
                       lambda p,t,k: k["global_shapes"]["pi_n"].errors.__setitem__(1, 99.)):
            r, candidate, spy, _ = self.run_comparison(modify)
            self.assertFalse(r["templates"]["pi_n"]["production_global_selection_closure"])
            self.assertEqual(r["status"], "unavailable")
            spy.assert_not_called()

    def test_invalid_template_never_falls_back(self):
        for change in (lambda p,t,k: t.pop("pi_n"),
                       lambda p,t,k: t["pi_n"]["evidence"].update(status="globally_inherited"),
                       lambda p,t,k: t["pi_n"]["evidence"].update(canonical_selection_closure=False),
                       lambda p,t,k: t["pi_n"]["evidence"].update(weighted_child_to_parent_closure=False),
                       lambda p,t,k: t["pi_n"]["evidence"].update(t_edges=[1., 2.]),
                       lambda p,t,k: t["pi_n"]["hist"].values.__setitem__(1, 999.),
                       lambda p,t,k: t["pi_n"]["hist"].errors.__setitem__(1, 999.)):
            r, candidate, spy, _ = self.run_comparison(change)
            self.assertEqual(r["status"], "unavailable")
            self.assertIsNone(candidate)
            spy.assert_not_called()

    def test_changed_candidate_alignment_config_or_inputs_rejected(self):
        for mutate in (lambda c: c["pion_component_alignment"].update(source="rescanned"),
                       lambda c: c["H_pion_control_input"].values.__setitem__(1, 999.),
                       lambda c: c.update(resolved_subtraction_config={})):
            def wrap(fitter):
                def run(*args, **kwargs):
                    candidate = fitter(*args, **kwargs)
                    mutate(candidate)
                    return candidate
                return run
            r, candidate, _, _ = self.run_comparison(fitter_wrap=wrap)
            self.assertEqual(r["status"], "unavailable")
            self.assertIsNone(candidate)

    def test_fitter_failure_explicit_and_inputs_still_frozen(self):
        def wrap(fitter):
            def fail(*args, **kwargs):
                args[0].values[1] = 100.
                raise RuntimeError("fixture solver failure")
            return fail
        r, _, _, _ = self.run_comparison(fitter_wrap=wrap)
        self.assertEqual(r["rejection_reasons"], ["candidate_fit_unavailable:fixture solver failure"])

    def test_deterministic_comparison_summary_serialization(self):
        a = self.run_comparison()[0]
        b = self.run_comparison()[0]
        self.assertEqual(comparison.digest(a), comparison.digest(b))
        self.assertEqual(json.loads(json.dumps(a, allow_nan=False)), a)

    def test_ownership_debug_ledger_preserved(self):
        ownership = comparison.ownership
        original_flag, original_records = ownership._DEBUG_ENABLED, list(ownership._DEBUG_RECORDS)
        try:
            ownership._DEBUG_ENABLED = True
            ownership._DEBUG_RECORDS[:] = [{"production": "retained"}]
            def wrap(fitter):
                def run(*args, **kwargs):
                    ownership.configure_particle_subtraction_root_ownership_debug(False)
                    return fitter(*args, **kwargs)
                return run
            self.run_comparison(fitter_wrap=wrap)
            self.assertTrue(ownership._DEBUG_ENABLED)
            self.assertEqual(ownership._DEBUG_RECORDS, [{"production": "retained"}])
        finally:
            ownership._DEBUG_ENABLED = original_flag
            ownership._DEBUG_RECORDS[:] = original_records


class DefaultAndAdapterTests(unittest.TestCase):
    def test_required_source_inventory_paths_exist(self):
        self.assertIn("src/cuts/data_coordinates.py", comparison.REQUIRED_SOURCE_FILES)
        self.assertNotIn("src/utility/data_coordinates.py", comparison.REQUIRED_SOURCE_FILES)
        for path in comparison.REQUIRED_SOURCE_FILES:
            with self.subTest(path=path):
                self.assertTrue((REPO / path).is_file(), path)
                identity = comparison.file_identity(REPO / path)
                self.assertEqual(identity["status"], "available")
                self.assertEqual(len(identity["sha256"]), 64)

    def test_observer_serializes_separate_setting_artifacts_without_changing_hist(self):
        parent, templates, options, fitter = fixture()
        parents = []
        raw = {c: dict(parents=[], selected={}, excluded={}) for c in comparison.COMPONENTS}
        for i in range(3):
            p = deepcopy(parent)
            p.update(t_bin_index=i, t_edges=[float(i), float(i+1)], analysis_scope="t_bin"+str(i+1))
            parents.append(p)
            for c in comparison.COMPONENTS:
                t = deepcopy(templates[c])
                t["evidence"].update(t_bin_index=i, t_edges=p["t_edges"])
                raw[c]["parents"].append(t)
        hist = dict(phi_setting="Center", _pion_t_bin_parent_results=tuple(parents))
        before = comparison.digest(hist)
        loader = Mock(return_value=raw)
        modules = {
            "pion_component_shapes": types.SimpleNamespace(load_raw_t_pion_comparison_shapes=loader),
            "pion_component_fits": types.SimpleNamespace(build_particle_subtraction_component_result=fitter,
                resolve_scope_component_shapes=lambda *a,**k: options["selected_shapes"],
                resolve_scope_single_shape=lambda *a,**k: Hist()),
        }
        def renderer(path, rows, title): path.write_bytes(b"fixture only; not a ROOT-rendered PDF")
        with tempfile.TemporaryDirectory() as tmp, patch.dict(sys.modules, modules), patch.dict(
                "os.environ", {comparison.ENABLE_ENV: "1"}), patch.object(comparison, "render_comparison_pdf", side_effect=renderer):
            result = comparison.observe_template_comparison(hist,
                dict(Q2="4p4", W="2p74", EPSSET="low", ParticleType="kaon"),
                t_edges=[0.,1.,2.,3.], phi_edges=[0.,90.,180.], component_payload={},
                signal_payload={}, sigma0_payload={}, outpath=tmp, hgcer_cutg="existing-cut",
                original_inputs={"data": __file__})
            self.assertEqual(result["status"], "available")
            self.assertEqual(result["source_identity_rejection_reasons"], [])
            paths = sorted(Path(tmp).iterdir())
            self.assertEqual([p.suffix for p in paths], [".json", ".pdf"])
            stored = json.loads(paths[0].read_text())
            self.assertEqual(stored["schema_version"], comparison.SCHEMA)
            self.assertTrue(all(p["authoritative_parent_unchanged"] for p in stored["parents"]))
            self.assertEqual(len(stored["source_HEAD"]), 40)
            self.assertFalse(stored["scientifically_accepted"])
            self.assertEqual(stored["original_analysis_input_identity"]["data"]["status"], "available")
            self.assertEqual(loader.call_args.kwargs["hgcer_cutg"], "existing-cut")
            original_identity = comparison.file_identity
            for missing_path in comparison.REQUIRED_SOURCE_FILES:
                with self.subTest(missing_required_source=missing_path):
                    def missing_identity(path):
                        if Path(path) == REPO / missing_path:
                            return dict(path=str(path), status="unavailable", reason="fixture missing required source")
                        return original_identity(path)
                    with patch.object(comparison, "file_identity", side_effect=missing_identity):
                        unavailable = comparison.observe_template_comparison(hist,
                            dict(Q2="4p4", W="2p74", EPSSET="low", ParticleType="kaon"),
                            t_edges=[0.,1.,2.,3.], phi_edges=[0.,90.,180.], component_payload={},
                            signal_payload={}, sigma0_payload={}, outpath=tmp, hgcer_cutg="existing-cut")
                    self.assertTrue(all(p["status"] == "available" for p in unavailable["parents"]))
                    self.assertEqual(unavailable["status"], "unavailable")
                    reason = "required_source_identity_unavailable:" + missing_path
                    self.assertEqual(unavailable["source_identity_rejection_reasons"], [reason])
                    persisted = json.loads(paths[0].read_text())
                    self.assertEqual(persisted["status"], "unavailable")
                    self.assertEqual(persisted["source_identity_rejection_reasons"], [reason])
            with patch.object(comparison, "render_comparison_pdf", side_effect=RuntimeError("fixture renderer failure")):
                with self.assertRaisesRegex(RuntimeError, "fixture renderer failure"):
                    comparison.observe_template_comparison(hist,
                        dict(Q2="4p4", W="2p74", EPSSET="low", ParticleType="kaon"),
                        t_edges=[0.,1.,2.,3.], phi_edges=[0.,90.,180.], component_payload={},
                        signal_payload={}, sigma0_payload={}, outpath=tmp, hgcer_cutg="existing-cut")
            failed = json.loads(paths[0].read_text())
            self.assertEqual(failed["status"], "failed")
            self.assertEqual(failed["pdf"]["reason"], "fixture renderer failure")
        self.assertEqual(before, comparison.digest(hist))

    def test_detached_pdf_draw_commands_do_not_touch_sources(self):
        clones = []
        class DrawHist(Hist):
            def Clone(self, name):
                cloned = DrawHist(self.values, self.errors, self.edges); clones.append(cloned); return cloned
            def __getattr__(self, key):
                if key.startswith("Set") or key == "Draw":
                    result = Mock(); self.__dict__[key] = result; return result
                raise AttributeError(key)
            def Add(self, other, scale):
                self.values = [a+scale*b for a,b in zip(self.values, other.values)]
        parent, _, _, _ = fixture()
        fit = parent["fit_result"]
        for side, data_key in (("pion", "H_pion_control_input"), ("kaon", "H_kaon_nosub_input")):
            fit[data_key] = DrawHist()
            for suffix in ("", "_refined"):
                fit["H_{}_fit_total{}".format(side, suffix)] = DrawHist()
                for c in comparison.COMPONENTS:
                    fit["H_{}_fit_{}_scaled{}".format(side, c, suffix)] = DrawHist()
        canvas, pad, directory = Mock(), Mock(), Mock()
        legends, panels = [], []
        def legend(*args):
            value = Mock(); legends.append(value); return value
        def panel(*args):
            value = Mock(); panels.append((args, value)); return value
        root = types.SimpleNamespace(TCanvas=Mock(return_value=canvas), TPad=Mock(side_effect=panel),
                                     TLegend=Mock(side_effect=legend), TPaveText=Mock(),
                                     kBlue=4, kRed=2, gPad=pad, gDirectory=directory)
        before = comparison.digest(fit)
        rows = [(dict(t_bin_index=i, rejection_reasons=[]), fit, deepcopy(fit)) for i in range(3)]
        with patch.dict(sys.modules, {"ROOT": root}):
            comparison.render_comparison_pdf(Path("detached.pdf"), rows, "fixture-only")
        self.assertEqual([c.args[0] for c in canvas.Print.call_args_list],
                         ["detached.pdf["] + ["detached.pdf"]*6 + ["detached.pdf]"])
        self.assertEqual(before, comparison.digest(fit))
        pad.cd.assert_called_once()
        directory.cd.assert_called_once()
        self.assertEqual(len(legends), 6)
        self.assertEqual(len(panels), 30)
        for legend in legends:
            legend.SetNColumns.assert_called_once_with(4)
            legend.SetTextSize.assert_called_once_with(.12)
            labels = [c.args[1] for c in legend.AddEntry.call_args_list]
            self.assertEqual(len(labels), len(set(labels)))
            self.assertIn("same frozen data", labels)
            self.assertIn("A baseline total", labels)
            self.assertIn("B local total", labels)
        canvas.Divide.assert_not_called()
        data_clones = [h for h in clones if "SetMarkerStyle" in h.__dict__]
        self.assertTrue(data_clones)
        for hist in data_clones:
            hist.SetStats.assert_called_once_with(False)
        titles = [call.args[0] for hist in clones for call in
                  hist.__dict__.get("SetTitle", Mock()).call_args_list if "NOT ACCEPTED" in call.args[0]]
        self.assertTrue(titles)
        for title in titles:
            self.assertEqual(title.split(";")[1:], ["MM (GeV)", "weighted content"])
            self.assertIn("fixed alignment, no subtraction, NOT ACCEPTED", title.split(";")[0])
        colors = [call.args[0] for hist in clones for call in hist.__dict__.get("SetLineColor", Mock()).call_args_list]
        self.assertEqual(set(colors), {root.kBlue, root.kRed})

    def test_layout_reserves_dedicated_disjoint_legend_and_data_area(self):
        layout = comparison.COMPARISON_LAYOUT
        boxes = [*layout["models"], *layout["residuals"], layout["legend"]]
        for box in boxes:
            x1, y1, x2, y2 = box
            self.assertTrue(0 <= x1 < x2 <= 1 and 0 <= y1 < y2 <= 1)
        for i, a in enumerate(boxes):
            for b in boxes[i+1:]:
                self.assertTrue(a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])
        self.assertTrue(all((r[2]-r[0])*layout["canvas"][0] >= 750 for r in layout["models"]))
        self.assertTrue(all((r[3]-r[1])*layout["canvas"][1] >= 450 for r in layout["models"]))

    def test_unavailable_candidate_notice_is_in_legend_band_and_globals_restore_on_error(self):
        class DrawHist(Hist):
            def __getattr__(self, key):
                if key.startswith("Set") or key == "Draw": return Mock()
                raise AttributeError(key)
            def Add(self, other, scale):
                self.values = [a+scale*b for a,b in zip(self.values, other.values)]
        fit = {key: DrawHist() for key in ("H_pion_control_input", "H_kaon_nosub_input")}
        for side in ("pion", "kaon"):
            for suffix in ("", "_refined"):
                fit["H_{}_fit_total{}".format(side, suffix)] = DrawHist()
        before = comparison.digest(fit)
        for fail in (False, True):
            canvas, previous_pad, directory = Mock(), Mock(), Mock()
            legend, text_boxes = [], []
            def make_legend(*args):
                obj = Mock(); legend.append(obj); return obj
            def make_text(*args):
                obj = Mock(); text_boxes.append((args, obj)); return obj
            root = types.SimpleNamespace(TCanvas=Mock(return_value=canvas), TPad=Mock(side_effect=lambda *a: Mock()),
                TLegend=Mock(side_effect=make_legend), TPaveText=Mock(side_effect=make_text),
                kBlue=4, kRed=2, gPad=previous_pad, gDirectory=directory)
            if fail: canvas.Print.side_effect = [None, RuntimeError("render failure"), None]
            with patch.dict(sys.modules, {"ROOT": root}):
                if fail:
                    with self.assertRaisesRegex(RuntimeError, "render failure"):
                        comparison.render_comparison_pdf(Path("fixture.pdf"), [(dict(t_bin_index=0,
                            rejection_reasons=["known unavailable reason"]), fit, None)])
                else:
                    comparison.render_comparison_pdf(Path("fixture.pdf"), [(dict(t_bin_index=0,
                        rejection_reasons=["known unavailable reason"]), fit, None)])
            for obj in legend: obj.SetNColumns.assert_called_once_with(2)
            notes = [obj for args, obj in text_boxes if args[:4] == (.50, .04, .99, .96)]
            self.assertTrue(notes)
            for note in notes:
                note.AddText.assert_any_call("B unavailable; no global fallback")
                note.AddText.assert_any_call("known unavailable reason")
            self.assertEqual(comparison.digest(fit), before)
            previous_pad.cd.assert_called_once()
            directory.cd.assert_called_once()
            canvas.Close.assert_called_once()

    def test_default_and_explicit_opt_in(self):
        self.assertFalse(comparison.enabled({}))
        for value in ("0", "true", "yes", ""):
            self.assertFalse(comparison.enabled({comparison.ENABLE_ENV: value}))
        self.assertTrue(comparison.enabled({comparison.ENABLE_ENV: "1"}))
        with patch.dict("os.environ", {}, clear=True):
            self.assertIsNone(comparison.observe_template_comparison(None, None, t_edges=None,
                phi_edges=None, component_payload=None, signal_payload=None, sigma0_payload=None, outpath=None))

    def test_out_of_scope_settings_return_without_import_or_write(self):
        with patch.dict("os.environ", {comparison.ENABLE_ENV: "1"}):
            for inp, setting in ((dict(Q2="3p0", W="2p74", EPSSET="low"), "Center"),
                                 (dict(Q2="4p4", W="2p74", EPSSET="high"), "Left"),
                                 (dict(Q2="4p4", W="2p74", EPSSET="low"), "Right")):
                self.assertIsNone(comparison.observe_template_comparison(dict(phi_setting=setting), inp,
                    t_edges=None, phi_edges=None, component_payload=None, signal_payload=None,
                    sigma0_payload=None, outpath=None))

    def test_hook_is_additive_with_no_existing_pdf_routing_changes(self):
        tree = ast.parse((REPO / "src/cuts/rand_sub.py").read_text())
        opt_in = [n for n in ast.walk(tree) if isinstance(n, ast.If)
                  and comparison.ENABLE_ENV in ast.unparse(n.test)]
        self.assertEqual(len(opt_in), 1)
        calls = [n for n in ast.walk(opt_in[0]) if isinstance(n, ast.Call)]
        diagnostic = [n for n in calls if isinstance(n.func, ast.Name) and n.func.id == "observe_template_comparison"]
        self.assertEqual(len(diagnostic), 1)
        self.assertNotIn("pdf_destinations", ast.unparse(opt_in[0]))
        launcher = (REPO / "run_Prod_Analysis.sh").read_text()
        self.assertIn('LT_ANALYSIS_DEBUG_LEFT_LOW= python3 main.py "${main_args_low[@]}"', launcher)
        self.assertNotIn("env -i", launcher)
        self.assertNotIn(comparison.ENABLE_ENV, launcher)

    def test_runtime_adapter_uses_raw_weight_selection_and_no_production_cache(self):
        # Execute only the additive adapter with narrow fake-ROOT dependencies.
        # No import substitution is used to claim a real fitter/runtime test.
        tree = ast.parse((REPO / "src/utility/pion_component_shapes.py").read_text())
        node = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                    and n.name == "load_raw_t_pion_comparison_shapes")
        def event(mm, t, phi, weight=None, selected=True):
            e = types.SimpleNamespace(missmass=mm, t=t, phipq=phi*math.pi/180., selected=selected,
                                      phgcer_x_det=0., phgcer_y_det=0.)
            if weight is not None: e.iter_weight = weight
            return e
        class Tree(list):
            def GetBranch(self, name): return name != "iter_weight"
            def GetLeaf(self, name): return None
            def GetEntries(self): return len(self)
        events = Tree([event(.5, -.5, 0., 2.), event(1.5, -.5, 110., 3.),
                       event(1.5, -.5, 110., 3.), event(.5, -.5, 0., 100., False),
                       event(.5, -1.5, 0.), event(1.5, -1.5, 180.)])
        f = Mock()
        f.IsZombie.return_value = False
        f.Get.return_value = events
        env = dict(math=math, deepcopy=deepcopy, resolve_simc_pion_component_files=lambda i,s: {c: __file__ for c in comparison.COMPONENTS},
                   resolve_simc_tree_name=lambda i: "h666", BG_OPT_MM_PLOT_NBINS=2,
                   BG_OPT_MM_PLOT_MIN=0., BG_OPT_MM_PLOT_MAX=2., set_val=Mock(),
                   apply_HGCer_hole_cut=lambda *a: types.SimpleNamespace(IsInside=lambda *a: False),
                   ROOT=types.SimpleNamespace(TFile=types.SimpleNamespace(Open=lambda *a: f), gDirectory=None),
                   _KAON_SIMC_REQUIRED_BRANCHES=("missmass", "t"),
                   SIMC_PION_COMPONENT_BACKGROUND_MAP=comparison.COMPONENTS,
                   apply_simc_sub_cuts=lambda e,*a: e.selected,
                   _make_component_hist=lambda *a: Hist(values=(0.,0.,0.,0.)))
        exec(compile(ast.Module(body=[node], type_ignores=[]), "adapter_fixture", "exec"), env)
        cut_state = types.SimpleNamespace(**{k: "fixture" for k in
                    ("W", "Q2", "EPSSET", "POL", "ParticleType", "tmin", "tmax", "cut_poly", "TMIN_RESOLUTION_THRESHOLD")})
        hole = types.SimpleNamespace(IsInside=lambda *a: False, GetN=lambda: 0)
        with patch.dict(sys.modules, {"apply_cuts": cut_state}):
            result = env[node.name](dict(mm_min=0., mm_max=2., ParticleType="kaon", Q2="4p4", W="2p74", EPSSET="low"),
                                    "Center", [0.,1.,2.], [0.,90.,180.], hgcer_cutg=hole)
        for c in comparison.COMPONENTS:
            first, second = result[c]["parents"]
            self.assertEqual(first["hist"].values, [0., .25, .75, 0.])
            self.assertEqual(first["evidence"]["raw"]["sumw2"], 22.)
            self.assertEqual(second["hist"].values, [0., .5, .5, 0.])
            self.assertEqual(second["evidence"]["raw"]["sumw2"], 2.)
        self.assertNotIn("_SIMC_COMPONENT_SHAPE_CACHE", ast.unparse(node))
        env["set_val"].assert_not_called()


if __name__ == "__main__":
    unittest.main()
