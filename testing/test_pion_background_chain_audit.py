"""Scalar/fake-ROOT evidence tests; no ROOT runtime or physics acceptance."""
from __future__ import annotations

import ast
import copy
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest.mock import Mock, patch

REPO = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPO / "src/cuts"), str(REPO / "src/utility")]
import pion_background_chain_audit as audit


class Hist:
    def __init__(self, values=(0., 0.), errors=None, edges=(.8, 1.03, 1.4), flows=(0., 0.)):
        self.edges = list(edges)
        self.values = [flows[0], *values, flows[1]]
        self.errors = list(errors) if errors is not None else [0.] * len(self.values)

    def GetNbinsX(self): return len(self.values) - 2
    def GetXaxis(self): return self
    def GetBinLowEdge(self, i): return self.edges[i-1]
    def GetXmin(self): return self.edges[0]
    def GetXmax(self): return self.edges[-1]
    def GetBinContent(self, i): return self.values[i]
    def GetBinError(self, i): return self.errors[i]
    def GetSumw2N(self): return len(self.values)
    def GetSumw2(self): return [e*e for e in self.errors]


def fixture():
    inp = dict(ParticleType="kaon", Q2="4p4", W="2p74", EPSSET="low",
               particle_subtraction_mode="simc_shape_components", pion_subtraction_scope="t_bin")
    hist = dict(phi_setting="Center", _authoritative_pion_control_source_cache=dict(
        coordinate_fingerprint="coord", definition="frozen signed source", by_t=[]))
    parents = []
    shapes = dict(zip(audit.COMPONENTS, ((.1, .9), (.4, .6), (.5, .5))))
    for index in range(3):
        t_edges = [.4 + .2 * index, .6 + .2 * index]
        records = [dict(source_label=s, coefficient=c, adj_MM=.9, adj_t=t_edges[0],
                        t_index=index, coordinate_fingerprint="coord", nommcuts=True, allcuts=False)
                   for s, c in zip(audit.SOURCES, (1., -2., -4., 3.))]
        records += [dict(records[0], coefficient=10., adj_MM=1.1, allcuts=True)]
        control = Hist((-2., 10.), errors=[0., math.sqrt(30), 10., 0.])
        cut = Hist((0., 10.), errors=[0., 0., 10., 0.])
        b = dict(pi_n=.0002021, pi_sidis=2., pi_delta=3.)
        a = dict(pi_n=1., pi_sidis=1., pi_delta=0.)
        p = [sum(b[c]*shapes[c][i] for c in audit.COMPONENTS) for i in range(2)]
        k = [sum(a[c]*shapes[c][i] for c in audit.COMPONENTS) for i in range(2)]
        weights = [0., *[num/den for num, den in zip(k, p)], 0.]
        applied = [x*y for x, y in zip(control.values[1:-1], weights[1:-1])]
        weighted_errors = [0., math.sqrt(30)*weights[1], 10.*weights[2], 0.]
        fit = dict(H_pion_control_input=control, H_kaon_nosub_input=Hist((5., 6.)),
            diagnostics=dict(pion=dict(staged_amplitudes_scaled=dict(b, pi_n=0.), refined_amplitudes=b,
                fit_order=list(audit.COMPONENTS), fit_min=.8, fit_max=1.4,
                chi2=12., ndf=2, chi2_ndf=6., fit_p_value=.002,
                joint_refinement=dict(success=True, message="fixture joint")),
                kaon=dict(staged_amplitudes_scaled=dict(a, pi_delta=2.), refined_amplitudes=a,
                    pi_delta_signal_protected_fit=dict(enabled=True, status="rejected", success=False,
                        reason="protected fixture failure", protected_applied_amplitudes=a))),
            B_n=b['pi_n'], B_sidis=b['pi_sidis'], B_delta=b['pi_delta'],
            A_n=a['pi_n'], A_sidis=a['pi_sidis'], A_delta=a['pi_delta'],
            fit_status_pion=0, fit_status_kaon=0)
        for c in audit.COMPONENTS:
            fit["H_simc_shape_"+c] = Hist(shapes[c])
            for side, amplitudes in (("pion", b), ("kaon", a)):
                fit["H_"+side+"_fit_"+c+"_scaled"] = Hist([amplitudes[c]*v for v in shapes[c]])
        for side in ("pion", "kaon"):
            fit['H_'+side+'_fit_step_overlays'] = [dict(pass_index=1, step_index=i+1, component_name=c,
                amplitude=b[c], H_residual_input=control, H_component_template=Hist(shapes[c]),
                H_component_scaled=Hist((1., 2.)), H_baseline_before=Hist((3., 4.)),
                H_cumulative_after=Hist((4., 6.)), fit_diagnostics=dict(fit_bin_indices=[1],
                    fit_y=[-2.], fit_sigma=[math.sqrt(30)], fit_template=[shapes[c][0]]))
                for i, c in enumerate(audit.COMPONENTS)]
            fit['diagnostics'][side]['staged_pass_history'] = [{c: dict(success=True, message="fixture solve") for c in audit.COMPONENTS}]
        proposal = dict(weights=weights, H_pion_control_model=Hist(p), H_kaon_pion_model=Hist(k),
            H_pion_weight_vs_MM=Hist(weights[1:-1]), H_pion_subtraction_template_MM_nosub=Hist(applied, weighted_errors),
            H_MM_nosub_before_pion_subtraction=Hist((5., 6.)),
            H_MM_nosub_after_pion_subtraction=Hist([5.-applied[0], 6.-applied[1]]),
            H_pion_subtraction_template_MM=Hist((0., applied[1]), [0., 0., weighted_errors[2], 0.]),
            H_MM_before_pion_subtraction=Hist((5., 6.)), H_MM_after_pion_subtraction=Hist((5., 6.-applied[1])),
            diagnostics=dict(denom_floor=1e-12, clip_min=0., clip_max=None, model_variant="final"))
        parent = dict(t_bin_index=index, t_edges=t_edges, pion_parent_id="parent-"+str(index),
            coordinate_fingerprint="coord", H_pion_control=control,
            pion_control_input_fingerprint=audit.fingerprint_histogram_content_error(control), fit_result=fit,
            proposed_diagnostic_application_result=proposal, final_diagnostic_application_result=copy.deepcopy(proposal),
            diagnostic_application_status=dict(final_status="applied_component"))
        parents.append(parent)
        hist['_authoritative_pion_control_source_cache']['by_t'].append(dict(t_index=index,
            t_edges=t_edges, coordinate_fingerprint="coord", H_pion_control=control, H_pion_control_cut=cut, records=records))
    return hist, inp, tuple(parents), {}


def fake_pdf(path, data):
    path.write_bytes(b"%PDF-1.4\nfixture; not a rendered physics document\n%%EOF")
    return len(audit.page_specs(data))


class ChainTests(unittest.TestCase):
    def test_model_algebra_passes_while_event_model_fails(self):
        result = audit.build_chain_audit(*fixture())
        for p in result['parents']:
            self.assertEqual(p['identity_closure'], 'PASS')
            self.assertEqual(p['source_closure']['content']['status'], 'PASS')
            self.assertEqual(p['source_closure']['error_squared']['status'], 'PASS')
            self.assertEqual(p['transfer']['closure']['model_algebra']['status'], 'PASS')
            self.assertEqual(p['transfer']['closure']['P_components']['status'], 'PASS')
            self.assertEqual(p['transfer']['closure']['K_components']['status'], 'PASS')
            self.assertEqual(p['proposal']['event_model_closure']['status'], 'FAIL')
            self.assertEqual(p['proposal']['residual_identity']['status'], 'PASS')
            self.assertEqual(p['proposal']['nommcuts_event_reproduction']['content']['status'], 'PASS')
            self.assertEqual(p['final']['nommcuts_after_closure']['status'], 'PASS')
            self.assertEqual(p['final']['allcuts_after_closure']['status'], 'PASS')

    def test_three_components_tiny_amplitude_protected_failure_stays_observed(self):
        p = audit.build_chain_audit(*fixture())['parents'][1]
        self.assertEqual(p['transfer']['selected_amplitudes']['pion']['pi_n'], .0002021)
        self.assertEqual(p['fit']['diagnostics']['pion']['staged_amplitudes_scaled']['pi_n'], 0.)
        self.assertEqual(p['transfer']['selected_amplitudes']['kaon']['pi_delta'], 0.)
        self.assertEqual(p['fit']['diagnostics']['kaon']['pi_delta_signal_protected_fit']['status'], 'rejected')
        self.assertEqual(p['fit']['steps']['pion'][0]['solver_outcome']['success'], True)
        self.assertEqual(p['fit']['steps']['pion'][0]['cumulative_arithmetic']['status'], 'PASS')

    def test_negative_source_and_distinct_cut_selection(self):
        p = audit.build_chain_audit(*fixture())['parents'][0]
        self.assertEqual(p['source_spectra']['combined']['signed'][1], -2.)
        self.assertEqual(p['source_spectra']['combined']['absolute'][1], 10.)
        self.assertEqual(p['source_spectra']['combined']['sumsq'][1], 30.)
        self.assertEqual(p['cut_source_closure']['content']['status'], 'PASS')
        self.assertNotEqual(p['proposal']['nommcuts_weighted_sources'], p['proposal']['allcuts_weighted_sources'])

    def test_mm_edges_flows_and_t_boundary(self):
        data = fixture(); parent = data[2][0]
        records = data[0]['_authoritative_pion_control_source_cache']['by_t'][0]['records']
        base = records[0]
        edge_records = [dict(base, coefficient=1., adj_MM=x) for x in (.79, .8, 1.029, 1.03, 1.399, 1.4)]
        ref = audit.histogram(parent['H_pion_control'])
        source = audit.source_spectra(edge_records, ref, parent)
        self.assertEqual(source['combined']['signed'], [1., 2., 2., 1.])
        weighted = audit.source_spectra(edge_records, ref, parent, weights=[99., 0., 2., 3.])
        self.assertEqual(weighted['combined']['signed'], [0., 0., 4., 3.])
        edge_records.append(dict(base, adj_t=parent['t_edges'][1]))
        self.assertEqual(audit.source_spectra(edge_records, ref, parent)['status'], 'unavailable')

    def test_only_final_canonical_endpoint_is_closed(self):
        data = fixture()
        last = data[0]['_authoritative_pion_control_source_cache']['by_t'][2]
        last['records'][0]['adj_t'] = last['t_edges'][1]
        self.assertEqual(audit.build_chain_audit(*data)['parents'][2]['source_spectra']['status'], 'available')
        interior = data[0]['_authoritative_pion_control_source_cache']['by_t'][0]
        interior['records'][0]['adj_t'] = interior['t_edges'][1]
        self.assertEqual(audit.build_chain_audit(*data)['parents'][0]['source_spectra']['status'], 'unavailable')

    def test_missing_and_mismatched_source_identity(self):
        for mutation in ('coefficient', 'coordinate_fingerprint', 't_index'):
            data = fixture(); record = data[0]['_authoritative_pion_control_source_cache']['by_t'][0]['records'][0]
            record.pop(mutation)
            p = audit.build_chain_audit(*data)['parents'][0]
            self.assertEqual(p['source_closure']['status'], 'unavailable')
        data = fixture(); data[2][0]['pion_control_input_fingerprint'] = 'wrong'
        self.assertEqual(audit.build_chain_audit(*data)['parents'][0]['identity_closure'], 'FAIL')

    def test_original_error_and_cache_mismatch(self):
        data = fixture(); data[2][0]['H_pion_control'].errors[1] += 1.
        p = audit.build_chain_audit(*data)['parents'][0]
        self.assertEqual(p['source_closure']['error_squared']['status'], 'FAIL')
        self.assertEqual(p['identity_closure'], 'FAIL')
        data = fixture(); data[0]['_authoritative_pion_control_source_cache']['coordinate_fingerprint'] = 'wrong'
        self.assertEqual(audit.build_chain_audit(*data)['parents'][0]['identity_closure'], 'FAIL')

    def test_corrupted_actual_template_after_and_cut_are_detected(self):
        for name, closure in [('H_pion_subtraction_template_MM_nosub','nommcuts_event_reproduction'),
                              ('H_MM_nosub_after_pion_subtraction','nommcuts_after_closure'),
                              ('H_MM_after_pion_subtraction','allcuts_after_closure')]:
            data = fixture(); final = data[2][0]['final_diagnostic_application_result']; final[name].values[2] += 1.
            result = audit.build_chain_audit(*data)['parents'][0]['final'][closure]
            self.assertEqual(result.get('status', result.get('content', {}).get('status')), 'FAIL')

    def test_unsupported_and_clip_qualify_identity_only(self):
        data = fixture(); payload = data[2][0]['proposed_diagnostic_application_result']
        payload['H_pion_control_model'].values[1] = 0.
        payload['weights'][1] = payload['H_pion_weight_vs_MM'].values[1] = 0.
        payload['weights'][2] = payload['H_pion_weight_vs_MM'].values[2] = .01
        payload['diagnostics']['clip_max'] = .01
        transfer = audit.build_chain_audit(*data)['parents'][0]['transfer']
        self.assertEqual(transfer['qualified_bins'], [])
        self.assertTrue(transfer['bins'][1]['unsupported_numerator'])
        self.assertFalse(transfer['bins'][2]['unclipped'])
        self.assertEqual(transfer['closure']['model_algebra']['status'], 'unavailable')

    def test_axis_signature_mismatch_not_hidden_by_equal_bin_counts(self):
        data = fixture(); data[2][0]['fit_result']['H_simc_shape_pi_delta'].edges[1] = 1.02
        self.assertEqual(audit.build_chain_audit(*data)['parents'][0]['transfer']['closure']['reason'], 'axis signature mismatch')

    def test_missing_steps_parent_and_nonfinite_strict_json(self):
        data = fixture(); data[2][0]['fit_result']['H_pion_fit_step_overlays'] = []
        data[2][1]['fit_result']['B_n'] = float('nan')
        result = audit.build_chain_audit(data[0], data[1], data[2][:2], data[3])
        self.assertEqual(result['parents'][0]['fit']['steps']['pion']['status'], 'unavailable')
        self.assertEqual(result['parents'][2]['status'], 'unavailable')
        self.assertIn('nonfinite', result['parents'][1]['fit']['amplitudes']['B_n']['reason'])
        json.loads(audit.encoded(result))
        self.assertEqual(sum('t3' in p['title'] for p in audit.page_specs(result)), 4)

    def test_fallback_final_is_not_relabelled_as_component_model(self):
        data = fixture(); final = data[2][0]['final_diagnostic_application_result']
        final.update(fallback_used=True, fallback_mode='single_scale', final_application_status='applied_fallback')
        final['weights'] = [1.] * 4
        final['H_pion_subtraction_template_MM_nosub'] = Hist((-2.,10.), [0.,math.sqrt(30),10.,0.])
        final['H_MM_nosub_after_pion_subtraction'] = Hist((7.,-4.))
        p = audit.build_chain_audit(*data)['parents'][0]
        self.assertEqual(p['final']['policy']['fallback_mode'], 'single_scale')
        self.assertEqual(p['final']['nommcuts_event_reproduction']['content']['status'], 'PASS')
        self.assertEqual(p['final']['nommcuts_after_closure']['status'], 'PASS')
        self.assertNotIn('event_model_closure', p['final'])

    def test_capture_does_not_mutate_any_input(self):
        data = fixture(); before = copy.deepcopy(data)
        first, second = audit.build_chain_audit(*data), audit.build_chain_audit(*data)
        self.assertEqual(audit.encoded(first), audit.encoded(second))
        self.assertEqual(audit.snapshot(data), audit.snapshot(before))

    def test_other_settings_modes_no_output(self):
        for key, value in [('ParticleType','pion'),('Q2','3p0'),('W','3p14'),('EPSSET','high'),
                           ('particle_subtraction_mode','single_scale'),('pion_subtraction_scope','setting')]:
            data = fixture(); data[1][key] = value
            with patch.object(audit,'write_chain_artifacts') as writer:
                self.assertIsNone(audit.observe_chain_audit(*data,'unused'))
                writer.assert_not_called()
        data = fixture(); data[0]['phi_setting'] = 'Left'
        self.assertIsNone(audit.build_chain_audit(*data))

    def test_atomic_publication_repeat_reuse_and_mismatched_collision(self):
        result = audit.build_chain_audit(*fixture())
        with tempfile.TemporaryDirectory() as directory:
            renderer = Mock(side_effect=fake_pdf)
            first = audit.write_chain_artifacts(directory,result,renderer)
            initial = {p.name:p.read_bytes() for p in Path(directory).iterdir()}
            repeat = audit.write_chain_artifacts(directory,result,renderer)
            self.assertTrue(repeat['reused']); self.assertEqual(renderer.call_count,1)
            self.assertEqual(initial,{p.name:p.read_bytes() for p in Path(directory).iterdir()})
            self.assertIn('Q4p4W2p74_Center_lowe', Path(first['json']).name)
            altered = copy.deepcopy(result); altered['limits'].append('collision')
            with self.assertRaises(FileExistsError): audit.write_chain_artifacts(directory,altered,renderer)
            Path(directory, first['pdf']['basename']).write_bytes(b'corrupt')
            with self.assertRaises(FileExistsError): audit.write_chain_artifacts(directory,result,renderer)

    def test_render_failure_publishes_only_explicit_unavailable_json(self):
        result = audit.build_chain_audit(*fixture())
        def fail(path, data):
            path.write_bytes(b'%PDF-partial'); raise RuntimeError('fixture renderer failure')
        with tempfile.TemporaryDirectory() as directory:
            receipt = audit.write_chain_artifacts(directory,result,fail)
            self.assertEqual(receipt['status'],'unavailable')
            self.assertEqual([p.suffix for p in Path(directory).iterdir()], ['.json'])
            self.assertIn('fixture renderer failure',json.loads(Path(receipt['json']).read_bytes())['pdf']['reason'])

    def test_altered_receipt_identity_is_not_reused(self):
        result = audit.build_chain_audit(*fixture())
        for field, value in [('schema_version','wrong'),('status','accepted')]:
            with tempfile.TemporaryDirectory() as directory:
                receipt = audit.write_chain_artifacts(directory,result,fake_pdf)
                path = Path(receipt['json']); changed=json.loads(path.read_bytes()); changed[field]=value
                path.write_bytes(audit.encoded(changed))
                with self.assertRaises(FileExistsError): audit.write_chain_artifacts(directory,result,fake_pdf)

    def test_event_weight_axis_mismatch_is_unavailable(self):
        data=fixture()
        data[2][0]['final_diagnostic_application_result']['H_pion_control_model'].edges[1]=1.02
        p=audit.build_chain_audit(*data)['parents'][0]
        self.assertEqual(p['final']['nommcuts_event_reproduction']['status'],'unavailable')

    def test_prior_pdf_collision_and_nonexistent_output_are_not_overwritten(self):
        result = audit.build_chain_audit(*fixture())
        with tempfile.TemporaryDirectory() as directory:
            stem = 'pion_background_chain_v1_Q4p4W2p74_Center_lowe_'+result['collection_identity']
            path = Path(directory,stem+'.pdf'); path.write_bytes(b'prior')
            with self.assertRaises(FileExistsError): audit.write_chain_artifacts(directory,result,fake_pdf)
            self.assertEqual(path.read_bytes(),b'prior')
            with self.assertRaises(FileNotFoundError): audit.write_chain_artifacts(Path(directory,'absent'),result,fake_pdf)

    def test_nonfinite_histogram_and_weight_cannot_claim_reproduction(self):
        for mode in ('histogram','weight'):
            data = fixture()
            if mode == 'histogram':
                data[2][0]['H_pion_control'].values[1] = float('nan')
            else:
                data[2][0]['proposed_diagnostic_application_result']['weights'][1] = float('inf')
            result = audit.build_chain_audit(*data)
            json.loads(audit.encoded(result))
            closure = result['parents'][0]['proposal']['nommcuts_event_reproduction']
            self.assertEqual(closure['status'],'unavailable')

    def test_zero_weights_are_observed_without_model_replacement(self):
        data = fixture(); payload = data[2][0]['proposed_diagnostic_application_result']
        payload['weights'] = [0.] * 4
        payload['H_pion_subtraction_template_MM_nosub'] = Hist()
        p = audit.build_chain_audit(*data)['parents'][0]
        self.assertEqual(p['proposal']['nommcuts_event_reproduction']['content']['status'],'PASS')
        self.assertEqual(p['proposal']['nommcuts_weighted_sources']['combined']['absolute'],[0.]*4)

    def test_strong_cancellation_retains_absolute_support(self):
        data = fixture(); parent = data[2][0]
        records = data[0]['_authoritative_pion_control_source_cache']['by_t'][0]['records']
        records[0]['coefficient'],records[1]['coefficient'] = 1e9,-1e9
        source = audit.source_spectra(records,audit.histogram(parent['H_pion_control']),parent)
        self.assertEqual(source['combined']['signed'][1],-1.)
        self.assertEqual(source['combined']['absolute'][1],2e9+7.)

    def test_zero_component_model_and_missing_weights_are_explicit(self):
        data = fixture(); fit = data[2][0]['fit_result']
        fit['diagnostics']['pion']['refined_amplitudes']['pi_n'] = 1e-30
        self.assertEqual(audit.build_chain_audit(*data)['parents'][0]['transfer']['selected_amplitudes']['pion']['pi_n'],1e-30)
        data[2][0]['final_diagnostic_application_result'].pop('weights')
        self.assertEqual(audit.build_chain_audit(*data)['parents'][0]['final']['nommcuts_event_reproduction']['status'],'unavailable')

    def test_observation_failure_does_not_raise_into_production(self):
        with patch.object(audit,'build_chain_audit',side_effect=ValueError('capture')), redirect_stdout(io.StringIO()) as output:
            result = audit.observe_chain_audit(*fixture(),'unused')
        self.assertEqual(result['status'],'unavailable'); self.assertIn('capture',output.getvalue())

    def test_parent_hooks_new_and_reuse_preserve_existing_anchor_and_pages(self):
        source = (REPO/'src/cuts/pion_t_bin_parents.py').read_text()
        tree = ast.parse(source)
        function = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='build_setting_t_bin_pion_parents')
        calls = [n for n in ast.walk(function) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name)]
        self.assertEqual(sum(n.func.id=='observe_chain_audit' for n in calls),2)
        self.assertEqual(sum(n.func.id=='observe_anchor_audit' for n in calls),2)
        # Execute the real reuse branch with scalar boundary doubles; no fit or renderer dependencies.
        data = fixture(); parents = data[2]; data[0]['_pion_t_bin_parent_results'] = parents
        contracts = [dict(identity=dict(pion_parent_id=p['pion_parent_id']),
                          proton_producer_fingerprint=None, pion_control_input_fingerprint=p['pion_control_input_fingerprint'],
                          coordinate_fingerprint='coord') for p in parents]
        env = dict(_is_component_t_bin_production=lambda inp: True,
            validate_kaon_data_coordinate_contract=lambda c,**kw:c,
            _parent_input_contract=Mock(side_effect=contracts),
            validate_authoritative_t_bin_pion_parent=Mock(), observe_anchor_audit=Mock(),
            observe_chain_audit=Mock(), _finalize_frozen_parent_collection=lambda h,i,p,t,**kw:p)
        exec(compile(ast.Module(body=[function],type_ignores=[]),'parent_builder','exec'),env)
        result = env['build_setting_t_bin_pion_parents'](data[0],data[1],t_bins=[.4,.6,.8,1.],
            parent_inputs=[{}]*3,pion_component_shape_payload={},kaon_signal_shape_payload={},
            kaon_sigma0_shape_payload={},parent_pion_alignment={},alignment_outpath='unused',mm_offset_data=0.,coordinate_contract={})
        self.assertIs(result,parents)
        env['observe_chain_audit'].assert_called_once_with(data[0],data[1],parents,{},'unused')
        env['observe_anchor_audit'].assert_called_once()
        self.assertNotIn('page_manifest',Path(audit.__file__).read_text())


class FakeRoot:
    def __init__(self): self.calls=[]; self.graphs=[]
    def TCanvas(self,*args): return self.Object(self,'canvas')
    def TGraph(self,n,x,y):
        self.graphs.append(list(y)); return self.Object(self,'graph')
    def TLegend(self,*args): return self.Object(self,'legend')
    def TLatex(self,*args): return self.Object(self,'latex')
    def TPaveText(self,*args): return self.Object(self,'text')
    def TLine(self,*args): return self.Object(self,'marker')
    class Object:
        def __init__(self,root,kind): self.root,self.kind=root,kind
        def __getattr__(self,name):
            def call(*args): self.root.calls.append((self.kind,name,args))
            return call


class RenderingTests(unittest.TestCase):
    def test_fake_root_panels_titles_axes_legends_count_and_immutability(self):
        result = audit.build_chain_audit(*fixture()); original = copy.deepcopy(result)
        root = FakeRoot(); plan = audit.page_specs(result)
        count = audit.render_companion(Path('fake.pdf'),result,root)
        self.assertEqual(count,len(plan)); self.assertEqual(result,original)
        titles = [args[0] for kind,name,args in root.calls if name=='SetTitle']
        self.assertTrue(any(';MM [GeV];weight' in title for title in titles))
        self.assertTrue(any(kind=='legend' and name=='AddEntry' for kind,name,args in root.calls))
        self.assertTrue(any(kind=='marker' and name=='Draw' for kind,name,args in root.calls))
        self.assertEqual(sum(name=='Print' for kind,name,args in root.calls),len(plan)+2)
        self.assertEqual(sum('Transfer' in p['title'] for p in plan),3)
        self.assertTrue(all(any(stage in p['title'] for p in plan) for stage in ('Control construction','SIMC','Application and outcome')))

    def test_mismatched_axes_are_not_plotted_on_control_axis(self):
        data=fixture(); data[2][0]['fit_result']['H_simc_shape_pi_n'].edges[1]=1.02
        plan=audit.page_specs(audit.build_chain_audit(*data))
        aligned=next(p for p in plan if 't1 - SIMC' in p['title'])['panels'][0]
        self.assertNotIn('pi_n',[name for name,values in aligned['series']])


if __name__ == '__main__':
    unittest.main()
