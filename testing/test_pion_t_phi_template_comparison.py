"""Deterministic tests. Fake histograms are not ROOT/farm validation."""
from copy import deepcopy
import ast
import importlib.util
import json
import math
from pathlib import Path
import sys
import tempfile
import types
import unittest
from contextlib import contextmanager, ExitStack
from unittest.mock import patch

REPO=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(REPO/'src/cuts'),str(REPO/'src/utility')]
import pion_t_phi_template_comparison as child
from testing.test_pion_t_template_comparison import Hist
from testing.test_pion_component_dynamic_alignment import fits


class H(Hist):
    def InheritsFrom(self, name): return name == 'TH1'
    def Add(self,other,factor=1.):
        self.values=[a+factor*b for a,b in zip(self.values,other.values)]
        self.errors=[math.sqrt(a*a+factor*factor*b*b) for a,b in zip(self.errors,other.errors)]
    def Scale(self,factor):
        self.values=[v*factor for v in self.values]; self.errors=[e*abs(factor) for e in self.errors]
    def Reset(self,*a):
        self.values=[0.]*len(self.values); self.errors=[0.]*len(self.errors)
    def Integral(self,first=None,last=None): return sum(self.values[1:-1] if first is None else self.values[first:last+1])
    def FindBin(self,x):
        from bisect import bisect_right
        return bisect_right(self.edges,x)
    def GetName(self): return 'test_hist'


def row(**kw):
    result=dict(source='prompt',entry=0,mm=.5,ti=0,pi=0,weight=1.,selected=True,allcuts=True,nommcuts=True)
    result.update(kw)
    return result


class ProjectionTests(unittest.TestCase):
    def test_signed_cancellation_retains_quadratic_sumw2(self):
        r=child.project([row(weight=2.),row(source='rand',weight=-2.)],[0.,1.,2.],1,2)
        self.assertEqual(r['parents'][0]['values'],[0.,0.,0.,0.])
        self.assertEqual(r['parents'][0]['variances'],[0.,8.,0.,0.])
        self.assertEqual(r['closure'],'PASS')

    def test_all_sources_flows_and_outside_phi_close(self):
        rows=[row(source=s,entry=i,weight=w,mm=mm,pi=pi) for i,(s,w,mm,pi) in enumerate([
            ('prompt',3.,-1.,0),('rand',-1.,0.,1),('dummy',-2.,1.,None),('dummy_rand',.5,2.,-1)])]
        r=child.project(rows,[0.,1.,2.],1,2)
        self.assertEqual(r['parents'][0]['values'],[3.,-1.,-2.,.5])
        self.assertEqual(r['parents'][0]['variances'],[9.,1.,4.,.25])
        self.assertEqual(r['outside_phi'][0]['count'],2)

    def test_invalid_event_failures(self):
        for bad,reason in ((row(mm=float('nan')),'nonfinite'),(row(weight=float('inf')),'nonfinite'),
                           (row(ti=-1),'parent'),(row(pi=9),'phi'),(row(source='other'),'source')):
            with self.subTest(bad=bad),self.assertRaisesRegex(ValueError,reason):
                child.project([bad],[0.,1.,2.],1,2)
        with self.assertRaisesRegex(ValueError,'duplicate'):
            child.project([row(),row()],[0.,1.,2.],1,2)

    def test_mask_exclusion_and_empty_cells(self):
        r=child.project([row(selected=False)],[0.,1.,2.],1,2)
        self.assertEqual(r['parents'][0]['count'],0)
        self.assertEqual(r['children'][0][1]['variances'],[0.]*4)

    def test_hist_reconstruction_keeps_negative_bins_and_flows(self):
        cell=dict(values=[2.,-3.,4.,5.],variances=[4.,9.,16.,25.],count=4)
        original=H()
        h=child.make_hist(original,cell,'new')
        self.assertEqual(h.values,cell['values']); self.assertEqual(h.errors,[2.,3.,4.,5.])
        self.assertEqual(original.values,[0.,2.,6.,0.])


class GuardTests(unittest.TestCase):
    def test_disabled_import_requires_no_ROOT(self):
        self.assertIsNone(child.observe({}, {}, [], [], '.', None, None))
        import subprocess
        subprocess.run([sys.executable,'-B','-c',
            "import sys;sys.path.insert(0,'src/cuts');import pion_t_phi_template_comparison;assert 'ROOT' not in sys.modules;assert 'numpy' not in sys.modules"],cwd=REPO,check=True)

    def test_setting_prepass_and_explicit_enable(self):
        inp=dict(Q2='4p4',W='2p74',EPSSET='low',ParticleType='kaon')
        env={child.ENABLE_ENV:'1'}
        self.assertTrue(child.requested(inp,'Center',env))
        for setting in ('Right','Other'):
            self.assertFalse(child.requested(inp,setting,env))
        for key,val in (('bg_opt_prepass_mode',True),('analysis_subtraction_context','bg_optimizer_pre_parent'),('EPSSET','high'),('ParticleType','pion')):
            self.assertFalse(child.requested(dict(inp,**{key:val}),'Center',env))
        for val in ('0','true',''):
            self.assertFalse(child.requested(inp,'Center',{child.ENABLE_ENV:val}))

    def test_integration_guard_and_unconditional_default_preserved(self):
        tree=ast.parse((REPO/'src/binning/calculate_yield.py').read_text())
        fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='calculate_yield_data')
        hook=next(n for n in fn.body if isinstance(n,ast.If) and 'observe_t_phi_comparison' in ast.unparse(n))
        text=ast.unparse(hook.test)
        for token in (child.ENABLE_ENV,'BG_OPT_PREPASS_CONTEXT','bg_opt_prepass_mode','simc_shape_components','t_bin','kaon'):
            self.assertIn(token,text)
        cache=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and '_yield_sub_event_cache' in ast.unparse(n))
        self.assertGreater(fn.body.index(hook),cache)


class ScopeAndFitterTests(unittest.TestCase):
    def test_actual_resolver_direct_t_setting_and_missing(self):
        direct=H(); fallback=H(values=(0.,0.,0.,0.)); setting=H()
        for local,global_,expected in ((direct,setting,'local_t_phi'),(fallback,setting,'setting_fallback'),(None,None,'unavailable')):
            payload=dict(binned_shapes={'t_bin1':{'phi_bin1':local}},setting_shape_full=global_)
            with patch.dict(sys.modules,{'pion_component_fits':fits}),patch.object(fits,'_sum_hist_list_to_unit_area',return_value=None):
                selected=fits._resolve_single_scope_hist(payload,'test',0,0)
                self.assertEqual(child.scope_evidence(payload,0,0,selected,fits._hist_has_usable_support),expected)
        payload=dict(binned_shapes={'t_bin1':{'phi_bin2':direct}},setting_shape_full=setting)
        with patch.dict(sys.modules,{'pion_component_fits':fits}),patch.object(fits,'_sum_hist_list_to_unit_area',return_value=direct):
            self.assertEqual(child.scope_evidence(payload,0,0,direct,fits._hist_has_usable_support),'local_t_aggregate_fallback')

    def test_cell_clones_all_inputs_and_passes_existing_procedure(self):
        pion,kaon=H(),H(); shapes={c:H() for c in child.COMPONENTS}
        parent=dict(t_edges=[.4,.6],fit_result=dict(template_mm_offset_data=.01,particle_subtraction_phi_setting='Center',pion_component_alignment={'parent':'map'}))
        inp={'physics':{'frozen':1}}
        def fitter(p,k,s,config,**kw):
            self.assertIsNot(p,pion); self.assertIsNot(k,kaon)
            self.assertEqual(p.values,pion.values); self.assertEqual(k.errors,kaon.errors)
            self.assertEqual(config,inp); self.assertEqual(kw['analysis_scope'],'t_bin1phi_bin2')
            self.assertEqual(kw['parent_alignment'],{'parent':'map'})
            self.assertNotIn('pion_component_alignment',kw)  # existing child resolver runs
            config['physics']['frozen']=99; s['pi_n'].values[1]=999
            return {'fresh_fit':True}
        r=child.fit_cell(pion,kaon,shapes,inp,ti=0,pi=1,parent=parent,signal=H(),sigma0=None,sigma0_diagnostics={},fitter=fitter)
        self.assertTrue(r['fresh_fit']); self.assertEqual(inp['physics']['frozen'],1)
        self.assertEqual(shapes['pi_n'].values[1],2.)

    def test_actual_public_builder_scope_equivalence_with_deterministic_inner_fits(self):
        # Executes the real public builder; only ROOT shift/clone plumbing and
        # the two inner solvers are fixtures. No claim of real ROOT numerics.
        from contextlib import ExitStack
        parent=dict(t_edges=[.4,.6],fit_result=dict(template_mm_offset_data=0.,particle_subtraction_phi_setting='Center'))
        pion,kaon=H(),H(); shapes={c:H() for c in child.COMPONENTS}
        def inner(side):
            h=H()
            return dict(**{side+'_n':1.,side+'_delta':2.,side+'_sidis':3.},
                B_sigma0=None,S_lambda=None,S_sigma0=None,fit_status=0,diagnostics={},
                pi_n_scaled_hist=h,pi_delta_scaled_hist=h,pi_sidis_scaled_hist=h,
                k_sigma0_scaled_hist=None,k_lambda_scaled_hist=None,fit_hist=h,pion_bg_fit_hist=h,residual_hist=h)
        with ExitStack() as stack:
            for name,value in {
                'configure_particle_subtraction_root_ownership_debug':lambda *a:None,
                'resolve_particle_subtraction_mode':lambda *a:'simc_shape_components',
                'get_pion_component_dynamic_alignment_config':lambda **kw:{'enabled':True},
                'resolve_pion_component_alignment':lambda *a,**kw:{'source':'fixture-alignment','components':{c:{'total_shift_gev':0.,'effective_windows':[[0.,2.]]} for c in child.COMPONENTS}},
                '_build_alignment_template':lambda h,*a,**kw:(deepcopy(h),{}),
                '_build_mm_shifted_hist':lambda h,*a,**kw:deepcopy(h),
                '_clone_hist':lambda h,*a,**kw:deepcopy(h),
                '_promote_sparse_zero_solution_if_applicable':lambda v,*a:v,
                'fit_pion_control_with_simc_shapes':lambda *a,**kw:inner('B'),
                'fit_kaon_nosub_with_simc_pion_shapes':lambda *a,**kw:inner('A'),
            }.items(): stack.enter_context(patch.object(fits,name,value))
            a=fits.build_particle_subtraction_component_result(pion,kaon,shapes,{},analysis_scope='t_bin1',kaon_signal_shape=H(),phi_setting='Center')
            b=child.fit_cell(pion,kaon,shapes,{},ti=0,pi=0,parent=parent,signal=H(),sigma0=None,sigma0_diagnostics={},fitter=fits.build_particle_subtraction_component_result)
        for key in ('A_n','A_delta','A_sidis','B_n','B_delta','B_sidis','fit_status_pion','fit_status_kaon','resolved_subtraction_config','pion_component_alignment'):
            self.assertEqual(a[key],b[key],key)
        self.assertEqual(a['H_pion_control_input'].values,b['H_pion_control_input'].values)


class CaptureRowsTests(unittest.TestCase):
    def hist(self):
        p=dict(source_label='prompt',coefficient=2.,coordinate_fingerprint='x',entry_index=4,adj_MM=.5,adj_t=.5,phi_degrees=0.,t_index=0,phi_index=0,allcuts=False,nommcuts=True)
        k=dict(adj_MM=[.5,.6],source_label=['dummy_prompt','rand'],source_entry_index=[4,5],t_index=[0,0],phi_index=[0,-1],allcuts=[True,False],nommcuts=[True,True],proton_cleaning_factor=[.5,.25],proton_rf_accept=[True,False])
        return dict(normfac_data=2.,normfac_dummy=3.,nWindows=2.,_authoritative_pion_control_source_cache={'records':[p],'coordinate_fingerprint':'x'},_yield_data_event_cache={'s':k})

    def test_source_normalization_proton_RF_and_cut_masks(self):
        rr,coeff=child.rows_from_caches(self.hist())
        self.assertEqual(coeff,dict(prompt=2.,rand=-1.,dummy=-3.,dummy_rand=1.5))
        self.assertEqual(rr['kaon'][0]['weight'],-1.5); self.assertEqual(rr['kaon'][1]['weight'],0.)
        self.assertFalse(rr['pion'][0]['allcuts']); self.assertTrue(rr['pion'][0]['nommcuts'])
        self.assertEqual(rr['pion'][0]['weight'],2.)

    def test_invalid_cache_and_provenance(self):
        for mutation,reason in ((lambda h:h.update(nWindows=0),'normalization'),
            (lambda h:h['_yield_data_event_cache']['s']['adj_MM'].pop(),'column'),
            (lambda h:h['_authoritative_pion_control_source_cache']['records'][0].update(coefficient=7.),'provenance')):
            h=self.hist(); mutation(h)
            with self.assertRaisesRegex(ValueError,reason): child.rows_from_caches(h)


class TransferTests(unittest.TestCase):
    def test_application_aggregation_quadratic_errors_and_missing_status(self):
        app=lambda v:dict(status='available',cuts={'nommcuts':{'values':[v,-v], 'variances':[v*v,v*v]}})
        cells=[{'applications':{'X':app(2.)}},{'applications':{'X':app(-1.)}}]
        aggregate=child.aggregate_applications(cells,'X','nommcuts')
        self.assertEqual(aggregate['values'],[1.,-1.]); self.assertEqual(aggregate['variances'],[5.,5.])
        self.assertEqual(child.aggregate_applications(cells,'missing','nommcuts')['status'],'unavailable')

    def test_actual_weight_application_and_paired_covariance(self):
        config=types.ModuleType('background_config')
        config.resolve_particle_subtraction_weight_denominator_floor=lambda inp:1e-12
        config.resolve_particle_subtraction_weight_clip_bounds=lambda inp:(0.,None)
        sub=types.ModuleType('pion_component_subtraction')
        sub.build_simc_shape_pion_control_weights=lambda fit,**kw:dict(weights=fit['w'],diagnostics=kw)
        cells={cut:dict(values=[0.,-2.,3.,0.],variances=[0.,4.,9.,0.]) for cut in ('allcuts','nommcuts')}
        inp=dict(mm_min=.25,mm_max=1.5)
        with patch.dict(sys.modules,{'background_config':config,'pion_component_subtraction':sub}):
            a=child.apply_transfer({'w':[0.,2.,3.,0.]},cells,inp)
            b=child.apply_transfer({'w':[0.,1.,2.,0.]},cells,inp)
            self.assertEqual(a['cuts']['allcuts']['values'],[0.,-4.,9.,0.])
            self.assertEqual(a['cuts']['nommcuts']['variances'],[0.,16.,81.,0.])
            paired=child.paired_fixed_weight({'A':a,'B':b},cells)['A__B']['nommcuts']
            self.assertEqual(paired['covariance'],[0.,8.,54.,0.])
            self.assertEqual(paired['difference_variance'],[0.,4.,9.,0.])
            for i,v in enumerate(paired['difference_variance']):
                self.assertAlmostEqual(v,a['cuts']['nommcuts']['variances'][i]+b['cuts']['nommcuts']['variances'][i]-2*paired['covariance'][i])
            with self.assertRaisesRegex(ValueError,'invalid_transfer'):
                child.apply_transfer({'w':[1.,2.,3.,0.]},cells,inp)


class ObserverTests(unittest.TestCase):
    @contextmanager
    def environment(self, mismatch=False, fit_error=False, empty=False, setting='Center',one_empty=False,
                    rejected=False, fallback=False):
        import pion_t_template_comparison as comparison
        import root_histogram_ownership as ownership
        fake_root=types.ModuleType('ROOT'); fake_root.gDirectory=None; fake_root.gPad=None
        cuts=types.ModuleType('apply_cuts'); cuts.science='original'
        shape_module=types.ModuleType('pion_component_shapes')
        hole_module=types.ModuleType('hgcer_hole'); hole_module.apply_HGCer_hole_cut=lambda *a:object()
        rows={'pion':[],'kaon':[]}
        for ti in range(3):
            for pi in range(9):
                for pop in rows:
                    rows[pop].append(row(ti=ti,pi=pi,entry=9*ti+pi,weight=0. if empty or (one_empty and pop=='kaon') else 1.))
        parents=[]
        proj=child.project(rows['pion'],[0.,1.,2.],3,9)
        for ti in range(3):
            h=child.make_hist(H(),proj['parents'][ti],'parent')
            kh=child.make_hist(H(),child.project(rows['kaon'],[0.,1.,2.],3,9)['parents'][ti],'kaon_parent')
            parents.append(dict(t_bin_index=ti,t_edges=[float(ti),float(ti+1)],analysis_scope=f't_bin{ti+1}',
                H_pion_control=h,H_proton_cleaned_final_rf=kh,fit_result=dict(template_mm_offset_data=0.,particle_subtraction_phi_setting=setting)))
        if mismatch: parents[0]['H_pion_control'].values[1]+=1
        shape_payload=dict(setting_shape_full=H(),binned_shapes={f't_bin{ti+1}':{f'phi_bin{pi+1}':H() for pi in range(9)} for ti in range(3)})
        components=dict(components={c:deepcopy(shape_payload) for c in child.COMPONENTS})
        if fallback:
            components['components']['pi_n']['binned_shapes'] = {}
        called=[]
        def loader(*a,**kw):
            cuts.science='temporary loader state'
            return components
        shape_module.load_setting_pion_component_shapes=loader
        shape_module.load_raw_t_pion_comparison_shapes=lambda *a,**kw:{c:{'parents':[{}]*3} for c in child.COMPONENTS}
        def fitter(p,k,s,inp,**kw):
            called.append((p,k,s,kw))
            self.assertEqual(cuts.science,'original')
            ownership._DEBUG_ENABLED=True; ownership._DEBUG_RECORDS[:]=['temporary']
            if fit_error: raise ValueError('deliberate child fitter failure')
            return dict(H_pion_control_input=p,H_kaon_nosub_input=k,analysis_scope=kw['analysis_scope'],
                particle_subtraction_mode='simc_shape_components',fit_status_pion='success',
                fit_status_kaon='failed' if rejected else 'success',fallback_used=False,
                diagnostics={pop:dict(validation=dict(accepted=not (rejected and pop=='kaon')))
                             for pop in ('pion','kaon')},
                **{key:1. for key in ('A_n','A_delta','A_sidis','B_n','B_delta','B_sidis')},
                **{'H_simc_shape_'+c:h for c,h in s.items()})
        inp=dict(Q2='4p4',W='2p74',EPSSET='low',ParticleType='kaon',mm_min=1.1,mm_max=1.13)
        with tempfile.TemporaryDirectory() as temp,ExitStack() as stack:
            inputfile=Path(temp)/'runtime.root'; inputfile.write_bytes(b'fixture runtime only')
            backing=types.SimpleNamespace(GetName=lambda:str(inputfile))
            hist=dict(phi_setting=setting,_pion_t_bin_parent_results=parents,
                      _authoritative_pion_control_source_cache={'kaon_data_coordinate':{'fixture':True}},
                      InFile_DATA=backing,InFile_DUMMY=backing)
            stack.enter_context(patch.dict(os_environ(),{child.ENABLE_ENV:'1'}))
            stack.enter_context(patch.dict(sys.modules,{'ROOT':fake_root,'apply_cuts':cuts,'pion_component_shapes':shape_module,
                'pion_component_fits':fits,'hgcer_hole':hole_module}))
            stack.enter_context(patch.object(child,'rows_from_caches',return_value=(rows,{})))
            stack.enter_context(patch.object(fits,'build_particle_subtraction_component_result',fitter))
            stack.enter_context(patch.object(fits,'normalize_hist_to_unit_area',lambda h,**kw:h.Scale(1./h.Integral()) if h.Integral()>0 else None))
            stack.enter_context(patch.object(comparison,'compare_parent',return_value=({'fixture':True},{})))
            stack.enter_context(patch.object(child,'apply_transfer',return_value={'status':'unavailable','reason':'fixture'}))
            prior=ownership._DEBUG_ENABLED,list(ownership._DEBUG_RECORDS)
            try:
                yield lambda:child.observe(hist,inp,[0.,1.,2.,3.],list(map(float,range(10))),temp,shape_payload,None),called,cuts,Path(temp),prior
            finally:
                self.assertEqual((ownership._DEBUG_ENABLED,list(ownership._DEBUG_RECORDS)),prior)
                self.assertEqual(cuts.science,'original')

    def test_all_27_independent_child_calls_capture_and_no_parent_mutation(self):
        with self.environment() as (run,called,cuts,temp,prior):
            result=run()
            self.assertEqual(len(called),27)
            self.assertEqual(len({k['analysis_scope'] for p,q,s,k in called}),27)
            self.assertTrue(all(p.values[1]==1. and q.values[1]==1. for p,q,s,k in called))
            self.assertTrue(result['production_parents_unchanged'])
            self.assertEqual(result['parent_input_closure'],'PASS')
            self.assertEqual(result['cell_status_counts'],{'fitted_local_templates':27})
            saved=json.loads(next(temp.glob('*.json')).read_text())
            self.assertEqual(saved['cells'][26]['fit']['analysis_scope'],'t_bin3phi_bin9')
            with self.assertRaisesRegex(ValueError,'destination_exists'): run()

    def test_parent_mismatch_stops_before_any_fit_and_records_failure(self):
        with self.environment(mismatch=True) as (run,called,cuts,temp,prior):
            with self.assertRaisesRegex(ValueError,'parent_observation_closure'): run()
            self.assertEqual(called,[])
            saved=json.loads(next(temp.glob('*.json')).read_text())
            self.assertEqual(saved['status'],'failed')

    def test_fit_failures_explicit_and_state_restored(self):
        with self.environment(fit_error=True) as (run,called,cuts,temp,prior):
            result=run()
            self.assertEqual(len(called),27)
            self.assertEqual(result['cell_status_counts'],{'fit_failed':27})
            self.assertTrue(all('deliberate' in c['reason'] for c in result['cells']))
            self.assertTrue(all(c['fit'] is None for c in result['cells']))

    def test_empty_cells_never_relabel_parent_coefficients(self):
        with self.environment(empty=True) as (run,called,cuts,temp,prior):
            result=run()
            self.assertEqual(called,[])
            self.assertTrue(all(c['reason']=='empty_observed_cell' for c in result['cells']))
            self.assertEqual(result['cell_status_counts'],{'empty_child':27})
            self.assertTrue(all(not c['fit_returned'] and not c['fit_accepted'] for c in result['cells']))

    def test_returned_rejected_gate_preserves_fit_and_excludes_child_weights(self):
        with self.environment(rejected=True) as (run,called,cuts,temp,prior):
            with patch.object(child,'apply_transfer',return_value={'status':'unavailable','reason':'A/B fixture'}) as transfer:
                result=run()
            self.assertEqual(len(called),27)
            self.assertEqual(transfer.call_count,54)  # only the separate A/B branches
            self.assertEqual(result['cell_status_counts'],{'fit_rejected':27})
            for cell in result['cells']:
                self.assertTrue(cell['fit_returned']); self.assertFalse(cell['fit_accepted'])
                self.assertIsNotNone(cell['fit'])
                self.assertIn('kaon validation rejected',cell['fit_acceptance']['reason'])
                self.assertFalse(cell['fit_acceptance']['diagnostics']['fit_validation_kaon'])
                self.assertEqual(cell['applications']['Local_t_phi']['status'],'unavailable')

    def test_accepted_fallback_and_genuine_local_statuses(self):
        for fallback,status in ((False,'fitted_local_templates'),(True,'fitted_with_template_fallback')):
            with self.subTest(fallback=fallback), self.environment(fallback=fallback) as (run,*_):
                result=run()
                self.assertEqual(result['cell_status_counts'],{status:27})
                self.assertTrue(all(c['fit_returned'] and c['fit_accepted'] for c in result['cells']))
                self.assertTrue(all(c['template_scopes']['pi_n']==('setting_fallback' if fallback else 'local_t_phi')
                                    for c in result['cells']))

    def test_owner_verifies_actual_observer_fixture_for_both_settings(self):
        from testing import run_pion_t_template_comparison_gate as owner
        for setting in ('Center','Left'):
            for rejected,empty in ((False,False),(True,False),(False,True)):
                with self.subTest(setting=setting,rejected=rejected,empty=empty),self.environment(
                        setting=setting,rejected=rejected,empty=empty) as (run,*_):
                    doc=run()
                    ab=dict(canonical_t_edges=doc['t_edges'],canonical_phi_edges=doc['phi_edges'],
                        original_analysis_input_identity=deepcopy(doc['runtime_inputs']),parents=[])
                    for ti in range(3):
                        frozen={}
                        for pop in ('pion','kaon'):
                            p=doc['observations'][pop]['nommcuts']
                            frozen[pop]=dict(edges=p['edges'],values=p['parents'][ti]['values'],
                                stored_sumw2=p['parents'][ti]['variances'])
                        ab['parents'].append(dict(frozen_inputs=frozen))
                    identities={p:{k:v[k] for k in ('bytes','sha256')} for p,v in doc['source_files'].items()}
                    outcome=owner.verify_child_payload(doc,setting,doc['source_HEAD'],identities,ab)
                    self.assertEqual(outcome['accepted_fits'],0 if rejected or empty else 27)

    def test_renderer_uses_correct_rejected_and_fallback_status_titles(self):
        from testing.render_pion_t_phi_template_comparison import render
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_pdf import PdfPages
        original=Figure.suptitle
        for kwargs,status in ((dict(rejected=True),'fit_rejected'),(dict(fallback=True),'fitted_with_template_fallback')):
            with self.subTest(status=status),self.environment(**kwargs) as (run,called,cuts,temp,prior):
                doc=run(); doc['cells']=doc['cells'][:1]
                with patch.object(Figure,'suptitle',autospec=True,side_effect=original) as title,patch.object(PdfPages,'savefig'):
                    render(doc,temp/'renderer-status.pdf')
                self.assertIn(status,title.call_args_list[0].args[1])

    def test_left_separate_and_one_empty_population_is_owned_by_existing_fitter(self):
        with self.environment(setting='Left',one_empty=True) as (run,called,cuts,temp,prior):
            result=run()
            self.assertEqual(result['phi_setting'],'Left')
            self.assertEqual(len(called),27)
            self.assertTrue(all(not any(k.values) for p,k,s,kw in called))


def os_environ():
    import os
    return os.environ


class RendererTests(unittest.TestCase):
    def test_window_sum_matches_actual_utility_source(self):
        from testing.render_pion_t_phi_template_comparison import window_sum
        node=next(n for n in ast.parse((REPO/'src/utility/utility.py').read_text()).body
                  if isinstance(n,ast.FunctionDef) and n.name=='integrate_hist_range')
        ns={'ROOT':types.SimpleNamespace(TH1F=H)}
        exec(compile(ast.Module(body=[node],type_ignores=[]),'<actual-integration>','exec'),ns)
        h=H()
        for low,high in ((-.1,3.),(.5,1.),(0.,2.),(2.,3.)):
            self.assertEqual(window_sum(h.values,h.edges,low,high),ns['integrate_hist_range'](h,low,high))

    def test_failed_or_wrong_schema_is_never_rendered(self):
        from testing.render_pion_t_phi_template_comparison import render
        with tempfile.TemporaryDirectory() as temp:
            for doc in ({},{'schema_version':child.SCHEMA,'status':'failed','production_parents_unchanged':True}):
                with self.assertRaises(ValueError): render(doc,Path(temp)/'uncreated.pdf')
            self.assertEqual(list(Path(temp).iterdir()),[])


if __name__=='__main__': unittest.main()
