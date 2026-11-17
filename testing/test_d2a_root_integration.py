"""ROOT-free owner/checker regressions; record fixtures are NOT ROOT execution."""
from copy import deepcopy
import csv
import ast
import os
import importlib.util
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch
import zlib

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("d2a_integration_owner", REPO/'testing/run_d2a_root_integration.py')
m = importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
SEED = m.load_json(m.SEED)


def histogram(raw, edges):
    scale = raw['normal_sumw'] if raw['normal_sumw'] > 0 else 1.
    variances = [v/(scale*scale) for v in raw['variances']]
    return dict(status='available', edges=deepcopy(edges), values=[v/scale for v in raw['values']], errors=[v**.5 for v in variances], stored_sumw2=variances)


def record_fixture(seed, plan):
    """Labelled record-shaped positive fixture. No fit or ROOT source is simulated."""
    payloads = {}
    bindings = {}
    cache = {}
    for name in m.ROLES:
        if plan['inputs'][name] is None:
            payloads[name] = dict(status='not_configured',scope_usable=False)
            continue
        support = deepcopy(seed['expected_pion_support'].get(name,seed['expected_pion_support']['pi_n']))
        identity = dict(schema=2,population_definition=3,reuse_safe=True,component=name,particle_type='kaon',phi_setting='Center',t_edges=seed['t_edges'],phi_edges=seed['phi_edges'],hgcer_selection={'points':seed['hgcer_points']},mm_axis=dict(nbins=100,min=.7,max=1.3,flows=True),mm_selection_window=[seed['config']['mm_min'],seed['config']['mm_max']],normalization='normal_bins_once_after_event_accumulation',weight='getattr(evt, iter_weight, 1.0)')
        item=plan['inputs'][name]
        token=m.file_token(Path(item['path']).resolve()) if Path(item['path']).is_file() else [item['path'],1,item['bytes'],1,1,1]
        canonical=str(Path(item['path']).resolve()) if Path(item['path']).is_file() else item['path']
        bindings[name]=dict(requested_path=item['path'],resolved_path=token[0],canonical_path=canonical,token=deepcopy(token),resolution={})
        identity.update(source=deepcopy(token),tree='h10',selection=m.selection_identity(seed['config']))
        payloads[name] = dict(status='available',scope_usable=True,identity=identity,support=support,setting=histogram(support['setting'],seed['mm_edges']),
            parents={k:histogram(v,seed['mm_edges']) for k,v in support['parents'].items()},children={k:{p:histogram(v,seed['mm_edges']) for p,v in children.items()} for k,children in support['children'].items()})
        cache[name] = dict(before_entries=0,cold_entries=1,warm_entries=1,probe_entries=2,warm_equal=True,clone_independent=True,probe_parent_equal=True,probe_phi_edges=[-180.,180.])
    fits = {}
    for mode in ['disabled','enabled']:
        snap={k:0. for k in ['A_n','A_delta','A_sidis','B_n','B_delta','B_sidis']}
        snap.update(H_kaon_nosub_input=deepcopy(seed['experimental']['kaon']),H_pion_control_input=deepcopy(seed['experimental']['pion']),pion_component_alignment=deepcopy(seed['alignment']),staged_amplitudes_raw={k:dict(pi_n=2.,pi_delta=2.,pi_sidis=2.) for k in ['kaon_nosub','pion_control']},resolved_subtraction_config=deepcopy(seed['resolved_subtraction_config']),fallback_used=False,fit_status_pion='success',fit_status_kaon='success',diagnostics={})
        for side,key in [('pion','H_pion_fit_step_overlays'),('kaon','H_kaon_fit_step_overlays')]:
            order=['pi_n','pi_sidis','pi_delta'];summary={};steps=[]
            for component in order:
                active=mode=='enabled' and side=='kaon' and component=='pi_sidis'
                bins=[61,63,64,65,66,67,68] if active else list(range(61,69))
                obj=dict(estimator='least_squares',fit_bin_indices=bins,fit_y=[2.]*len(bins),fit_sigma=[1.]*len(bins),fit_template=[1.]*len(bins),weighted_numerator=2.*len(bins),weighted_denominator=float(len(bins)))
                if active:obj['d1a_original_support']=dict(candidate='D1.a',original_support_schema='kaon_sidis_original_support/v1',minimum_fit_bins=2,excluded_bin_indices=[62])
                residual=dict(status='available',edges=deepcopy(seed['mm_edges']),values=[2.]*102,errors=[1.]*102,stored_sumw2=[1.]*102)
                steps.append(dict(component_name=component,pass_index=1,amplitude=2.,n_fit_bins=len(bins),anchor_windows=[[1.06,1.11]],amplitude_windows=[],excluded_windows=[],fit_diagnostics=obj,H_residual_input=residual))
                summary[component]=dict(success=True,amplitude=2.,n_fit_bins=len(bins))
            snap[key]=steps;snap['diagnostics'][side]=dict(success=True,fallback_used=False,fit_order=order,n_passes=1,staged_pass_history=[summary],fit_min=.7,fit_max=1.3)
            if side=='kaon':snap['diagnostics'][side]['pi_delta_signal_protected_fit']=dict(enabled=True,status='success',protected_fit_succeeded=True,fit_variant='lambda_only_protected_fallback',fallback_used=True)
        fits[mode]=dict(d1a_enabled=mode=='enabled',input_unchanged=True,snapshot=snap)
    before=dict(source={'scientific_source_LF_sha256':seed['scientific_source_LF_sha256']},inputs={k:dict(status='available',path=bindings[k]['canonical_path'],bytes=v['bytes'],sha256=v['sha256'],stat_token=deepcopy(bindings[k]['token'])) if v else dict(status='not_configured') for k,v in plan['inputs'].items()})
    return dict(schema=m.SCHEMA,scope=m.SCOPE,seed_sha256=m.digest(seed),plan_digest=m.digest(plan),runtime=dict(kind='PyROOT',ROOT_version='unit record fixture - NOT EXECUTED'),input_bindings=bindings,preservation=dict(before=before,after=deepcopy(before)),templates=payloads,cache_checks=cache,fits=fits,
        d1a_boundary=dict(disabled_bins=list(range(61,69)),enabled_bins=[61,63,64,65,66,67,68],excluded_bins=[62],disabled_sidis=0.,enabled_sidis=.218867157746124,physical_histograms_unchanged=True),cut_state_restored=True,templates_unchanged_after_fits=True,full_analysis_executed=False)


def plan_fixture(seed):return dict(schema=m.SCHEMA,scope=deepcopy(m.SCOPE),seed_sha256=m.digest(seed),inputs=deepcopy(seed['inputs']))


def minimal_figures(result,seed,out):
    """Valid blank test images; never scientific/runtime output evidence."""
    f=out/'figures';f.mkdir()
    def chunk(kind,data):return struct.pack('>I',len(data))+kind+data+struct.pack('>I',zlib.crc32(kind+data)&0xffffffff)
    png=b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',800,400,8,2,0,0,0))+chunk(b'IDAT',zlib.compress((b'\0'+b'\xff'*2400)*400))+chunk(b'IEND',b'')
    names=[c+s for c in m.PIONS for s in ['_direct_t3','_closure']]+['kaon_D1a_modes','pion_D1a_modes']
    for name in names:
        (f/(name+'.png')).write_bytes(png);(f/(name+'.svg')).write_text('<svg xmlns="http://www.w3.org/2000/svg"><text>UNIT FIXTURE ONLY</text></svg>')
    with (out/'fit_scalars.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['D1a_mode','A_n','A_delta','A_sidis','B_n','B_delta','B_sidis']);writer.writerows([[mode,*([0.]*6)] for mode in ['disabled','enabled']])


class CheckerTests(unittest.TestCase):
    def setUp(self):
        self.seed=deepcopy(SEED);self.plan=plan_fixture(self.seed);self.result=record_fixture(self.seed,self.plan)
    def checked(self):return m.check_result(self.result,self.plan,self.seed,m.digest(self.seed))
    def test_retained_record_shaped_fixture_not_ROOT_execution(self):self.assertEqual(self.checked()['status'],'PASS')
    def test_wrong_scope(self):
        self.result['scope']=dict(m.SCOPE,t_index=0)
        with self.assertRaises(ValueError):self.checked()
    def test_wrong_plan_hash(self):
        self.result['plan_digest']='0'*64
        with self.assertRaises(ValueError):self.checked()
    def test_plan_rejects_wrong_input_identity(self):
        self.plan['inputs']['pi_n']['sha256']='0'*64
        with self.assertRaises(ValueError):m.plan_checks(self.plan,self.seed,m.digest(self.seed))
    def test_relative_input_path(self):
        self.plan['inputs']['pi_n']['path']='relative.root'
        with self.assertRaises(ValueError):m.plan_checks(self.plan,self.seed,m.digest(self.seed))
    def test_unknown_optional_source_without_explicit_identity(self):
        self.plan['inputs']['k_sigma0_signal']={'path':'/existing/root','sha256':'0'*64,'bytes':12}
        with self.assertRaises(ValueError):m.plan_checks(self.plan,self.seed,m.digest(self.seed))
    def test_real_retained_support_has_empty_physical_bin(self):
        self.assertEqual(self.seed['experimental']['kaon']['values'][62],0.)
        self.assertTrue(all(x['entries'][62]==0 for x in self.seed['d1a_fixture']['support']['sources'].values()))
    def test_mutation_rejections(self):
        mutations=[('variance',lambda r:r['templates']['pi_n']['support']['parents']['t_bin3']['variances'].__setitem__(1,9.)),
          ('count',lambda r:r['templates']['pi_n']['support']['parents']['t_bin3'].__setitem__('count',-1)),
          ('raw child',lambda r:r['templates']['pi_n']['support']['children']['t_bin3']['phi_bin1']['values'].__setitem__(20,7.)),
          ('normalization',lambda r:r['templates']['pi_n']['parents']['t_bin3']['values'].__setitem__(80,8.)),
          ('axis',lambda r:r['templates']['pi_n']['parents']['t_bin3']['edges'].__setitem__(10,.71)),
          ('geometry',lambda r:r['templates']['pi_n']['identity']['hgcer_selection']['points'].__setitem__(0,[200.,200.])),
          ('definition',lambda r:r['templates']['pi_n']['identity'].__setitem__('population_definition',2)),
          ('weight semantics',lambda r:r['templates']['pi_n']['identity'].__setitem__('weight','renormalized_weight')),
          ('MM selection',lambda r:r['templates']['pi_n']['identity'].__setitem__('mm_selection_window',[.7,1.3])),
          ('cache',lambda r:r['cache_checks']['pi_n'].__setitem__('probe_entries',1)),
          ('clone',lambda r:r['cache_checks']['pi_n'].__setitem__('clone_independent',False)),
          ('unsupported',lambda r:r['templates']['pi_n'].__setitem__('scope_usable',False)),
          ('Sigma substitution',lambda r:r['templates']['k_sigma0_signal'].__setitem__('scope_usable',True)),
          ('amplitude',lambda r:r['fits']['enabled']['snapshot'].__setitem__('A_n',float('nan'))),
          ('fit history',lambda r:r['fits']['enabled']['snapshot'].__setitem__('H_pion_fit_step_overlays',[])),
          ('config',lambda r:r['fits']['enabled']['snapshot']['resolved_subtraction_config']['kaon_nosub'].__setitem__('fit_mode','another_model')),
          ('alignment',lambda r:r['fits']['enabled']['snapshot']['pion_component_alignment'].__setitem__('accepted',False)),
          ('D1 bins',lambda r:r['d1a_boundary'].__setitem__('enabled_bins',list(range(61,69)))),
          ('D1 physical',lambda r:r['d1a_boundary'].__setitem__('physical_histograms_unchanged',False)),
          ('full D1 flag',lambda r:r['fits']['enabled']['snapshot'].__setitem__('H_kaon_fit_step_overlays',[{'component_name':'pi_sidis'}])),
          ('source template mutation',lambda r:r.__setitem__('templates_unchanged_after_fits',False)),
          ('cut state',lambda r:r.__setitem__('cut_state_restored',False)),
          ('full main',lambda r:r.__setitem__('full_analysis_executed',True)),
          ('preservation',lambda r:r['preservation']['after']['source'].__setitem__('HEAD','changed'))]
        for label,mutate in mutations:
            with self.subTest(label=label):
                self.result=record_fixture(self.seed,self.plan);mutate(self.result)
                with self.assertRaises((ValueError,KeyError)):self.checked()
    def test_unknown_retained_config_remains_unknown(self):
        m.compatible({'a':17},{'a':{'status':'unavailable','reason':'not retained'}})
        with self.assertRaises(ValueError):m.compatible({'a':1},{'a':True})
    def test_duplicate_and_nonfinite_JSON_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bad.json'
            for text in ['{"a":1,"a":2}','{"a":NaN}']:
                p.write_text(text)
                with self.assertRaises(ValueError):m.load_json(p)


class BackendWiringTests(unittest.TestCase):
    """Exercise actual backend control flow with labelled fake modules, not ROOT."""
    def fake_modules(self,unexpected_sigma=False):
        seed=deepcopy(SEED);temp=tempfile.TemporaryDirectory();self.addCleanup(temp.cleanup)
        for name,item in seed['inputs'].items():
            if item:
                f=Path(temp.name)/(name+'.root');f.write_bytes(('UNIT '+name).encode());item.update(path=str(f),bytes=f.stat().st_size,sha256=m.file_hash(f))
        plan=plan_fixture(seed);record=record_fixture(seed,plan);calls=[]
        class Hist:
            def __init__(self,name,title,n,edges):self.edges=list(edges);self.values=[0.]*(n+2);self.errors=[0.]*(n+2)
            def SetDirectory(self,*a):pass
            def Sumw2(self):pass
            def SetBinContent(self,i,v):self.values[i]=v
            def SetBinError(self,i,v):self.errors[i]=v
            def SetEntries(self,*a):pass
            def GetBinContent(self,i):return self.values[i]
        def hist(h):return dict(status='available',edges=deepcopy(h.edges),values=deepcopy(h.values),errors=deepcopy(h.errors),stored_sumw2=[x*x for x in h.errors])
        def restored(saved):
            h=Hist('unit','unit',100,saved['edges']);h.values=deepcopy(saved['values']);h.errors=deepcopy(saved['errors']);return h
        def snapshot(x):
            if isinstance(x,Hist):return hist(x)
            if isinstance(x,dict):return {k:snapshot(v) for k,v in x.items()}
            if isinstance(x,(list,tuple)):return [snapshot(v) for v in x]
            return x
        root=types.ModuleType('ROOT');root.__file__=__file__;root.TH1D=Hist;root.gROOT=types.SimpleNamespace(SetBatch=lambda *a:None,GetVersion=lambda:'UNIT TEST - NOT ROOT')
        cuts=types.ModuleType('apply_cuts');cuts.current='experimental'
        shapes=types.ModuleType('pion_component_shapes');shapes._SIMC_COMPONENT_SHAPE_CACHE={}
        def loader(path,cfg,phi,particle,component,**kw):
            cuts.current='SIMC temporary'
            key=(path,component,tuple(kw['phi_bins']))
            if key not in shapes._SIMC_COMPONENT_SHAPE_CACHE:
                payload=record['templates'][component]
                obj=dict(event_population_identity=deepcopy(payload['identity']),event_population_support=deepcopy(payload['support']),diagnostics={},setting_shape_full=restored(payload['setting']),t_shapes={k:restored(v) for k,v in payload['parents'].items()},binned_shapes={k:{p:restored(v) for p,v in children.items()} for k,children in payload['children'].items()})
                obj['event_population_identity']['phi_edges']=kw['phi_bins']
                shapes._SIMC_COMPONENT_SHAPE_CACHE[key]=obj
            return deepcopy(shapes._SIMC_COMPONENT_SHAPE_CACHE[key])
        shapes._get_root_file_cache_token=m.file_token
        shapes.resolve_simc_tree_name=lambda cfg:'h10'
        shapes.resolve_kaon_simc_signal_root_filename=lambda path,*a:(path,{})
        shapes.load_pion_simc_component_shape=loader
        shapes.load_kaon_simc_signal_shape=lambda path,cfg,phi,**kw:loader(path,cfg,phi,'kaon','k_lambda_signal',**kw)
        shapes.load_kaon_simc_sigma0_shape=lambda path,*a,**kw:dict(diagnostics={})
        shapes.resolve_kaon_simc_sigma0_root_filename=lambda path,*a:('undeclared.root',{}) if unexpected_sigma else (path,{})
        fits=types.ModuleType('pion_component_fits')
        fits.resolve_scope_single_shape=lambda p,*a,**kw:(p.get('t_shapes') or {}).get('t_bin3')
        fits.resolve_scope_component_shapes=lambda p,*a,**kw:{c:v['t_shapes']['t_bin3'] for c,v in p['components'].items()}
        def fitter(hp,hk,templates,inp,scope,**kw):
            calls.append(dict(scope=scope,mode=inp['d1a_kaon_sidis_support_enabled'],alignment=deepcopy(kw['pion_component_alignment']),support=deepcopy(kw['kaon_experimental_support'])))
            return deepcopy(record['fits']['enabled' if inp['d1a_kaon_sidis_support_enabled'] else 'disabled']['snapshot'])
        fits.build_particle_subtraction_component_result=fitter
        def boundary(*a,**kw):
            enabled=kw['kaon_sidis_experimental_support'] is not None
            diag=dict(fit_bin_indices=[61,63,64,65,66,67,68] if enabled else list(range(61,69)))
            if enabled:diag['d1a_original_support']={'excluded_bin_indices':[62]}
            return ({'pi_n':.182202244063977,'pi_sidis':.218867157746124 if enabled else 0.},{},None,[{},dict(fit_diagnostics=diag)])
        fits._run_staged_component_pass=boundary
        audit=types.ModuleType('pion_background_chain_audit');audit.snapshot=snapshot;audit.histogram=hist
        numpy=types.ModuleType('numpy');numpy.__version__='UNIT';scipy=types.ModuleType('scipy');scipy.__version__='UNIT'
        modules=dict(ROOT=root,apply_cuts=cuts,pion_component_shapes=shapes,pion_component_fits=fits,pion_background_chain_audit=audit,numpy=numpy,scipy=scipy)
        return seed,plan,modules,calls
    def test_actual_backend_glue_calls_scoped_fitter_and_both_D1a_modes(self):
        seed,plan,modules,calls=self.fake_modules()
        with patch.dict(sys.modules,modules):result=m.compute_root(plan,seed)
        result.update(schema=m.SCHEMA,scope=m.SCOPE,seed_sha256=m.digest(seed),plan_digest=m.digest(plan),preservation=record_fixture(seed,plan)['preservation'])
        self.assertEqual(m.check_result(result,plan,seed,m.digest(seed))['status'],'PASS')
        self.assertEqual([x['mode'] for x in calls],[False,True]);self.assertTrue(all(x['scope']=='t_bin3' for x in calls))
        self.assertTrue(all(x['alignment']==seed['alignment'] and x['support']==seed['d1a_fixture']['support'] for x in calls))
    def test_substituted_Lambda_is_rejected_before_any_loader_or_fit(self):
        seed,plan,modules,calls=self.fake_modules()
        modules['pion_component_shapes'].resolve_kaon_simc_signal_root_filename=lambda *a:(plan['inputs']['pi_n']['path'],{})
        with patch.dict(sys.modules,modules):
            with self.assertRaisesRegex(ValueError,'resolved_input_substitution'):m.compute_root(plan,seed)
        self.assertEqual(calls,[]);self.assertEqual(len(modules['pion_component_shapes']._SIMC_COMPONENT_SHAPE_CACHE),0)
    def test_undeclared_Sigma0_is_rejected_before_any_loader_or_fit(self):
        seed,plan,modules,calls=self.fake_modules(True)
        with patch.dict(sys.modules,modules):
            with self.assertRaisesRegex(ValueError,'undeclared_Sigma0'):m.compute_root(plan,seed)
        self.assertEqual(calls,[]);self.assertEqual(len(modules['pion_component_shapes']._SIMC_COMPONENT_SHAPE_CACHE),0)


class OwnerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name);self.repo=self.root/'repo';self.repo.mkdir();self.seed=deepcopy(SEED)
        for name in m.ROLES:
            if self.seed['inputs'][name] is not None:
                p=self.root/(name+'.root');p.write_bytes(('UNIT INPUT '+name).encode());self.seed['inputs'][name]=dict(path=str(p),bytes=p.stat().st_size,sha256=m.file_hash(p))
        self.seed_path=self.root/'seed.json';m.write_json(self.seed_path,self.seed);self.plan=plan_fixture(self.seed);self.plan['seed_sha256']=m.file_hash(self.seed_path);self.plan_path=self.root/'plan.json';m.write_json(self.plan_path,self.plan)
        self.out=self.root/'output'
    def state(self,*a):return {'HEAD':'UNIT','scientific_source_LF_sha256':deepcopy(self.seed['scientific_source_LF_sha256'])}
    def backend(self,plan,seed):return record_fixture(seed,plan)
    def run_owner(self,backend=None,renderer=minimal_figures):
        return m.run_owner(self.plan_path,m.file_hash(self.plan_path),self.out,self.seed_path,self.repo,backend or self.backend,renderer,self.state)
    def test_success_owner_and_independent_checker_with_injected_backend(self):
        self.run_owner();self.assertEqual(m.verify_directory(self.out,self.seed_path)['status'],'PASS')
        self.assertEqual(m.load_json(self.out/'owner.json')['status'],'passed')
    def test_existing_output_never_overwritten(self):
        self.out.mkdir();(self.out/'sentinel').write_text('original')
        with self.assertRaises(ValueError):self.run_owner()
        self.assertEqual((self.out/'sentinel').read_text(),'original')
    def test_inside_repository_output_rejected(self):
        self.out=self.repo/'diagnostic'
        with self.assertRaises(ValueError):self.run_owner()
        self.assertFalse(self.out.exists())
    def test_missing_input_before_backend(self):
        Path(self.plan['inputs']['pi_n']['path']).unlink()
        with self.assertRaises(FileNotFoundError):self.run_owner()
        self.assertEqual(m.load_json(self.out/'owner.json')['status'],'failed');self.assertFalse((self.out/'receipt.json').exists())
    def test_input_mutation_during_backend(self):
        def mutate(p,s):
            Path(p['inputs']['pi_n']['path']).write_bytes(b'changed');return record_fixture(s,p)
        with self.assertRaises(ValueError):self.run_owner(mutate)
        self.assertEqual(m.load_json(self.out/'owner.json')['status'],'failed');self.assertFalse((self.out/'receipt.json').exists())
    def test_backend_failure_stops_success(self):
        def fail(*a):raise RuntimeError('original ROOT failure')
        with self.assertRaisesRegex(RuntimeError,'original ROOT failure'):self.run_owner(fail)
        status=m.load_json(self.out/'owner.json');self.assertEqual(status['status'],'failed');self.assertEqual(status['preservation']['before'],status['preservation']['after'])
    def test_render_failure_stops_success(self):
        def fail(*a):raise RuntimeError('figure failure')
        with self.assertRaises(RuntimeError):self.run_owner(renderer=fail)
        self.assertEqual(m.load_json(self.out/'owner.json')['status'],'failed');self.assertFalse((self.out/'receipt.json').exists())
    def test_checker_rejects_tampered_artifact(self):
        self.run_owner();(self.out/'result.json').write_text('{}')
        with self.assertRaises(ValueError):m.verify_directory(self.out,self.seed_path)
    def test_rehashed_wrong_numbers_still_rejected(self):
        self.run_owner();r=m.load_json(self.out/'result.json');r['d1a_boundary']['enabled_bins']=[62];m.write_json(self.out/'result.json',r)
        manifest=m.load_json(self.out/'SHA256SUMS.json');manifest['result.json']=dict(bytes=(self.out/'result.json').stat().st_size,sha256=m.file_hash(self.out/'result.json'));m.write_json(self.out/'SHA256SUMS.json',manifest)
        with self.assertRaises(ValueError):m.verify_directory(self.out,self.seed_path)
    def test_complete_true_alone_never_passes(self):
        self.out.mkdir();m.write_json(self.out/'owner.json',{'complete':True})
        with self.assertRaises(FileNotFoundError):m.verify_directory(self.out,self.seed_path)
    def test_failure_after_receipt_is_written_removes_success_status(self):
        with patch.object(m,'verify_directory',side_effect=ValueError('final artifact failure')):
            with self.assertRaises(ValueError):self.run_owner()
        self.assertEqual(m.load_json(self.out/'receipt.json')['status'],'failed');self.assertEqual(m.load_json(self.out/'owner.json')['status'],'failed')
    def test_cut_globals_restored_on_exception(self):
        import types
        cut=types.ModuleType('cuts');cut.current='experimental';before=dict(vars(cut))
        with self.assertRaises(RuntimeError):
            with m.preserved_cuts(cut):cut.current='temporary';cut.new_field=1;raise RuntimeError()
        self.assertEqual(dict(vars(cut)),before)
    def test_cli_help_does_not_import_ROOT(self):
        p=subprocess.run([sys.executable,'-B',str(REPO/'testing/run_d2a_root_integration.py'),'--help'],capture_output=True,text=True)
        self.assertEqual(p.returncode,0,p.stderr);self.assertIn('prepare',p.stdout)
    def test_prepare_plan_no_ROOT_or_input_read(self):
        p=self.root/'prepared.json'
        with patch.object(m,'SEED',self.seed_path),patch.object(m,'REPO',self.repo):
            self.assertEqual(m.main(['prepare','--output',str(p)]),0)
        self.assertTrue(p.is_file())



class ReviewGapRegressions(unittest.TestCase):
    setUp=CheckerTests.setUp
    checked=CheckerTests.checked
    """Exact rejection probes from the failed actual-diff review."""
    def test_failed_fit_status_and_fallback_never_pass(self):
        for key,value in [('fit_status_pion','fallback'),('fit_status_kaon','fallback'),('fallback_used',True)]:
            with self.subTest(key=key):
                self.result=record_fixture(self.seed,self.plan);self.result['fits']['enabled']['snapshot'][key]=value
                with self.assertRaises(ValueError):self.checked()
    def test_component_failure_missing_history_or_unavailable_raw_rejected(self):
        mutations=[lambda s:s['diagnostics']['kaon']['staged_pass_history'][0]['pi_delta'].__setitem__('success',False),
                   lambda s:s['staged_amplitudes_raw']['kaon_nosub'].__setitem__('pi_sidis',{'status':'unavailable'}),
                   lambda s:s['H_kaon_fit_step_overlays'].pop(),
                   lambda s:s['diagnostics']['pion'].__setitem__('success',False)]
        for fn in mutations:
            self.result=record_fixture(self.seed,self.plan);fn(self.result['fits']['enabled']['snapshot'])
            with self.assertRaises((ValueError,KeyError)):self.checked()
    def test_empty_D1_marker_or_marker_on_neutron_rejected(self):
        for component in ['pi_n','pi_sidis']:
            self.result=record_fixture(self.seed,self.plan)
            step=next(x for x in self.result['fits']['enabled']['snapshot']['H_kaon_fit_step_overlays'] if x['component_name']==component)
            step['fit_diagnostics']['d1a_original_support']={}
            with self.assertRaises((ValueError,KeyError)):self.checked()
    def test_wrong_D1_exclusion_and_objective_arithmetic_rejected(self):
        for field,value in [('weighted_denominator',99.),('fit_sigma',[0.]*7),('fit_bin_indices',[62]*7)]:
            self.result=record_fixture(self.seed,self.plan);obj=self.result['fits']['enabled']['snapshot']['H_kaon_fit_step_overlays'][1]['fit_diagnostics'];obj[field]=value
            with self.assertRaises(ValueError):self.checked()
        self.result=record_fixture(self.seed,self.plan);self.result['fits']['enabled']['snapshot']['H_kaon_fit_step_overlays'][1]['fit_diagnostics']['d1a_original_support']['excluded_bin_indices']=[]
        with self.assertRaises(ValueError):self.checked()
    def test_successful_Lambda_only_policy_allowed_but_failed_Delta_rejected(self):
        self.assertEqual(self.checked()['status'],'PASS')
        protected=self.result['fits']['enabled']['snapshot']['diagnostics']['kaon']['pi_delta_signal_protected_fit']
        protected.update(status='lambda_preservation_rejected',protected_fit_succeeded=False,fit_variant='zero_pi_delta_failure')
        with self.assertRaisesRegex(ValueError,'failed_protected_Delta_policy'):self.checked()
    def test_wrong_loader_source_tree_selection_rejected(self):
        for key,value in [('source',['/unrelated.root',1,1,1,1,1]),('tree','wrong'),('selection',{'Q2':'wrong'})]:
            self.result=record_fixture(self.seed,self.plan);self.result['templates']['k_lambda_signal']['identity'][key]=value
            with self.assertRaises(ValueError):self.checked()
    def test_wrong_resolved_path_or_stat_token_rejected(self):
        self.result['input_bindings']['k_lambda_signal']['canonical_path']='/unrelated.root'
        with self.assertRaises(ValueError):self.checked()
        self.result=record_fixture(self.seed,self.plan);self.result['templates']['pi_n']['identity']['source'][1]=999
        with self.assertRaises(ValueError):self.checked()


class ResolverBindingRegressions(unittest.TestCase):
    def actual_resolver(self):
        path=REPO/'src/utility/pion_component_shapes.py';tree=ast.parse(path.read_text())
        node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='resolve_kaon_simc_signal_root_filename')
        ns={'os':os};exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),ns);return ns[node.name]
    def fixture(self,directory):
        seed=deepcopy(SEED)
        for name,item in seed['inputs'].items():
            if item:
                p=directory/(name+'.root');p.write_bytes(('UNIT '+name).encode());item.update(path=str(p),bytes=p.stat().st_size,sha256=m.file_hash(p))
        shapes=types.SimpleNamespace(resolve_kaon_simc_signal_root_filename=self.actual_resolver(),resolve_kaon_simc_sigma0_root_filename=lambda path,*a:(path,{}),resolve_simc_tree_name=lambda cfg:'h10',_get_root_file_cache_token=m.file_token)
        return seed,plan_fixture(seed),shapes
    def test_actual_resolver_substitution_rejected_before_any_loader(self):
        with tempfile.TemporaryDirectory() as d:
            directory=Path(d);seed,plan,shapes=self.fixture(directory)
            declared=directory/'Prod_Coin_same_hash_model.root';declared.write_bytes(b'UNIT declared model')
            alternate=directory/'Center_kaon_Simc_Q4p4W2p74_highe.root';alternate.write_bytes(b'UNIT different model')
            plan['inputs']['k_lambda_signal']['path']=str(declared)
            resolved,_=shapes.resolve_kaon_simc_signal_root_filename(str(declared),seed['config'],'Center');self.assertEqual(Path(resolved),alternate)
            with self.assertRaisesRegex(ValueError,'resolved_input_substitution'):m.resolve_input_bindings(shapes,plan,seed['config'])
    def test_same_file_canonical_alias_permitted(self):
        with tempfile.TemporaryDirectory() as d:
            seed,plan,shapes=self.fixture(Path(d));original=plan['inputs']['k_lambda_signal']['path']
            (Path(d)/'alias').mkdir()
            shapes.resolve_kaon_simc_signal_root_filename=lambda *a:(str(Path(d)/'alias'/'..'/Path(original).name),{})
            bindings=m.resolve_input_bindings(shapes,plan,seed['config']);self.assertEqual(bindings['k_lambda_signal']['canonical_path'],str(Path(original).resolve()))
    def test_successful_explicit_optional_input_bound(self):
        with tempfile.TemporaryDirectory() as d:
            seed,plan,shapes=self.fixture(Path(d));p=Path(d)/'optional.root';p.write_bytes(b'UNIT optional')
            plan['inputs']['k_sigma0_signal']=dict(path=str(p),bytes=p.stat().st_size,sha256=m.file_hash(p),explicit_existing_source=True)
            self.assertIn('k_sigma0_signal',m.resolve_input_bindings(shapes,plan,seed['config']))


class FailedFitEvidenceTests(unittest.TestCase):
    setUp=OwnerTests.setUp
    state=OwnerTests.state
    backend=OwnerTests.backend
    run_owner=OwnerTests.run_owner
    def test_rejected_fit_record_retained_without_accepted_receipt(self):
        def bad(plan,seed):
            result=record_fixture(seed,plan);result['fits']['enabled']['snapshot']['fit_status_kaon']='fallback';return result
        with self.assertRaises(ValueError):self.run_owner(bad)
        self.assertEqual(m.load_json(self.out/'result.json')['fits']['enabled']['snapshot']['fit_status_kaon'],'fallback')
        self.assertEqual(m.load_json(self.out/'owner.json')['status'],'failed');self.assertFalse((self.out/'receipt.json').exists())


if __name__=='__main__':unittest.main()
