"""D2.a deterministic production-function tests with ROOT adapters; no runtime proof."""
from pathlib import Path
import ast,copy,importlib.util,math,sys,types,unittest,tempfile,os
from unittest.mock import patch
import numpy as np
REPO=Path(__file__).resolve().parents[1]
if not (REPO/'src').exists():REPO=Path('C:/Users/trott/Documents/Programs/lt_analysis')
sys.path[:0]=[str(REPO/'src/cuts'),str(REPO/'src/utility'),str(REPO/'src/binning')]
from testing.test_pion_component_dynamic_alignment import fits

class Hist:
    def __init__(self,name='h',title='',n=100,lo=.7,hi=1.3):
        self.name=name;self.edges=np.linspace(lo,hi,n+1);self.values=np.zeros(n+2);self.var=np.zeros(n+2);self.entries=0;self.sumw2=False
    def Sumw2(self):self.sumw2=True
    def SetDirectory(self,*a):pass
    def Fill(self,x,w=1):
        assert self.sumw2,'Sumw2 must precede fills'
        i=int(np.searchsorted(self.edges,x,side='right'));self.values[i]+=w;self.var[i]+=w*w;self.entries+=1;return i
    def GetNbinsX(self):return len(self.edges)-1
    def GetBinContent(self,i):return self.values[i]
    def GetBinError(self,i):return math.sqrt(self.var[i])
    def GetEntries(self):return self.entries
    def GetName(self):return self.name
    def GetXaxis(self):return self
    def GetXmin(self):return self.edges[0]
    def GetXmax(self):return self.edges[-1]
    def Integral(self):return float(self.values[1:-1].sum())
    def Scale(self,f):self.values*=f;self.var*=f*f
    def Clone(self,name):h=copy.deepcopy(self);h.name=name;return h

class Joint:
    def __init__(self,*a):self.rows=[];self.sumw2=False
    def SetDirectory(self,*a):pass
    def Sumw2(self):self.sumw2=True
    def Fill(self,x,y,w):assert self.sumw2;self.rows.append((x,y,w))
    def GetZaxis(self):return self
    def SetRangeUser(self,*a):pass

def event(mm=1.,t=.45,phi=0.,weight=1.,selected=True):
    d=dict(missmass=mm,t=-t,phipq=phi*math.pi/180.,iter_weight=weight,selected=selected,
           ssdelta=0.,ssxptar=0.,ssyptar=0.,hsdelta=0.,hsxptar=0.,hsyptar=0.,Q2=4.4,W=2.74,phgcer_x_det=100.,phgcer_y_det=100.)
    return types.SimpleNamespace(**d)

def utility_functions():
    tree=ast.parse((REPO/'src/utility/utility.py').read_text(encoding='utf-8-sig'))
    names=['normalize_hist_to_unit_area','prune_hist','TH1D_to_TH2D']
    ns=dict(math=math,np=np,ROOT=types.SimpleNamespace(TH2D=Joint),_find_outlier_bins_by_frac_err=lambda *a,**k:[1])
    exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names],type_ignores=[]),'actual_utility_functions','exec'),ns);return ns

class PopulationHole:
    """Detached geometry adapter; default is outside already-selected records."""
    def __init__(self,low=1e100,high=2e100):
        self.points=[(low,-200.),(high,-200.),(high,200.),(low,200.),(low,-200.)]
    def GetN(self):return len(self.points)
    def GetX(self):return [p[0] for p in self.points]
    def GetY(self):return [p[1] for p in self.points]
    def IsInside(self,x,y):return self.points[0][0]<=x<=self.points[1][0] and -200.<=y<=200.


def loader_module(records,missing=()):
    class Tree:
        def GetEntries(self):return len(records)
        def __iter__(self):return iter(records)
        def GetBranch(self,name):return name not in missing
        GetLeaf=GetBranch
    class File:
        def IsZombie(self):return False
        def Get(self,*a):return Tree()
        def Close(self):pass
    root=types.ModuleType('ROOT');root.TH1D=Hist;root.TFile=types.SimpleNamespace(Open=lambda *a:File())
    lt=types.ModuleType('ltsep');lt.Misc=types.SimpleNamespace(progressBar=lambda *a,**k:None)
    utility=types.ModuleType('utility');utility.normalize_hist_to_unit_area=utility_functions()['normalize_hist_to_unit_area']
    cuts=types.ModuleType('apply_cuts');cuts.set_val=lambda *a:None;cuts.apply_simc_sub_cuts=lambda e,*a:e.selected;cuts.apply_simc_cuts=lambda e,lo,hi:e.selected and lo<=e.missmass<hi
    hole=types.ModuleType('hgcer_hole');hole.apply_HGCer_hole_cut=lambda *a:PopulationHole()
    with patch.dict(sys.modules,{'ROOT':root,'ltsep':lt,'utility':utility,'apply_cuts':cuts,'hgcer_hole':hole}):
        spec=importlib.util.spec_from_file_location('d2a_actual_shapes',REPO/'src/utility/pion_component_shapes.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

class EventPopulationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.path=Path(self.tmp.name)/'simc.root';self.path.write_bytes(b'original source identity')
        self.inp=dict(Q2='4p4',W='2p74',EPSSET='high',POL='1',ParticleType='kaon',mm_min=.7,mm_max=1.3,tmin=.4,tmax=.9,cut_poly=[[4.,2.],[5.,2.],[5.,3.]],bg_opt_mm_plot_nbins=100)
    def tearDown(self):self.tmp.cleanup()
    def load(self,m,c='pi_n',**kw):
        return m.load_pion_simc_component_shape(str(self.path),self.inp,'Center','kaon',c,t_bins=[.4,.6,.9],phi_bins=[-180.,0.,180.],**kw)
    def test_all_five_channels_direct_event_parent_and_child_moments(self):
        records=[event(mm=1.0,phi=-90.) for _ in range(9)]+[event(mm=1.2,phi=90.)]+[event(mm=.6,phi=-90.,weight=-2),event(mm=1.4,phi=180.,weight=3),event(t=.6,phi=90.,weight=2),event(t=.9),event(selected=False)]
        for channel in ['pi_n','pi_delta','pi_sidis','k_lambda_signal','k_sigma0_signal']:
            with self.subTest(channel=channel):
                m=loader_module(records);p=self.load(m,channel);sup=p['event_population_support'];a=sup['parents']['t_bin1'];children=sup['children']['t_bin1'];outside=sup['outside_phi']['t_bin1']
                np.testing.assert_allclose(a['values'],np.sum([h['values'] for h in children.values()],axis=0)+outside['values']);np.testing.assert_allclose(a['variances'],np.sum([h['variances'] for h in children.values()],axis=0)+outside['variances'])
                self.assertEqual(a['count'],12);self.assertEqual(a['absolute'],15.);self.assertEqual(a['sumw2'],23.);self.assertEqual(a['values'][0],-2);self.assertEqual(a['values'][-1],3);self.assertEqual(outside['count'],1);self.assertEqual(sup['outside_t']['count'],1);np.testing.assert_allclose(sup['setting']['values'],np.array([x['values'] for x in sup['parents'].values()]).sum(0)+np.array(sup['outside_t']['values']))
                h=fits._resolve_component_scope_hist(p,channel,'t_bin1',0);np.testing.assert_allclose(h.values,np.array(a['values'])/10.)
                self.assertIs(h,fits._resolve_single_scope_hist(p,'t_bin1',0));self.assertIsNone(fits._resolve_single_scope_hist(p,'t_bin2phi_bin1',1,0));self.assertIsNone(fits._resolve_single_scope_hist(p,'missing',4,0))
                self.assertIs(fits._resolve_single_scope_hist(p,'setting-wide'),p['setting_shape_full'])
    def test_unequal_children_cannot_manufacture_parent(self):
        m=loader_module([event(mm=1.,phi=-90.) for _ in range(9)]+[event(mm=1.2,phi=90.)]);p=self.load(m);a=p['t_shapes']['t_bin1'];children=p['binned_shapes']['t_bin1'];wrong=sum((h.values for h in children.values()))/2
        self.assertGreater(np.abs(a.values-wrong).sum(),.7)
        with self.assertRaisesRegex(RuntimeError,'cannot_be_reconstructed'):fits._sum_hist_list_to_unit_area(list(children.values()),'invalid')
    def test_signed_cancellation_zero_weight_and_sparse_not_replaced(self):
        p=self.load(loader_module([event(weight=2),event(weight=-2),event(weight=0),event(t=.7)]));sup=p['event_population_support']['parents']['t_bin1'];self.assertEqual(sup['count'],3);self.assertEqual(sup['normal_sumw'],0);self.assertAlmostEqual(sup['normal_sumw2'],8);self.assertIsNone(fits._resolve_single_scope_hist(p,'t_bin1',0));self.assertIsNotNone(fits._resolve_single_scope_hist(p,'setting-wide'))
        tiny=self.load(loader_module([event(weight=1e-14)]))
        self.assertEqual(tiny['event_population_support']['parents']['t_bin1']['status'],'available')
        self.assertIsNotNone(fits._resolve_single_scope_hist(tiny,'t_bin1',0))

    def test_cache_selection_exact_edges_source_and_clone_isolation(self):
        m=loader_module([event(),event(t=.7)]);p=self.load(m);q=self.load(m);p['t_shapes']['t_bin1'].values[2]=999;p['event_population_support']['parents']['t_bin1']['count']=999
        self.assertNotEqual(q['t_shapes']['t_bin1'].values[2],999);self.assertNotEqual(q['event_population_support']['parents']['t_bin1']['count'],999)
        n=len(m._SIMC_COMPONENT_SHAPE_CACHE);self.inp['poly_points']=[[4.,2.],[5.,2.],[5.,3.]];self.load(m);self.assertEqual(len(m._SIMC_COMPONENT_SHAPE_CACHE),n+1);self.inp['poly_points'][0][0]=3.;self.load(m);self.assertEqual(len(m._SIMC_COMPONENT_SHAPE_CACHE),n+2)
        self.assertNotEqual(m._freeze_bin_edges_for_cache([0.,.500000001,1.]),m._freeze_bin_edges_for_cache([0.,.500000002,1.]))
        old=q['event_population_identity']['source'];self.path.write_bytes(b'changed source identity payload');r=self.load(m);self.assertNotEqual(old,r['event_population_identity']['source'])
    def test_invalid_edges_missing_branches_nonfinite_and_missing_sources(self):
        m=loader_module([event(weight=float('nan')),event(),event(mm=float('nan'))]);p=self.load(m);self.assertEqual(p['event_population_support']['setting']['count'],1)
        for edges in [[.4,.4,.9],[.9,.4],[.4,float('nan'),.9]]:
            with self.assertRaisesRegex(ValueError,'edges'):m.load_pion_simc_component_shape(str(self.path),self.inp,'Center','kaon','pi_n',t_bins=edges)
        p=self.load(loader_module([event()],missing=('phipq',)));self.assertIsNone(fits._resolve_single_scope_hist(p,'t_bin1',0));self.assertIn('phipq',p['diagnostics']['missing_required_branches'])
        self.path.unlink();p=self.load(m);self.assertIsNone(fits._resolve_single_scope_hist(p,'t_bin1',0))
    def test_parent_only_scope_and_absent_provenance(self):
        m=loader_module([event()]);p=m.load_pion_simc_component_shape(str(self.path),self.inp,'Center','kaon','pi_n',t_bins=[.4,.9]);self.assertIsNotNone(fits._resolve_single_scope_hist(p,'t_bin1',0));self.assertIsNone(fits._resolve_single_scope_hist(dict(p,event_population_identity=None),'t_bin1',0))
        with self.assertRaisesRegex(RuntimeError,'indices_required'):fits._resolve_single_scope_hist(p,'t_bin1')
        with self.assertRaisesRegex(RuntimeError,'t_index'):fits._resolve_single_scope_hist(p,'scope',None,0)
    def test_pruning_never_changes_signed_physical_population(self):
        u=utility_functions();h=Hist();h.Sumw2();h.Fill(1.,-2);h.Fill(1.,2);h.Fill(1.2,-1);before=copy.deepcopy(h.__dict__);self.assertEqual(u['prune_hist'](h,threshold=10000),1)
        np.testing.assert_equal(h.values,before['values']);np.testing.assert_equal(h.var,before['var']);self.assertEqual(h.entries,before['entries'])
    def test_joint_population_requires_real_pairs_not_marginals(self):
        u=utility_functions();h=Hist()
        with self.assertRaisesRegex(ValueError,'paired'):u['TH1D_to_TH2D'](h,h)
        pairs=[(.8,1.2,2.),(1.1,.9,-1.),(1.4,.6,0.)];out=u['TH1D_to_TH2D'](h,h,paired_events=pairs);self.assertEqual(out.rows,pairs)
        with self.assertRaisesRegex(ValueError,'nonfinite'):u['TH1D_to_TH2D'](h,h,paired_events=[(float('nan'),1.,1.)])
    def test_shared_fit_entry_rejects_missing_scoped_population(self):
        with self.assertRaisesRegex(RuntimeError,'unsupported_event_population'):
            fits.build_particle_subtraction_component_result(None,None,{},
                {'particle_subtraction_mode':'simc_shape_components'},'t_bin1')

    def test_parent_identity_preserves_phi_only_adaptation(self):
        tree=ast.parse((REPO/'src/cuts/pion_t_bin_parents.py').read_text())
        node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_parent_template_population_identity')
        ns={'deepcopy':copy.deepcopy};exec(compile(ast.Module(body=[node],type_ignores=[]),'actual_parent_identity','exec'),ns)
        fn=ns[node.name];p={'event_population_identity':{'schema':2,'source':['source.root',1], 't_edges':[.4,.9],'phi_edges':[-180.,0.,180.], 'weight':'iter_weight'}}
        a=fn({'components':{'pi_n':p}},p,None);q=copy.deepcopy(p);q['event_population_identity']['phi_edges']=[-180.,180.]
        self.assertEqual(a,fn({'components':{'pi_n':q}},q,None))
        q['event_population_identity']['t_edges']=[.4,.6,.9];self.assertNotEqual(a,fn({'components':{'pi_n':q}},q,None))
        self.assertEqual(p['event_population_identity']['phi_edges'],[-180.,0.,180.])

    def test_effective_hole_mm_grid_and_particle_semantics_bind_reuse(self):
        tree=ast.parse((REPO/'src/cuts/pion_t_bin_parents.py').read_text())
        helper=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_parent_template_population_identity')
        ns={'deepcopy':copy.deepcopy};exec(compile(ast.Module(body=[helper],type_ignores=[]),'actual_parent_identity','exec'),ns)
        identity=lambda p:ns[helper.name]({'components':{c:p for c in ['pi_n','pi_delta','pi_sidis']}},p,None)
        reuse_test=next(n.test for n in ast.walk(tree) if isinstance(n,ast.If) and isinstance(n.test,ast.BoolOp) and any(isinstance(x,ast.Subscript) and isinstance(x.value,ast.Name) and x.value.id=='template_population_identity' for x in ast.walk(n.test)))
        code=compile(ast.Expression(body=reuse_test),'actual_parent_reuse_condition','eval')
        def reusable(old,new):return eval(code,dict(template_population_identity=identity(new),inp_dict={},existing_parents=({'template_population_identity':identity(old)},),expected_count=1))
        rows=[event(mm=1.),event(mm=1.2)];rows[0].phgcer_x_det=-100.
        m=loader_module(rows)
        a=self.load(m,hgcer_cutg=PopulationHole(-200.,0.));b=self.load(m,hgcer_cutg=PopulationHole(0.,200.))
        self.assertEqual(a['event_population_support']['parents']['t_bin1']['count'],1)
        self.assertEqual(b['event_population_support']['parents']['t_bin1']['count'],1)
        self.assertAlmostEqual(np.abs(a['t_shapes']['t_bin1'].values-b['t_shapes']['t_bin1'].values).sum(),2.)
        self.assertNotEqual(identity(a),identity(b));self.assertFalse(reusable(a,b));self.assertTrue(reusable(a,a))
        self.inp['bg_opt_mm_plot_nbins']=200;c=self.load(m,hgcer_cutg=PopulationHole(-200.,0.))
        self.assertEqual(c['t_shapes']['t_bin1'].GetNbinsX(),200);self.assertFalse(reusable(a,c))
        for key,value in [('bg_opt_mm_plot_min',.6),('bg_opt_mm_plot_max',1.6)]:
            self.inp[key]=value;d=self.load(m,hgcer_cutg=PopulationHole(-200.,0.));self.assertFalse(reusable(c,d))
        unknown=self.load(m,hgcer_cutg=types.SimpleNamespace(IsInside=lambda *a:False))
        self.assertFalse(identity(unknown)['reuse_safe']);self.assertFalse(reusable(unknown,unknown))
        legacy=copy.deepcopy(a);legacy['event_population_identity'].pop('population_definition')
        self.assertFalse(reusable(legacy,legacy))
        phi=copy.deepcopy(a);phi['event_population_identity']['phi_edges']=[-180.,-90.,0.,90.,180.]
        self.assertTrue(reusable(a,phi))
        pion=m.load_pion_simc_component_shape(str(self.path),self.inp,'Center','pion','pi_n',t_bins=[.4,.6,.9])
        self.assertEqual(pion['event_population_identity']['hgcer_selection'],{'applied':False,'selection':'no_hgcer_hole'})
        self.assertFalse(reusable(a,pion))

    def test_auto_hole_geometry_changes_invalidate_loader_cache(self):
        m=loader_module([event()]);p=self.load(m);n=len(m._SIMC_COMPONENT_SHAPE_CACHE)
        m.apply_HGCer_hole_cut=lambda *a:PopulationHole(90.,110.)
        q=self.load(m);self.assertEqual(len(m._SIMC_COMPONENT_SHAPE_CACHE),n+1)
        self.assertNotEqual(p['event_population_identity'],q['event_population_identity'])
        self.assertEqual(q['event_population_support']['setting']['count'],0)
        self.assertEqual(p['event_population_support']['setting']['count'],1)
        m.apply_HGCer_hole_cut=lambda *a:types.SimpleNamespace(IsInside=lambda *a:False)
        self.load(m);self.load(m);self.assertEqual(len(m._SIMC_COMPONENT_SHAPE_CACHE),n+1)

    def test_scoped_simc_loads_restore_live_cut_state_on_success_and_failure(self):
        tree=ast.parse((REPO/'src/cuts/rand_sub.py').read_text())
        block=next(n for n in ast.walk(tree) if isinstance(n,ast.Try) and any(isinstance(x,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='scope_payload' for t in x.targets) for x in n.body))
        for fail in [False,True]:
            cut=types.ModuleType('apply_cuts');cut.science='experimental authority'
            saved=dict(vars(cut))
            def load(*args,**kw):
                cut.science='temporary SIMC state';cut.new_field='temporary'
                if fail:raise RuntimeError('source load failure')
                return {'diagnostics':{}}
            ns=dict(_simc_scope_cuts=cut,_simc_scope_cut_state=dict(vars(cut)),inpDict={},phi_setting='Center',ParticleType='kaon',frozen_t_bins=[.4,.9],frozen_phi_bins=[-180.,180.],hgcer_cutg=None,kaon_signal_shape_payload=None,kaon_sigma0_shape_payload=None,load_setting_pion_component_shapes=load,load_kaon_simc_signal_shape=load,load_kaon_simc_sigma0_shape=load)
            code=compile(ast.Module(body=[block],type_ignores=[]),'actual_scoped_loader_finally','exec')
            if fail:
                with self.assertRaisesRegex(RuntimeError,'source load failure'):exec(code,ns)
            else:exec(code,ns)
            self.assertEqual(dict(vars(cut)),saved)

    def test_reuse_and_production_edge_handoff_are_population_bound(self):
        source=(REPO/'src/cuts/pion_t_bin_parents.py').read_text();self.assertIn('parent.get("template_population_identity") == template_population_identity',source)
        source=(REPO/'src/cuts/rand_sub.py').read_text();self.assertIn('rand_sub_direct_event_signal_scopes',source);self.assertIn('rand_sub_direct_event_sigma_scopes',source)
        self.assertNotIn('TH1D_to_TH2D(', (REPO/'src/plotting/data_vs_simc.py').read_text())

if __name__=='__main__':unittest.main()
