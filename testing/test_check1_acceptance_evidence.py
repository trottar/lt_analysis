"""Pure deterministic Check 1 evidence tests; no ROOT or farm execution."""
from copy import deepcopy
from contextlib import contextmanager, redirect_stderr
import ast
import io
import json
import math
import os
import shutil
from pathlib import Path
import sys
import tempfile
import time
import types
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"src/utility"),str(ROOT/"src/cuts"),str(ROOT/"src/binning")]
import check1_acceptance_evidence as ev
import pion_t_phi_template_comparison as observer
from pion_background_chain_audit import snapshot
from testing.test_pion_t_phi_template_comparison import H
from testing import run_pion_t_template_comparison_gate as owner


def inp(eps="lowe"):
    return dict(EPSSET=eps,ParticleType="kaon",Q2="4p4",W="2p74",tmin=.4,tmax=.9,mm_min=1.10,mm_max=1.16)


def meta(kind="data",setting="Center",eps="lowe"):
    populations = ("simc",) if kind == "simc" else ("kaon","pion")
    bindings = {p+":"+s:dict(tree="h10" if p == "simc" else p+"_"+s,
                original=dict(sha256="a"*64,bytes=12),runtime=dict(sha256="b"*64,bytes=13))
                for p in populations for s in (("simc",) if p == "simc" else ev.SOURCES)}
    if kind=="simc":
        roles=ev.simc_derivation(setting,eps);original=dict(sha256="a"*64,bytes=12)
        authority=dict(original_inputs={roles[k]:original for k in ("archive","archive_hist")},
                       runtime_before={roles[k]:original for k in ("archive","archive_hist")},
                       simc_derivations={roles["model"]:roles})
        bindings=ev.source_bindings(authority,setting,eps,"simc",{"simc":"h10"},
            {roles["model"]:dict(sha256="b"*64,bytes=13),roles["model_hist"]:original})
    return ev.metadata(inp(eps),setting,kind,
        dict(t=np.linspace(.4,.9,101),mm=[.7,1.1,1.16,1.3],native_mm=[1.1,1.13,1.16],
             canonical_t=[.4,.5666666666666667,.7333333333333334,.9],phi=list(range(0,361,40))),bindings,
        source_sha="c"*40,attempt="controlled-attempt",executables={p:dict(sha256="d"*64,bytes=42) for p in ev.EXECUTABLES})


def row(pop="kaon",source=0,entry=0,weight=2.,t=.5,mm=1.12,phi=20.,allcuts=True):
    m=meta("simc" if pop == "simc" else "data")
    r=dict(source=source,entry=entry,tree="h10" if pop == "simc" else pop+"_"+ev.SOURCES[source],
           t=t,mm=mm,phi=phi,ti=ev.member(t,m["axes"]["canonical_t"],pop=="pion"),
           pi=ev.member(phi,m["axes"]["phi"],pop=="pion"),allcuts=allcuts,nommcuts=True,
           coefficient=weight,weight=weight,sumw2=weight*weight)
    if pop == "kaon": r.update(proton_factor=1.,rf_accept=True)
    elif pop == "pion": r.update(transfer=1.,transfer_bin=ev.bisect_right(m["axes"]["mm"],mm),transfer_fingerprint="e"*64)
    else: r.update(iter_weight=weight/3.,normfac=3.,coefficient=3.)
    return r


def packet():
    m=meta()
    tables={"kaon":ev.columns([row(weight=2.),row(source=1,weight=-2.),row(entry=1,weight=0.)],"kaon"),
            "pion":ev.columns([row("pion",weight=1.)],"pion")}
    return ev.build_packet(m,tables)


def five_fixture(directory):
    """Controlled rows through real projections, exporters and owner validators."""
    from testing import test_run_pion_t_template_comparison_gate as old
    from pion_background_chain_audit import digest as old_digest
    from data_coordinates import build_kaon_data_coordinate_contract
    out=Path(directory)/"artifacts";out.mkdir()
    legacy={p:dict(sha256=ev.file_hash(ROOT/p),bytes=(ROOT/p).stat().st_size) for p in owner.child_source_inventory(ROOT)}
    identities=dict(legacy,**{p:dict(sha256=ev.file_hash(ROOT/p),bytes=(ROOT/p).stat().st_size) for p in ev.EXECUTABLES})
    originals={};copies={};runtime={};ref=dict(children={},parents={},summaries={},simc={})
    for eps in ("lowe","highe"):
        ref["summaries"][eps]=dict(original_inputs=originals)
        ref["simc"][eps]=dict(mm_edges=np.linspace(1.1,1.16,101),t_vs_tmin_y_edges=np.linspace(.4,.9,101))
    for setting,eps in ev.SETTINGS:
        for name in [f'{setting}_{p}_Analysed_{r}_Q4p4W2p74_{eps}.root' for p in ('kaon','pion') for r in ('Data','Dummy')]+[
                f'Prod_Coin_Q4p4W2p74{setting.lower()}_{eps}'+ext for ext in ('.root','.hist')]:
            path=out/name;path.write_bytes(('controlled input '+name).encode())
            identity=owner.file_identity(path);originals[name]=dict(original=identity)
            copies[name]=dict(original=identity,runtime_container_before=identity,source=str(path));runtime[name]=identity
    # Execute the real producer's copy/name/iteration-call function without
    # importing main or ROOT. The controlled iteration callback only marks the
    # model bytes, proving archive/model roles are distinct; no physics is run.
    main=ast.parse((ROOT/'src/main.py').read_bytes().decode())
    prepare=next(n for n in main.body if isinstance(n,ast.FunctionDef) and n.name=='_prepare_iter_weighted_kaon_simc')
    def controlled_iteration(parameters,model,inputs,setting):
        path=Path(model)
        if path.name!=ev.simc_derivation(setting,inputs['EPSSET']+'e')['model']:
            raise AssertionError('producer passed an unexpected iteration target')
        path.write_bytes(path.read_bytes()+b' controlled iteration-zero model branch')
    for setting,eps in ev.SETTINGS:
        namespace=dict(Q2='4p4',W='2p74',EPSSET=eps[:-1],OUTPATH=str(out),os=os,shutil=shutil,
            initial_param_file='controlled-parameters',inpDict={'EPSSET':eps[:-1]},output_file_lst=[],
            iter_weight=controlled_iteration)
        exec(compile(ast.Module(body=[prepare],type_ignores=[]),'src/main.py','exec'),namespace)
        namespace['_prepare_iter_weighted_kaon_simc'](setting)
    runtime=owner.check1_runtime_models(ev,out,runtime)
    broad=np.linspace(.7,1.3,101).tolist();te=np.linspace(.4,.9,101).tolist()
    canonical=[.4,.5666666666666667,.7333333333333334,.9];phi=list(range(-180,181,40));mm=1.115
    b=ev.bisect_right(broad,mm)
    def expand(value):
        if isinstance(value,dict):
            result={k:expand(v) for k,v in value.items()}
            if result.get('edges')==[0.,1.]:result['edges']=broad
            for k in ('values','variances','stored_sumw2'):
                if isinstance(result.get(k),list) and len(result[k])==3:
                    result[k]=[result[k][0]]+[result[k][1] if i==b else 0. for i in range(1,101)]+[result[k][-1]]
            return result
        if isinstance(value,list):return [expand(v) for v in value]
        return value
    def pdf_bytes():
        objects=[b'<< /Type /Catalog /Pages 2 0 R >>',
                 b'<< /Type /Pages /Count 6 /Kids [3 0 R 4 0 R 5 0 R 6 0 R 7 0 R 8 0 R] >>']
        objects += [b'<< /Type /Page /Parent 2 0 R /MediaBox [0 0 72 72] /Resources << >> >>']*6
        data=b'%PDF-1.4\n';offsets=[0]
        for i,obj in enumerate(objects,1):
            offsets.append(len(data));data+=str(i).encode()+b' 0 obj\n'+obj+b'\nendobj\n'
        start=len(data);data+=b'xref\n0 9\n0000000000 65535 f \n'
        data+=b''.join(('%010d 00000 n \n'%offset).encode() for offset in offsets[1:])
        return data+b'trailer\n<< /Size 9 /Root 1 0 R >>\nstartxref\n'+str(start).encode()+b'\n%%EOF\n'
    rows={};histograms={};histlist={}
    for setting,eps in ev.SETTINGS:
        key=setting+':'+eps;stem=f'{setting}_kaon_Q4p4W2p74_{eps}_'
        pdf=out/(stem+'pion-t-template-comparison.pdf');pdf.write_bytes(pdf_bytes())
        with patch.object(old,'IDENTITIES',legacy):ab=expand(old.payload(setting,pdf))
        ab['canonical_t_edges']=canonical;ab['canonical_phi_edges']=phi;ab['setting']['EPSSET']=eps[:-1]
        ab['original_analysis_input_identity']={r:dict(status='available',path=str(out/f'{setting}_kaon_Analysed_{r.title()}_Q4p4W2p74_{eps}.root'),
             **{k:runtime[f'{setting}_kaon_Analysed_{r.title()}_Q4p4W2p74_{eps}.root'][k] for k in ('sha256','bytes')}) for r in ('data','dummy')}
        if eps=='highe':ab['capture_stage']='final_selected_rerun'
        for ti,p in enumerate(ab['parents']):
            p['t_edges']=canonical[ti:ti+2]
            for c in owner.COMPONENTS:p['templates'][c]['t_edges']=p['t_edges']
        with patch.object(old,'IDENTITIES',legacy):child=expand(old.child_payload(setting,ab))
        child['t_edges']=canonical;child['phi_edges']=phi
        if eps=='highe':child['setting']=deepcopy(ab['setting'])
        child['coordinate']=build_kaon_data_coordinate_contract(setting,{'shift':0.},{'shift':0.})
        for p in ab['parents']:p['coordinate_fingerprint']=child['coordinate']['coordinate_fingerprint']
        child['selection_provenance']=dict(physical_pion_control_mask_fingerprint='f'*64)
        child['source_coefficients']=dict(zip(ev.SOURCES,(1.,-1.,-1.,1.)))
        rr=[dict(source='prompt',entry=ti*9+pi,mm=mm,ti=ti,pi=pi,allcuts=True,nommcuts=True,weight=1.)
            for ti in range(3) for pi in range(9)]
        child['observations']={p:{cut:observer.project([dict(r,selected=r[cut]) for r in rr],broad,3,9) for cut in ev.CUTS} for p in ('kaon','pion')}
        child['selected_row_fingerprints']={p:old_digest(rr) for p in ('kaon','pion')}
        child['parents']=deepcopy(ab['parents'])
        for c in child['cells']:
            c['t_edges']=canonical[c['t_index']:c['t_index']+2];c['phi_edges']=phi[c['phi_index']:c['phi_index']+2]
            c['applications']['Local_t']=dict(status='available',weights=[0.]+[1.]*100+[0.],
                cuts={cut:{k:child['observations']['pion'][cut]['children'][c['t_index']][c['phi_index']][k] for k in ('values','variances')} for cut in ev.CUTS})
        ref['children'][key]=child;ref['parents'][key]=ab
        (out/(stem+'pion-t-template-comparison.json')).write_bytes(ev.canonical_json(ab))
        (out/(stem+'pion-t-phi-template-comparison.json')).write_bytes(ev.canonical_json(child))
        pc=dict(records=[],coordinate_fingerprint=child['coordinate']['coordinate_fingerprint'],source_accounting={})
        sections={s:{k:[] for k in ('adj_MM','adj_t','t_index','phi_index','phi_shift_deg','allcuts','nommcuts',
            'source_label','source_entry_index','proton_cleaning_factor','proton_rf_accept')} for s in ev.SOURCES}
        for s in ev.SOURCES:pc['source_accounting'][s]=dict(tree_name='pion_'+s)
        for r in rr:
            ti,pi=r['ti'],r['pi'];tv=(canonical[ti]+canonical[ti+1])/2;pv=(phi[pi]+phi[pi+1])/2
            pc['records'].append(dict(source_label='prompt',entry_index=r['entry'],adj_t=tv,phi_degrees=pv,
                source_tree_name='pion_prompt',coordinate_fingerprint=pc['coordinate_fingerprint'],adj_MM=mm,
                t_index=ti,phi_index=pi,allcuts=True,nommcuts=True,coefficient=1.))
            for k,v in dict(adj_MM=mm,adj_t=tv,t_index=ti,phi_index=pi,phi_shift_deg=pv,allcuts=True,nommcuts=True,
                source_label='prompt',source_entry_index=r['entry'],proton_cleaning_factor=1.,proton_rf_accept=True).items():sections['prompt'][k].append(v)
        histlist[key]=dict(phi_setting=setting,normfac_data=1.,normfac_dummy=1.,nWindows=1.,
            _authoritative_pion_control_source_cache=pc,_yield_data_event_cache=sections,
            H_t_DATA=H(values=[0.]*102,errors=[0.]*102,edges=te))
        rows[key]={p:deepcopy(rr) for p in ('kaon','pion')}
        hists={k:[[H(values=[0.]*102,errors=[0.]*102,edges=te if k=='t' else ref['simc'][eps]['mm_edges']) for _ in range(9)] for _ in range(3)] for k in ('t','mm')}
        sink=ev.SimcSink();native=np.zeros((3,9,100));fine=np.zeros((3,9,100,100))
        for r in rr:
            ti,pi=r['ti'],r['pi'];tv=(canonical[ti]+canonical[ti+1])/2;pv=(phi[pi]+phi[pi+1])/2
            sink.observe(r['entry'],mm,tv,pv,2.,True,True,ti,pi)
            mb=ev.bisect_right(ref['simc'][eps]['mm_edges'],mm);tb=ev.bisect_right(te,tv)
            for k,index in (('t',tb),('mm',mb)):hists[k][ti][pi].values[index]=6.;hists[k][ti][pi].errors[index]=6.
            native[ti,pi,mb-1]=6.;fine[ti,pi,0,tb-1]=6.
        ref['simc'][eps]['simc_mm_'+setting.lower()+'_values']=native
        ref['simc'][eps]['simc_mm_'+setting.lower()+'_errors']=native.copy()
        ref['simc'][eps]['simc_t_vs_tmin_'+setting.lower()+'_values']=fine
        histograms[key]=(sink,hists)
    auth=owner.check1_authority(ev,ref,out,copies,identities,old.SHA,'controlled-five-setting-attempt')
    ap=Path(directory)/'authority.json';ap.write_bytes(ev.canonical_json(auth))
    # Data's tree names are read from its actual producer decision; use the
    # existing prompt-tree resolver rather than redefining RF ownership.
    from prompt_trees import get_prompt_tree_name,get_rand_tree_name
    with patch.dict(os.environ,{ev.AUTHORITY_ENV:str(ap),ev.AUTHORITY_HASH_ENV:ev.file_hash(ap)}):
        for setting,eps in ev.SETTINGS:
            key=setting+':'+eps;hist=histlist[key]
            ev.export_data(hist,inp(eps),ref['children'][key],rows[key],out)
            sink,hists=histograms[key]
            ev.export_simc(sink,types.SimpleNamespace(GetName=lambda s=setting,e=eps:str(out/f'{s}_kaon_Simc_Q4p4W2p74_{e}.root')),
                3.,inp(eps),setting,canonical,phi,hists,out)
    return out,identities,runtime,ref,auth,histograms


class OptimizerReferencePublicationTests(unittest.TestCase):
    """Real caller/yield/export sequence with controlled events and ROOT histograms.

    Numerical loop, normalization, native closure, authority and packet writer
    execute unchanged; the small histogram adapter replaces ROOT, not validators.
    """
    def sequence(self, eps, enabled=True, ordinary=False):
        from collections import defaultdict
        from binning_helpers import find_2d_bin_indices
        from background_config import (BG_OPT_PREPASS_CONTEXT,
            resolve_particle_subtraction_mode, resolve_pion_subtraction_scope)

        class Hist(H):
            def __init__(self, name, title, n, lo, hi, *second):
                edges=np.linspace(lo,hi,n+1).tolist()
                shape=(n+2,) if not second else (n+2,second[0]+2)
                self.values=np.zeros(shape);self.errors=np.zeros(shape);self.edges=edges
                self.second=np.linspace(second[1],second[2],second[0]+1) if second else None
            def Fill(self, *args):
                indices=(ev.bisect_right(self.edges,args[0]),)
                if self.second is not None:indices+=(ev.bisect_right(self.second,args[1]),)
                weight=args[len(indices)] if len(args)>len(indices) else 1.
                self.values[indices]+=weight
                self.errors[indices]=math.hypot(self.errors[indices],weight)
            def Scale(self, scale):
                self.values*=scale;self.errors*=abs(scale)
            def GetBinCenter(self, index):return (self.edges[index-1]+self.edges[index])/2

        names={'find_yield_simc','calculate_yield_simc','bin_simc','process_hist_simc',
            '_process_yield_simc_tree','_init_hist_group_matrices','_init_ave_simc_event_cache',
            '_append_ave_simc_event','_freeze_ave_simc_event_cache','_noop_progress_bar',
            '_assert_canonical_t_edges','integral_with_stat_error'}
        module=types.ModuleType('calculate_yield')
        ns=module.__dict__
        def unavailable(*args):raise ValueError('controlled ROOT-independent provenance unavailable')
        ns.update(np=np,math=math,os=os,sys=sys,defaultdict=defaultdict,TH1D=Hist,TH2D=Hist,
            BG_OPT_PREPASS_CONTEXT=BG_OPT_PREPASS_CONTEXT,
            resolve_particle_subtraction_mode=resolve_particle_subtraction_mode,
            resolve_pion_subtraction_scope=resolve_pion_subtraction_scope,
            find_2d_bin_indices=find_2d_bin_indices,
            calculate_theta_cm_deg=lambda *a:45.,_get_simc_true_theta_cm_deg=lambda *a:45.,
            _fill_t_vs_tmin_hist=lambda *a,**k:None,
            _clone_hist_for_plot=lambda h,**k:deepcopy(h),_build_f6_3_simc_provenance=unavailable)
        tree=ast.parse((ROOT/'src/binning/calculate_yield.py').read_text())
        exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)
            and n.name in names],type_ignores=[]),'actual_calculate_yield','exec'),ns)
        optimizer=dict(np=np,sys=sys)
        tree=ast.parse((ROOT/'src/utility/bg_optimization.py').read_text())
        exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)
            and n.name in ('_build_simc_reference','_prepare_candidate_hist')],type_ignores=[]),
            'actual_bg_optimization','exec'),optimizer)
        cuts=types.ModuleType('apply_cuts');cuts.set_val=lambda d:None
        cuts.apply_simc_cuts=lambda e,lo,hi:lo<=e.missmass<hi
        cuts.apply_simc_sub_cuts=lambda e,lo,hi:True
        hole=types.ModuleType('hgcer_hole')
        hole.apply_HGCer_hole_cut=lambda *a:types.SimpleNamespace(IsInside=lambda x,y:False)
        class Tree(list):
            def GetEntries(self):return len(self)
        t=np.linspace(.4,.9,4);phi=np.linspace(-180.,180.,10)
        events=Tree([types.SimpleNamespace(t=-(t[j]+t[j+1])/2,
            phipq=math.radians((phi[k]+phi[k+1])/2),missmass=1.115,iter_weight=2.,
            W=2.74,Q2=4.4,epsilon=.5,ssxptar=0.,ssyptar=0.,hsxptar=0.,hsyptar=0.,
            phgcer_x_det=0.,phgcer_y_det=0.) for j in range(3) for k in range(9)])
        with tempfile.TemporaryDirectory() as temp:
            out=Path(temp);ns['OUTPATH']=str(out)
            identities={p:dict(sha256=ev.file_hash(ROOT/p),bytes=(ROOT/p).stat().st_size) for p in ev.EXECUTABLES}
            roles=ev.simc_derivation('Center',eps);originals={}
            for role in ('archive','archive_hist','model','model_hist'):
                path=out/roles[role]
                path.write_bytes(b'controlled model' if role=='model' else b'controlled archive')
                if role in ('archive','archive_hist'):
                    originals[path.name]=dict(sha256=ev.file_hash(path),bytes=path.stat().st_size)
            auth=dict(schema_version='check1_owner_authority/v1',reference_sha256=ev.REFERENCE_SHA256,
                settings=[list(s) for s in ev.SETTINGS],private_output=str(out.resolve()),
                executables=identities,original_inputs=originals,runtime_before=originals,
                simc_derivations={roles['model']:roles},source_sha='a'*40,attempt='controlled-sequence',
                broad_mm={eps:np.linspace(.7,1.3,101).tolist()},reference_source_sha={'Center:'+eps:'b'*40})
            ap=out/'authority.json';ap.write_bytes(ev.canonical_json(auth))
            inputs=inp(eps);inputs.update(OutFilename='controlled',POL='',Q2min=4.,Q2max=5.,
                Wmin=2.,Wmax=3.,yield_emit_plots=False,yield_show_progress=False,
                particle_subtraction_mode='simc_shape_components',pion_subtraction_scope='t_bin')
            if ordinary:inputs['pion_subtraction_scope']='setting_wide'
            before=deepcopy(inputs)
            hist=dict(InFile_SIMC=types.SimpleNamespace(Get=lambda n:events,
                GetName=lambda:str(out/roles['model'])),normfac_simc=3.,phi_setting='Center',t_bins=t,phi_bins=phi)
            env={'LT_ANALYSIS_CHECK1_EVIDENCE':'1' if enabled else '0',
                ev.AUTHORITY_ENV:str(ap),ev.AUTHORITY_HASH_ENV:ev.file_hash(ap)}
            with patch.dict(sys.modules,{'calculate_yield':module,'apply_cuts':cuts,'hgcer_hole':hole}),patch.dict(os.environ,env):
                for _ in range(2):
                    reference,support=optimizer['_build_simc_reference'](hist,inputs,t,phi)
                    self.assertEqual(list(out.glob('*check1*')),[], 'Step 4 must never publish')
                    self.assertEqual(inputs,before)
                final=ns['find_yield_simc']([hist],inputs)
                self.assertEqual(reference['binned_SIMC']['Center']['yield'],final['binned_SIMC']['Center']['yield'])
                self.assertTrue(all(v=={'yield':6.,'yield_err':6.} for v in final['binned_SIMC']['Center']['yield'].values()))
                for name in ('mm','t','Q2','W'):
                    for j in range(3):
                        for k in range(9):
                            np.testing.assert_array_equal(support[name][j][k].values,hist['_xsect_support_simc'][name][j][k].values)
                            np.testing.assert_array_equal(support[name][j][k].errors,hist['_xsect_support_simc'][name][j][k].errors)
                packets=list(out.glob('*check1*.json'))
                if enabled and not ordinary:
                    self.assertEqual(len(packets),1)
                    doc,arrays=ev.read_packet(packets[0]);ev.validate_packet(doc,arrays)
                    saved={p.name:ev.file_hash(p) for p in out.glob('*check1*')}
                    self.assertEqual(len(saved),2)
                    with self.assertRaisesRegex(ValueError,'check1:destination_exists'):
                        ns['find_yield_simc']([hist],inputs)
                    self.assertEqual(saved,{p.name:ev.file_hash(p) for p in out.glob('*check1*')})
                else:
                    self.assertEqual(packets,[])
                    ns['find_yield_simc']([hist],inputs)
                    self.assertEqual(list(out.glob('*check1*')),[])

    def test_step4_then_step6_low_epsilon(self):self.sequence('lowe')
    def test_step4_then_step6_high_epsilon(self):self.sequence('highe')
    def test_default_off_both_epsilon_modes(self):
        for eps in ('lowe','highe'):
            with self.subTest(eps=eps):self.sequence(eps,enabled=False)
    def test_ordinary_scope_both_epsilon_modes(self):
        for eps in ('lowe','highe'):
            with self.subTest(eps=eps):self.sequence(eps,ordinary=True)


    def test_adjacent_repeated_call_guards(self):
        from background_config import (BG_OPT_PREPASS_CONTEXT,build_bg_optimization_prepass_config,
            resolve_particle_subtraction_mode,resolve_pion_subtraction_scope)
        tree=ast.parse((ROOT/'src/binning/calculate_yield.py').read_text())
        data=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='calculate_yield_data')
        simc=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='process_hist_simc')
        data_guard=next(n.test for n in ast.walk(data) if isinstance(n,ast.If)
            and any(isinstance(i,ast.ImportFrom) and i.module=='pion_t_phi_template_comparison' for i in n.body))
        simc_guard=next(n.test for n in ast.walk(simc) if isinstance(n,ast.If)
            and any(isinstance(i,ast.ImportFrom) and i.module=='check1_acceptance_evidence' for i in n.body))
        def allowed(guard,inputs,iteration=0):
            return eval(compile(ast.Expression(guard),'actual_export_guard','eval'),dict(os=os,
                inpDict=inputs,iteration=iteration,phi_setting='Center',
                BG_OPT_PREPASS_CONTEXT=BG_OPT_PREPASS_CONTEXT,
                resolve_particle_subtraction_mode=resolve_particle_subtraction_mode,
                resolve_pion_subtraction_scope=resolve_pion_subtraction_scope))
        with patch.dict(os.environ,{'LT_ANALYSIS_CHECK1_EVIDENCE':'1',
                                   'LT_ANALYSIS_PION_T_PHI_TEMPLATE_COMPARISON':'1'}):
            for eps in ('lowe','highe'):
                inputs=dict(inp(eps),particle_subtraction_mode='simc_shape_components',pion_subtraction_scope='t_bin')
                self.assertTrue(allowed(data_guard,inputs));self.assertTrue(allowed(simc_guard,inputs))
                prepass=build_bg_optimization_prepass_config(inputs)
                for _ in range(2):self.assertFalse(allowed(data_guard,prepass))
                self.assertFalse(allowed(simc_guard,inputs,iteration=True))
                reference=dict(inputs,check1_export_context='optimizer_simc_reference')
                self.assertFalse(allowed(simc_guard,reference))
                # This context suppresses only the Check 1 SIMC publisher.
                self.assertTrue(allowed(data_guard,reference))


class ScientificFailureDiagnosticsTests(unittest.TestCase):
    invariant = ('check1:scientific_identity:Left:lowe.cells[0].fit.'
                 'pion_component_alignment.applied_components.pi_delta.parent_alignment_hash')

    @staticmethod
    def documents():
        from data_coordinates import build_kaon_data_coordinate_contract
        parent = dict(parent_setting_key='Q4p4W2p74:Left:lowe', common_setting_shift_gev=.0012744878163235551,
            resolved_configuration_hash='c'*64, accepted=True, source='bin_local_scan', analysis_scope='t_bin',
            persistence_status='persisted', persistence_rejection_reasons=[], components={name:dict(
                residual_shift_gev=.001, window_expansion_gev=.002, effective_windows=[[1.1,1.16]],
                immutable_source_hist_identifier={'component':name,'label':'\u03bb'}, immutable_source_hist_checksum='d'*64)
                for name in ('pi_n','pi_delta','pi_sidis')})
        parent['parent_alignment_hash']=ev.digest(ev.alignment_parent_preimage(parent))
        parent['applied_components']=deepcopy(parent['components'])
        child=deepcopy(parent)
        child.update(accepted=False,source='current_common_shift_fallback',analysis_scope='t_phi_bin',common_setting_shift_gev=0.)
        for value in child['applied_components'].values():
            value.update(parent_alignment_hash=parent['parent_alignment_hash'],source='parent_setting_fallback',
                component_fallback_reason='parent alignment unavailable, stale, incompatible, or invalid')
        child = dict(applied_components=child['applied_components'],
                     **{k:v for k,v in child.items() if k != 'applied_components'})
        # Keep first traversal at the observed child hash, before embedded parents.
        document=dict(cells=[dict(t_index=0,phi_index=0,fit=dict(pion_component_alignment=child))],
                      parents=[dict(baseline_A=dict(pion_component_alignment=parent))],
                      coordinate=build_kaon_data_coordinate_contract('Left',{'shift':0.},{'shift':0.}))
        return deepcopy(document),document

    @staticmethod
    def change_parent(doc, field='immutable_source_hist_checksum', value='e'*64):
        parent=doc['parents'][0]['baseline_A']['pion_component_alignment']
        parent['components']['pi_delta'][field]=value
        parent['parent_alignment_hash']=ev.digest(ev.alignment_parent_preimage(parent))
        child=doc['cells'][0]['fit']['pion_component_alignment']
        child['parent_alignment_hash']=parent['parent_alignment_hash']
        child['applied_components']['pi_delta']['parent_alignment_hash']=parent['parent_alignment_hash']

    def receipt(self, actual, expected):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp); private=root/'private';private.mkdir()
            path=private/'Left_kaon_Q4p4W2p74_lowe_pion-t-phi-template-comparison.json'
            path.write_bytes(ev.canonical_json(actual))
            status=owner.GateStatus(root,root/'unpublished.zip','a'*40,high=True,check1=True)
            status.update('verify_artifacts',analysis_completed=True,analysis_returncode=0)
            stream=io.StringIO()
            with redirect_stderr(stream), self.assertRaises(ValueError) as caught:
                owner._verify_check1_packets(ev,private,'a'*40,{p:{} for p in owner.CHECK1_EXECUTABLES},1,
                    dict.fromkeys(owner.names(check1=True)),{},dict(children={'Left:lowe':expected}),{},
                    dict.fromkeys(owner.names()|owner.names(True)),{},diagnostic_status=status)
            # The real owner's JSON receipt survives removal of its private inputs.
            shutil.rmtree(private)
            raw=status.path.read_bytes().decode('utf-8')
            marker='"check1_scientific_failure": '
            start=raw.index(marker)+len(marker)
            diagnostic,end=json.JSONDecoder().raw_decode(raw,start)
            # Measure the actual nested value from the real writer, preserving
            # CRLF on Windows rather than normalizing via read_text().
            persisted_bytes=len(raw[start:end].encode('utf-8'))
            self.assertLessEqual(persisted_bytes,64*1024)
            self.assertGreaterEqual(owner._check1_diagnostic_persisted_size(diagnostic),persisted_bytes)
            saved=json.loads(raw)
            self.assertEqual(saved['first_failed_invariant'],str(caught.exception))
            self.assertEqual(saved['failure_reason'],str(caught.exception))
            self.assertEqual(saved['failed_stage'],'verify_artifacts')
            self.assertFalse(saved['artifact_verification_completed'])
            self.assertFalse(saved['collection_completed']);self.assertFalse(saved['zip_verification_completed'])
            self.assertFalse(saved['evidence_admissible']);self.assertFalse(saved['scientifically_accepted'])
            self.assertFalse((root/'unpublished.zip').exists())
            return saved['check1_scientific_failure'],str(caught.exception),stream.getvalue()

    def test_review_oversized_fixture_is_bounded_in_real_persisted_receipt(self):
        actual,expected=self.documents()
        for document in (actual,expected):
            alignment=document['parents'][0]['baseline_A']['pion_component_alignment']
            for component in alignment['components'].values():
                for field in ('immutable_source_hist_identifier','effective_windows'):
                    component[field]=[[[0]*10 for _ in range(5)] for _ in range(4)]
        self.change_parent(actual)
        original=ev.scientific_failure_diagnostic(actual,expected,'a'*40,ev.REFERENCE_SHA256,'Left:lowe',self.invariant)
        self.assertLess(len(ev.canonical_json(original)),64*1024)
        self.assertGreater(owner._check1_diagnostic_persisted_size(original),64*1024)
        diag,error,stderr=self.receipt(actual,expected)
        self.assertEqual(error,self.invariant);self.assertIn('bounded_fallback',diag)
        self.assertEqual(diag['first_difference'],original['first_difference'])
        self.assertEqual(diag['first_parent_operand_difference'],original['first_parent_operand_difference'])
        for side in ('actual','expected'):
            self.assertIn(diag['first_difference'][side]['value'],stderr)
            summary=diag[side]['parent_hash_operands']['components']['pi_delta']['effective_windows']
            self.assertNotIn('value',summary)
            self.assertEqual(summary['sha256'],ev.digest([[[0]*10 for _ in range(5)] for _ in range(4)]))
            self.assertEqual(diag[side]['parent']['selection'],original[side]['parent']['selection'])
            self.assertEqual(diag[side]['parent_hash_consistency'],original[side]['parent_hash_consistency'])

    def test_oversized_ancillary_context_fallback_keeps_first_hash_and_operand(self):
        actual,expected=self.documents();self.change_parent(actual)
        build=ev.scientific_failure_diagnostic
        original=build(actual,expected,'a'*40,ev.REFERENCE_SHA256,'Left:lowe',self.invariant)
        def oversized(*args,**kwargs):
            diagnostic=build(*args,**kwargs)
            diagnostic['actual']['extra_context']=['ancillary']*20000
            return diagnostic
        with patch.object(ev,'scientific_failure_diagnostic',side_effect=oversized):
            diag,error,_=self.receipt(actual,expected)
        self.assertEqual(error,self.invariant)
        self.assertEqual(diag['first_difference'],original['first_difference'])
        self.assertEqual(diag['first_parent_operand_difference'],original['first_parent_operand_difference'])
        self.assertIn('ancillary context summarized',diag['bounded_fallback'])
        self.assertIn('sha256',diag['actual']);self.assertNotIn('diagnostic_error',diag)

    def test_short_operand_windows_survive_oversized_ancillary_context(self):
        actual,expected=self.documents()
        for doc in (actual,expected):
            parent=doc['parents'][0]['baseline_A']['pion_component_alignment']
            for component in parent['components'].values():
                for field in ('immutable_source_hist_identifier','effective_windows'):
                    component[field]=[[[0]*10 for _ in range(5)] for _ in range(4)]
            parent['components']['pi_delta']['effective_windows']=[[1.1,1.16]]
        self.change_parent(actual,'effective_windows',[[1.1,1.16],[1.16,1.18]])
        original=ev.scientific_failure_diagnostic(actual,expected,'a'*40,ev.REFERENCE_SHA256,'Left:lowe',self.invariant)
        self.assertGreater(owner._check1_diagnostic_persisted_size(original),65536)
        diag,error,_=self.receipt(actual,expected)
        self.assertEqual(error,self.invariant)
        self.assertEqual(diag['first_difference'],original['first_difference'])
        self.assertEqual(diag['first_parent_operand_difference'],original['first_parent_operand_difference'])
        difference=diag['first_parent_operand_difference']
        self.assertEqual(difference['path'],'parent_hash_operands.components.pi_delta.effective_windows')
        self.assertEqual(difference['actual']['value'],[[1.1,1.16],[1.16,1.18]])
        self.assertEqual(difference['expected']['value'],[[1.1,1.16]])

    def test_short_scientific_first_difference_survives_context_only_fallback(self):
        actual,expected=self.documents()
        actual['cells'][0]['fit']['pion_component_alignment']['applied_components']['pi_delta']['effective_windows']=[[1.1,1.16],[1.16,1.18]]
        build=ev.scientific_failure_diagnostic
        def oversized(*args,**kwargs):
            diagnostic=build(*args,**kwargs)
            diagnostic['actual']['extra_context']=['ancillary']*20000
            return diagnostic
        with patch.object(ev,'scientific_failure_diagnostic',side_effect=oversized):
            diag,error,_=self.receipt(actual,expected)
        self.assertEqual(error,'check1:scientific_length:Left:lowe.cells[0].fit.pion_component_alignment.applied_components.pi_delta.effective_windows')
        self.assertEqual(diag['first_difference']['actual']['value'],[[1.1,1.16],[1.16,1.18]])
        self.assertEqual(diag['first_difference']['expected']['value'],[[1.1,1.16]])
        self.assertIn('ancillary context summarized',diag['bounded_fallback'])

    def test_only_overbudget_first_operand_is_summarized_with_complete_digest(self):
        actual,expected=self.documents()
        windows=[1.1,1.16]
        for _ in range(210):windows=[windows]
        self.change_parent(actual,'effective_windows',windows)
        original=ev.scientific_failure_diagnostic(actual,expected,'a'*40,ev.REFERENCE_SHA256,'Left:lowe',self.invariant)
        operand=original['first_parent_operand_difference']
        self.assertIn('value',operand['actual'])
        self.assertGreater(owner._check1_diagnostic_persisted_size({'first_parent_operand_difference':operand}),65536)
        diag,error,_=self.receipt(actual,expected)
        self.assertEqual(error,self.invariant)
        self.assertEqual(diag['first_difference'],original['first_difference'])
        saved=diag['first_parent_operand_difference']
        self.assertNotIn('value',saved['actual'])
        self.assertEqual(saved['actual']['sha256'],ev.digest(operand['actual']['value']))
        self.assertEqual(saved['expected'],operand['expected'])
        self.assertEqual(saved['expected']['value'],[1.1,1.16])

    def test_owner_hash_failure_receipt_contains_exact_operands_and_selection(self):
        actual,expected=self.documents();self.change_parent(actual)
        # Other child scalars stay identical; inherited hash alone still fails.
        actual['cells'][0]['fit']['pion_component_alignment']['parent_alignment_hash']=expected['cells'][0]['fit']['pion_component_alignment']['parent_alignment_hash']
        diag,error,stderr=self.receipt(actual,expected)
        self.assertEqual(error,self.invariant)
        self.assertEqual(diag['first_difference']['path'],self.invariant.split(':',2)[2])
        self.assertEqual(diag['first_parent_operand_difference']['path'],
            'parent_hash_operands.components.pi_delta.immutable_source_hist_checksum')
        for side in ('actual','expected'):
            full_hash=diag['first_difference'][side]['value']
            self.assertEqual(len(full_hash),64);self.assertIn(full_hash,stderr)
            self.assertEqual(diag[side]['parent_hash_consistency']['status'],'matches')
            self.assertTrue(diag[side]['parent']['selection']['accepted']['value'])
            self.assertFalse(diag[side]['child']['selection']['accepted']['value'])
        self.assertEqual(diag['actual']['lineage_consistency'],'inconsistent')
        self.assertEqual(diag['expected']['lineage_consistency'],'consistent')
        self.assertEqual(diag['source_sha'],'a'*40);self.assertEqual(diag['reference_sha256'],ev.REFERENCE_SHA256)
        self.assertLess(len(ev.canonical_json(diag)),64*1024)

    def test_exact_shift_operand_even_below_scientific_tolerance(self):
        actual,expected=self.documents();self.change_parent(actual,'residual_shift_gev',.001+1e-14)
        diag,error,_=self.receipt(actual,expected)
        self.assertEqual(error,self.invariant)
        diff=diag['first_parent_operand_difference']
        self.assertEqual(diff['path'],'parent_hash_operands.components.pi_delta.residual_shift_gev')
        self.assertEqual(diff['actual']['value'],.001+1e-14)
        self.assertEqual(diff['expected']['value'],.001)

    def test_missing_parent_and_corrupt_stored_hash_are_explicit(self):
        for missing in (False,True):
            actual,expected=self.documents();self.change_parent(actual)
            if missing:actual['parents']=[]
            else:actual['parents'][0]['baseline_A']['pion_component_alignment']['parent_alignment_hash']='0'*64
            diag,error,_=self.receipt(actual,expected)
            self.assertEqual(error,self.invariant)
            self.assertEqual(diag['actual']['parent_hash_consistency']['status'],'not_verified' if missing else 'mismatch')
            self.assertEqual(diag['actual']['lineage_consistency'],'missing' if missing else 'inconsistent')

    def test_diagnostic_failure_never_replaces_original_exception(self):
        actual,expected=self.documents();self.change_parent(actual)
        with patch.object(ev,'scientific_failure_diagnostic',side_effect=RuntimeError('private content must not leak')):
            diag,error,stderr=self.receipt(actual,expected)
        self.assertEqual(error,self.invariant);self.assertEqual(diag['diagnostic_error'],'RuntimeError')
        self.assertNotIn('private content',json.dumps(diag)+stderr)
        original=ValueError(self.invariant)
        with patch.object(ev,'compare_science',side_effect=original),redirect_stderr(io.StringIO()):
            with self.assertRaises(ValueError) as caught:
                owner._compare_check1_science(ev,{}, {},'Left:lowe',{}, {},'a'*40,ev.REFERENCE_SHA256)
        self.assertIs(caught.exception,original)

    def test_unrelated_differences_preserve_comparator_semantics(self):
        cases=[({'x':1.},{'x':2.},'scientific_value:Left:lowe.x'),
               (None,'hash','scientific_identity:Left:lowe'),
               ({'x':1},{'y':1},'scientific_keys:Left:lowe'),
               ([1],[1,2],'scientific_length:Left:lowe'),
               ('Center','Left','scientific_identity:Left:lowe')]
        for actual,expected,error in cases:
            with self.subTest(error=error),redirect_stderr(io.StringIO()):
                with self.assertRaises(ValueError) as caught:
                    owner._compare_check1_science(ev,actual,expected,'Left:lowe',{}, {},'a'*40,ev.REFERENCE_SHA256)
                self.assertEqual(str(caught.exception),'check1:'+error)
                diag=ev.scientific_failure_diagnostic({}, {},'a'*40,ev.REFERENCE_SHA256,'Left:lowe',str(caught.exception),actual,expected)
                self.assertEqual(diag['first_difference']['category'],error.split(':')[0])

    def test_success_is_silent_and_leaves_owner_status_bytes_unchanged(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);status=owner.GateStatus(root,root/'none.zip','a'*40,high=True,check1=True)
            before=status.path.read_bytes();stream=io.StringIO()
            with patch.object(ev,'scientific_failure_diagnostic') as diagnostic,redirect_stderr(stream):
                owner._compare_check1_science(ev,{'x':1.+1e-10},{'x':1.},'Left:lowe',{}, {},'a'*40,ev.REFERENCE_SHA256,status)
            diagnostic.assert_not_called();self.assertEqual(stream.getvalue(),'');self.assertEqual(status.path.read_bytes(),before)

    def test_status_write_and_stderr_failure_do_not_mask_scientific_failure(self):
        actual,expected=self.documents();self.change_parent(actual)
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);status=owner.GateStatus(root,root/'none.zip','a'*40,high=True,check1=True)
            with patch.object(status,'update',side_effect=OSError('closed status')),patch.object(owner.sys,'stderr') as stream:
                stream.write.side_effect=OSError('closed stderr')
                with self.assertRaises(ValueError) as caught:
                    owner._compare_check1_science(ev,ev.scientific_projection(actual),ev.scientific_projection(expected),
                        'Left:lowe',actual,expected,'a'*40,ev.REFERENCE_SHA256,status)
            self.assertEqual(str(caught.exception),self.invariant)
            self.assertEqual(status.record['first_failed_invariant'],self.invariant)
            self.assertEqual(status.record['diagnostic_write_failure'],'OSError')

    def test_large_identifiers_keep_complete_digests_and_exact_first_scalar(self):
        actual,expected=self.documents()
        for doc in (actual,expected):
            p=doc['parents'][0]['baseline_A']['pion_component_alignment']
            p['components']['pi_delta']['immutable_source_hist_identifier']={'bins':list(range(10000))}
        actual['parents'][0]['baseline_A']['pion_component_alignment']['components']['pi_delta']['immutable_source_hist_identifier']['bins'][37]=-37
        actual['cells'][0]['fit']['pion_component_alignment']['applied_components']['pi_delta']['parent_alignment_hash']='e'*64
        diag,_,_=self.receipt(actual,expected)
        diff=diag['first_parent_operand_difference']
        self.assertEqual(diff['path'],'parent_hash_operands.components.pi_delta.immutable_source_hist_identifier.bins[37]')
        self.assertEqual(diff['actual']['value'],-37)
        summary=diag['actual']['parent_hash_operands']['components']['pi_delta']['immutable_source_hist_identifier']
        self.assertNotIn('value',summary)
        self.assertLess(len(ev.canonical_json(diag)),64*1024)

    def test_preimage_matches_unchanged_producer_and_large_values_are_bounded(self):
        # Execute the exact frozen producer function without importing ROOT.
        tree=ast.parse((ROOT/'src/cuts/pion_component_fits.py').read_bytes())
        function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='_alignment_parent_hash')
        scope=dict(_alignment_hash=ev.digest,COMPONENT_NAMES=('pi_n','pi_delta','pi_sidis'))
        exec(compile(ast.Module(body=[function],type_ignores=[]),'<frozen-parent-hash>','exec'),scope)
        actual,_=self.documents();alignment=actual['parents'][0]['baseline_A']['pion_component_alignment']
        self.assertEqual(ev.digest(ev.alignment_parent_preimage(alignment)),scope['_alignment_parent_hash'](
            alignment['parent_setting_key'],alignment['common_setting_shift_gev'],alignment['components'],alignment['resolved_configuration_hash']))
        alignment['components']['pi_delta']['immutable_source_hist_identifier']=['large']*10000
        summary=ev.diagnostic_value(alignment['components']['pi_delta']['immutable_source_hist_identifier'])
        self.assertNotIn('value',summary);self.assertEqual(summary['sha256'],ev.digest(['large']*10000))
        self.assertIsNone(ev.first_exact_difference({'x':1.},{'x':1.}))
        self.assertIsNotNone(ev.first_exact_difference({'x':-0.},{'x':0.}))


class PacketTests(unittest.TestCase):
    def test_signed_cancellation_counts_and_quadratic_sumw2(self):
        doc,a=packet()
        self.assertEqual(a["allcuts/before/sumw_joint"].sum(),0.)
        self.assertEqual(a["allcuts/before/count_joint"].sum(),3)
        self.assertEqual(a["allcuts/before/support_count_joint"].sum(),2)
        self.assertEqual(a["allcuts/before/sumw2_joint"].sum(),8.)
        self.assertEqual(a["allcuts/after/sumw_joint"].sum(),-1.)
        self.assertEqual(a["allcuts/after/sumw2_joint"].sum(),9.)
        self.assertEqual(a["allcuts/after/absw_joint"].sum(),5.)
        self.assertEqual(doc["slots"][1]["status"],"empty")
        ev.validate_packet(doc,a)

    def test_actual_roundtrip_and_collision(self):
        doc,a=packet()
        with tempfile.TemporaryDirectory() as temp:
            saved=ev.write_packet(temp,doc,a)
            path=Path(temp)/(ev.basename(doc)+".json")
            loaded,b=ev.read_packet(path)
            self.assertEqual(saved,loaded)
            self.assertEqual(list(b),sorted(b))
            self.assertEqual(b["rows/kaon/entry"].dtype,np.dtype("<i8"))
            self.assertEqual(loaded["units"]["t"],"GeV^2")
            with self.assertRaisesRegex(ValueError,"destination_exists"): ev.write_packet(temp,doc,a)
            target=path.with_suffix(".npz"); raw=target.read_bytes();target.write_bytes(raw+b"tampered")
            with self.assertRaisesRegex(ValueError,"npz_identity"): ev.read_packet(path)

    def test_pickled_arrays_rejected_even_with_matching_file_hash(self):
        d,a=packet()
        with tempfile.TemporaryDirectory() as temp:
            saved=ev.write_packet(temp,d,a)
            path=Path(temp)/(ev.basename(d)+".json");npz=path.with_suffix(".npz")
            a["rows/kaon/entry"]=a["rows/kaon/entry"].astype(object)
            with npz.open("wb") as stream:np.savez_compressed(stream,**a)
            saved["npz"].update(sha256=ev.file_hash(npz),bytes=npz.stat().st_size)
            path.write_bytes(ev.canonical_json(saved))
            with self.assertRaisesRegex(ValueError,"Object arrays"):ev.read_packet(path)

    def test_existing_dangling_sidecar_symlink_is_not_replaced(self):
        d,a=packet()
        with tempfile.TemporaryDirectory() as temp:
            destination=Path(temp)/(ev.basename(d)+'.json')
            try:destination.symlink_to(Path(temp)/'absent')
            except OSError as exc:self.skipTest('local symlink capability unavailable: '+str(exc))
            with self.assertRaisesRegex(ValueError,'destination_exists'):ev.write_packet(temp,d,a)
            self.assertTrue(destination.is_symlink())

    def test_combined_profile_has_exact_five_settings_and_35_unique_artifacts(self):
        from testing import collect_pion_hgcer_validation_bundle as collector
        profile=collector.load_validation_profile(ROOT/"testing/pion_check1_acceptance_evidence_profile_v1.json")
        declarations=profile["artifacts"]["global"]
        self.assertEqual(len(declarations),35)
        self.assertEqual(len({r["key"] for r in declarations}),35)
        self.assertEqual(len({r["basename_template"] for r in declarations}),35)
        self.assertEqual(set(collector.resolve_settings(None,None,profile)),set(ev.SETTINGS))
        self.assertEqual(profile["artifacts"]["settings"],[])
        for setting,eps in ev.SETTINGS:
            names=[r["basename_template"] for r in declarations if r["basename_template"].startswith(setting+"_")
                   and "_"+eps+"_" in r["basename_template"]]
            self.assertEqual(len(names),7)
            self.assertEqual(sum("_check1_" in n for n in names),4)
        self.assertTrue(all(r["required"] for r in declarations))

    def test_all_five_settings_and_default_off_wrong_modes(self):
        for setting,eps in ev.SETTINGS:
            self.assertEqual(meta(setting=setting,eps=eps)["setting"],setting)
            self.assertFalse(ev.requested(inp(eps),setting,environ={}))
            self.assertTrue(ev.requested(inp(eps),setting,environ={ev.ENABLE_ENV:"1"}))
            self.assertFalse(ev.requested(inp(eps),setting,1,{ev.ENABLE_ENV:"1"}))
            for k,v in (("bg_opt_prepass_mode",True),("ParticleType","pion"),("Q2","3p0"),
                        ("analysis_subtraction_context","bg_optimizer_pre_parent")):
                self.assertFalse(ev.requested(dict(inp(eps),**{k:v}),setting,environ={ev.ENABLE_ENV:"1"}))

    def test_exporter_has_no_fitter_root_traversal_or_histogram_mutators(self):
        source=ast.parse((ROOT/'src/utility/check1_acceptance_evidence.py').read_text())
        calls={n.func.id if isinstance(n.func,ast.Name) else n.func.attr
               for n in ast.walk(source) if isinstance(n,ast.Call) and isinstance(n.func,(ast.Name,ast.Attribute))}
        self.assertFalse(calls&{'set_val','Fill','Scale','SetBinContent','SetBinError','GetEntry','GetEntries',
            'build_particle_subtraction_component_result','get_kaon_proton_cleaning_event_payload','Popen','run'})

    def test_actual_simc_construction_guard_excludes_disabled_prepass_iteration_and_wrong_scope(self):
        source=ast.parse((ROOT/'src/binning/calculate_yield.py').read_text())
        function=next(n for n in source.body if isinstance(n,ast.FunctionDef) and n.name=='process_hist_simc')
        guarded=next(n for n in function.body if isinstance(n,ast.If) and any(isinstance(v,ast.Constant) and v.value==ev.ENABLE_ENV for v in ast.walk(n.test)))
        expression=compile(ast.Expression(guarded.test),'actual_simc_guard','eval')
        state=dict(os=types.SimpleNamespace(environ={ev.ENABLE_ENV:'1'}),inpDict=inp(),phi_setting='Center',iteration=0,
            BG_OPT_PREPASS_CONTEXT='bg_optimizer_pre_parent',resolve_particle_subtraction_mode=lambda d:'simc_shape_components',
            resolve_pion_subtraction_scope=lambda d:'t_bin')
        self.assertTrue(eval(expression,state))
        for changes in ({'iteration':1},{'os':types.SimpleNamespace(environ={})},{'inpDict':dict(inp(),bg_opt_prepass_mode=True)},
                        {'inpDict':dict(inp(),analysis_subtraction_context='bg_optimizer_pre_parent')},
                        {'inpDict':dict(inp(),ParticleType='pion')},{'phi_setting':'Right'},
                        {'resolve_pion_subtraction_scope':lambda d:'global'}):
            self.assertFalse(eval(expression,{**state,**changes}))
        self.assertTrue(any(isinstance(n,ast.ImportFrom) and n.module=='check1_acceptance_evidence' for n in guarded.body))

    def test_boundaries_and_nextafter_owner_conventions(self):
        m=meta(); edges=m["axes"]["canonical_t"]
        for i,edge in enumerate(edges):
            for value in (np.nextafter(edge,-math.inf),edge,np.nextafter(edge,math.inf)):
                canonical=ev.member(value,edges,True); halfopen=ev.member(value,edges)
                self.assertEqual(canonical,halfopen if value!=edges[-1] else 2)
                if 0<=halfopen<3:
                    t=ev.columns([row(t=float(value))],"kaon");ev.validate_rows(t,"kaon",m)
        self.assertEqual(ev.member(edges[-1],edges),-1)
        t=ev.columns([row("pion",t=.9,phi=360.,mm=1.3,allcuts=False)],"pion")
        ev.validate_rows(t,"pion",m)
        self.assertEqual(t["ti"][0],2);self.assertEqual(t["pi"][0],8)
        self.assertEqual(t["transfer_bin"][0],4)
        rows=[row(mm=.6,allcuts=False),row(entry=1,mm=1.3,allcuts=False),row(entry=2,phi=-1.)]
        d,a=ev.build_packet(m,{"kaon":ev.columns(rows,"kaon"),"pion":ev.columns([],"pion")})
        self.assertEqual(a["nommcuts/before/count_joint"][:,:,0].sum(),1)
        self.assertEqual(a["nommcuts/before/count_joint"][:,:,-1].sum(),1)
        self.assertEqual(d["rows"]["kaon"],3)

    def test_nonfinite_object_duplicate_wrong_axis_and_identity_rejected(self):
        for pop in ("kaon","pion","simc"):
            m=meta("simc" if pop=="simc" else "data")
            for key in ("mm","t","phi","weight","sumw2"):
                for value in (float("nan"),float("inf"),-float("inf")):
                    tab=ev.columns([row(pop)],pop);tab[key][0]=value
                    with self.assertRaisesRegex(ValueError,"nonfinite"):ev.validate_rows(tab,pop,m)
            tab=ev.columns([row(pop),row(pop)],pop)
            with self.assertRaisesRegex(ValueError,"identity"):ev.validate_rows(tab,pop,m)
            tab=ev.columns([row(pop)],pop);tab["entry"]=tab["entry"].astype(object)
            with self.assertRaisesRegex(ValueError,"dtype"):ev.validate_rows(tab,pop,m)
            tab=ev.columns([row(pop)],pop);tab["tree"][0]="wrong_tree"
            with self.assertRaisesRegex(ValueError,"source_binding"):ev.validate_rows(tab,pop,m)
        d,a=packet();key="allcuts/before/sumw_joint";a[key]=a[key].transpose(0,2,1)
        with self.assertRaisesRegex(ValueError,"array_type"):ev.validate_packet(d,a)
        with self.assertRaisesRegex(ValueError,"axis"):ev.axis([1.,.5],"bad")
        with self.assertRaisesRegex(ValueError,"accepted_interval"):
            m=meta();m["axes"]["canonical_t"][-1]=.8
            ev.metadata(inp(),"Center","data",m["axes"],m["source_bindings"],source_sha="a"*40,attempt="x",executables=m["executables"])

    def test_normalization_once_and_deliberate_double_norm_rejected(self):
        sink=ev.SimcSink();sink.observe(0,1.12,.5,20.,2.,True,True,0,0)
        table=sink.table(3.)
        self.assertEqual(table["weight"].tolist(),[6.])
        self.assertEqual(table["sumw2"].tolist(),[36.])
        ev.validate_rows(table,"simc",meta("simc"))
        table["weight"]*=3.;table["sumw2"]=table["weight"]**2
        with self.assertRaisesRegex(ValueError,"normalization_once"):ev.validate_rows(table,"simc",meta("simc"))

    def test_injected_write_and_link_failures_leave_no_complete_pair(self):
        d,a=packet()
        for target in ("np.savez_compressed","os.link"):
            with tempfile.TemporaryDirectory() as temp:
                module,attr=target.split(".")
                with patch.object(getattr(ev,module),attr,side_effect=OSError("injected")):
                    with self.assertRaisesRegex(OSError,"injected"):ev.write_packet(temp,d,a)
                self.assertEqual(list(Path(temp).iterdir()),[])
        with tempfile.TemporaryDirectory() as temp:
            actual=ev.os.link;calls=[]
            def fail_second(source,target):
                calls.append(target)
                if len(calls)==2: raise OSError("metadata publication")
                return actual(source,target)
            with patch.object(ev.os,"link",side_effect=fail_second):
                with self.assertRaisesRegex(OSError,"metadata publication"):ev.write_packet(temp,d,a)
            self.assertEqual(list(Path(temp).iterdir()),[])

    def test_metadata_and_joint_tampering_fail(self):
        for key,value in (("iteration",1),("branch","Global"),("reference_sha256","0"*64),
                          ("scientifically_accepted",True),("units",{})):
            d,a=packet();d[key]=value
            with self.assertRaises(ValueError):ev.validate_packet(d,a)
        d,a=packet();a["allcuts/after/sumw2_joint"][0,0,0]=1.
        with self.assertRaisesRegex(ValueError,"closure"):ev.validate_packet(d,a)
        d,a=packet();a["allcuts/after/sumw2_joint"][0,0,0]=-1e-15
        with self.assertRaisesRegex(ValueError,"negative_moment"):ev.validate_packet(d,a)

    def test_postpublication_context_cleanup_failure_rolls_back_only_own_pair(self):
        d,a=packet();temporary=ev.tempfile.TemporaryDirectory
        for mode in ('ordinary','replaced_metadata','rollback_unlink_error'):
            with self.subTest(mode=mode),temporary() as temp:
                root=Path(temp);unrelated=root/'unrelated';unrelated.write_bytes(b'preserve')
                failure=OSError('exact post-publication cleanup failure');seen=[]
                cleanup=temporary.cleanup
                def cleanup_failure(pending):
                    jp=root/(ev.basename(d)+'.json')
                    saved,_=ev.read_packet(jp)  # both links really exist and validate
                    seen.append(saved['status'])
                    if mode=='replaced_metadata':
                        # Retain another hardlink so a replacement cannot
                        # accidentally reuse the unlinked inode.
                        os.link(jp,root/'retained-old-json');jp.unlink();jp.write_bytes(b'competitor')
                    cleanup(pending)
                    raise failure
                unlink=Path.unlink
                def injected_unlink(path,*args,**kwargs):
                    if mode=='rollback_unlink_error' and path.name.endswith('_data.json'):
                        raise OSError('secondary rollback error')
                    return unlink(path,*args,**kwargs)
                with patch.object(temporary,'cleanup',cleanup_failure),patch.object(Path,'unlink',injected_unlink):
                    with self.assertRaises(OSError) as caught:ev.write_packet(root,d,a)
                self.assertIs(caught.exception,failure)
                self.assertEqual(seen,['captured_pending_owner_verification'])
                self.assertEqual(unrelated.read_bytes(),b'preserve')
                self.assertFalse((root/(ev.basename(d)+'.npz')).exists())
                jp=root/(ev.basename(d)+'.json')
                if mode=='ordinary':self.assertFalse(jp.exists())
                elif mode=='replaced_metadata':self.assertEqual(jp.read_bytes(),b'competitor')
                else:self.assertTrue(failure.check1_rollback_failures)

    def test_publication_failure_precedes_context_cleanup_failure(self):
        d,a=packet();temporary=ev.tempfile.TemporaryDirectory
        with temporary() as temp:
            failure=OSError('first publication failure');cleanup=temporary.cleanup;link=ev.os.link
            def fail_link(source,target):
                if str(target).endswith('.json'):raise failure
                return link(source,target)
            def fail_cleanup(pending):
                cleanup(pending)
                raise OSError('later context cleanup failure')
            with patch.object(temporary,'cleanup',fail_cleanup),patch.object(ev.os,'link',fail_link):
                with self.assertRaises(OSError) as caught:ev.write_packet(temp,d,a)
            self.assertIs(caught.exception,failure)
            self.assertEqual(list(Path(temp).iterdir()),[])

    def test_real_projection_transfer_and_serializer_closure(self):
        edges=meta()["axes"]["mm"]
        h=H(values=[0.,1.,1.,1.,0.],errors=[0.,1.,1.,1.,0.],edges=edges)
        fit=dict(particle_subtraction_mode="simc_shape_components",fit_status_pion="success",
                 fit_status_kaon="success",fallback_used=False,
                 diagnostics={p:dict(validation=dict(accepted=True)) for p in ("pion","kaon")},
                 **{key:float(value) for key,value in (("A_n",2),("A_delta",2),("A_sidis",2),
                                                     ("B_n",1),("B_delta",1),("B_sidis",1))},
                 **{"H_simc_shape_"+c:deepcopy(h) for c in observer.COMPONENTS})
        rows={p:[row(p,source=s,entry=j,weight=w,mm=mm,allcuts=ac,phi=phi)
                 for j,(s,w,mm,ac,phi) in enumerate(((0,2.,1.12,True,20.),(1,-2.,1.12,True,20.),
                    (2,-1.,.6,False,-1.),(3,0.,1.3,False,20.)))] for p in ("kaon","pion")}
        old={p:[dict(source=ev.SOURCES[r["source"]],entry=r["entry"],mm=r["mm"],ti=r["ti"],pi=r["pi"],
                     allcuts=r["allcuts"],nommcuts=r["nommcuts"],weight=r["weight"]) for r in rr] for p,rr in rows.items()}
        doc=dict(observations={p:{c:observer.project([dict(r,selected=r[c]) for r in rr],edges,3,9)
                                 for c in ev.CUTS} for p,rr in old.items()},cells=[])
        for ti in range(3):
            for pi in range(9):
                applications=observer.apply_transfer(fit,{c:doc["observations"]["pion"][c]["children"][ti][pi]
                                                        for c in ev.CUTS},inp())
                doc["cells"].append(dict(t_index=ti,phi_index=pi,applications={"Local_t":applications},
                                        status="fit_rejected" if pi==1 else "empty_child" if pi==2 else "fitted_local_templates"))
        weights=doc["cells"][0]["applications"]["Local_t"]["weights"]
        self.assertEqual(weights,[0.,2.,2.,2.,0.])
        for r in rows["pion"]:
            r["transfer"]=weights[r["transfer_bin"]]
            r["transfer_fingerprint"]=ev.digest(weights)
            r["weight"]=r["coefficient"]*r["transfer"];r["sumw2"]=r["weight"]**2
        tables={p:ev.columns(rr,p) for p,rr in rows.items()}
        serialized=json.loads(json.dumps(snapshot(doc),allow_nan=False))
        ev.validate_data_closure(tables,serialized)
        bad=deepcopy(serialized);bad["cells"][0]["applications"]["Local_t"]["cuts"]["allcuts"]["variances"][2]+=1.
        with self.assertRaisesRegex(ValueError,"Local_t"):ev.validate_data_closure(tables,bad)
        bad=deepcopy(serialized);bad["cells"][0]["applications"]["Local_t"]={"status":"unavailable","reason":"fit_missing"}
        with self.assertRaisesRegex(ValueError,"populated_transfer"):ev.validate_data_closure(tables,bad)
        # Empty support retains unavailable fit rather than inventing a transfer.
        serialized["cells"][1]["applications"]["Local_t"]={"status":"unavailable","reason":"fit_missing"}
        ev.validate_data_closure(tables,serialized)

    def test_actual_cache_join_exact_identities_flags_and_proton_rf(self):
        hist=dict(normfac_data=2.,normfac_dummy=1.,nWindows=1.,_yield_data_event_cache={},
                  _authoritative_pion_control_source_cache=dict(coordinate_fingerprint="coordinate",records=[]))
        for si,s in enumerate(ev.SOURCES):
            weight=(2.,-2.,-1.,1.)[si]
            hist["_authoritative_pion_control_source_cache"]["records"].append(dict(source_label=s,entry_index=0,
                adj_MM=1.12,adj_t=.5,t_index=0,phi_index=0,phi_degrees=20.,allcuts=True,nommcuts=True,
                coordinate_fingerprint="coordinate",coefficient=weight,source_tree_name="pion_"+s))
            hist["_yield_data_event_cache"][s]=dict(adj_MM=[1.12],adj_t=[.5],t_index=[0],phi_index=[0],
                phi_shift_deg=[20.],allcuts=[True],nommcuts=[True],source_label=["dummy_prompt" if s=="dummy" else s],
                source_entry_index=[0],proton_cleaning_factor=[.5],proton_rf_accept=[si!=3])
        before=snapshot(hist)
        rows,coeff=observer.rows_from_caches(hist)
        edges=meta()["axes"]["mm"]
        doc=dict(source_coefficients=coeff,observations={p:{c:observer.project([dict(r,selected=r[c]) for r in rr],edges,3,9)
                 for c in ev.CUTS} for p,rr in rows.items()},cells=[])
        for ti in range(3):
            for pi in range(9):
                doc["cells"].append(dict(t_index=ti,phi_index=pi,applications={"Local_t":dict(status="available",
                    weights=[0.,2.,2.,2.,0.],cuts={c:dict(values=[2.*v for v in doc["observations"]["pion"][c]["children"][ti][pi]["values"]],
                     variances=[4.*v for v in doc["observations"]["pion"][c]["children"][ti][pi]["variances"]]) for c in ev.CUTS})}))
        tables=ev.data_tables(hist,doc,rows,{s:"kaon_"+s for s in ev.SOURCES})
        self.assertEqual(tables["kaon"]["weight"].tolist(),[1.,-1.,-.5,0.])
        self.assertEqual(tables["pion"]["weight"].tolist(),[4.,-4.,-2.,2.])
        self.assertEqual(tables["pion"]["sumw2"].tolist(),[16.,16.,4.,4.])
        self.assertEqual(before,snapshot(hist))
        for pop in ("kaon","pion"):ev.validate_rows(tables[pop],pop,meta())
        bad=deepcopy(rows);bad["kaon"][0]["entry"]=9
        with self.assertRaisesRegex(ValueError,"observer_cache_identity"):
            ev.data_tables(hist,doc,bad,{s:"kaon_"+s for s in ev.SOURCES})
        for field,value in (("mm",1.13),("allcuts",False),("weight",10.)):
            bad=deepcopy(rows);bad["kaon"][0][field]=value
            with self.assertRaisesRegex(ValueError,"cache_scalar_or_flag_mismatch"):
                ev.data_tables(hist,doc,bad,{s:"kaon_"+s for s in ev.SOURCES})

    def test_real_simc_loop_sink_and_disabled_preservation(self):
        from binning_helpers import find_2d_bin_indices
        source=ast.parse((ROOT/"src/binning/calculate_yield.py").read_text())
        functions=[n for n in source.body if isinstance(n,ast.FunctionDef) and n.name in
                   ("_process_yield_simc_tree","_init_ave_simc_event_cache","_append_ave_simc_event")]
        ns=dict(math=math,np=np,find_2d_bin_indices=find_2d_bin_indices,
                calculate_theta_cm_deg=lambda *a:45.,_get_simc_true_theta_cm_deg=lambda *a:45.,
                _fill_t_vs_tmin_hist=lambda *a,**k:None)
        exec(compile(ast.Module(body=functions,type_ignores=[]),"actual_calculate_yield_seam","exec"),ns)
        class Hist:
            def __init__(self):self.fills=[]
            def Fill(self,*args):self.fills.append(args)
        class Tree(list):
            def GetEntries(self):return len(self)
        events=Tree([types.SimpleNamespace(t=-.5,phipq=math.radians(20.),missmass=mm,iter_weight=w,
                      W=2.74,Q2=4.4,epsilon=.5,ssxptar=0.,ssyptar=0.,hsxptar=0.,hsyptar=0.,
                      phgcer_x_det=hole,phgcer_y_det=0.) for mm,w,hole in
                     ((1.12,2.,0.),(.8,-1.,0.),(1.12,4.,1.))])
        def run(sink):
            group={k:[[Hist() for _ in range(9)] for _ in range(3)] for k in
                   ("t","Q2","W","q2_w","theta_cm","theta_cm_true","ssxptar","ssyptar",
                    "hsxptar","hsyptar","t_vs_tmin","mm","mm_unweighted")}
            cache=ns["_init_ave_simc_event_cache"]();calls=[]
            def allcut(e,lo,hi):calls.append(("all",e.missmass));return lo<=e.missmass<hi
            def nocut(e,lo,hi):calls.append(("no",e.missmass));return True
            ns["_process_yield_simc_tree"](events,group,cache,"kaon","",lambda x,y:bool(x),
                allcut,1.1,1.16,meta()["axes"]["canonical_t"],meta()["axes"]["phi"],lambda *a,**k:None,
                **({"check1_sink":sink,"apply_simc_sub_cuts":nocut} if sink is not None else {}))
            return {k:[[h.fills for h in rr] for rr in v] for k,v in group.items()},cache,calls
        old,cache,calls=run(None)
        self.assertTrue(all(c[0]=="all" for c in calls))
        sink=ev.SimcSink();new,new_cache,new_calls=run(sink)
        self.assertEqual(old,new);self.assertEqual(cache,new_cache)
        self.assertEqual(len([c for c in new_calls if c[0]=="no"]),3)
        self.assertEqual([r["mm"] for r in sink.rows],[1.12,.8])
        self.assertEqual([r["allcuts"] for r in sink.rows],[True,False])
        self.assertEqual(sink.table(3.)["weight"].tolist(),[6.,-3.])
        self.assertEqual(new["mm"][0][0],[(1.12,2.)])

    def test_native_histogram_closure_and_export_failure_preserves_objects(self):
        m=meta("simc");sink=ev.SimcSink()
        sink.observe(0,1.12,.5,20.,2.,True,True,0,0)
        table=sink.table(3.)
        mm_edges=m["axes"]["native_mm"];t_edges=m["axes"]["t"]
        hists={"mm":[[H(values=[0.]*4,errors=[0.]*4,edges=mm_edges) for _ in range(9)] for _ in range(3)],
               "t":[[H(values=[0.]*102,errors=[0.]*102,edges=t_edges) for _ in range(9)] for _ in range(3)]}
        for key,value in (("mm",1.12),("t",.5)):
            b=ev.bisect_right(mm_edges if key=="mm" else t_edges,value)
            hists[key][0][0].values[b]=6.;hists[key][0][0].errors[b]=6.
        before=snapshot(hists)
        ev.native_simc_closure(table,m["axes"],hists)
        bad=deepcopy(table);bad["weight"]*=3.;bad["sumw2"]*=9.
        with self.assertRaisesRegex(ValueError,"native_simc"):ev.native_simc_closure(bad,m["axes"],hists)
        self.assertEqual(before,snapshot(hists))
        hists["mm"][0][0].edges=[1.1,1.14,1.16]
        with self.assertRaisesRegex(ValueError,"native_mm_axis"):ev.native_simc_closure(table,m["axes"],hists)


class FiveOwnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory()
        cls.base=Path(cls.temp.name)
        cls.out,cls.identities,cls.runtime,cls.reference,cls.authority,cls.histograms=five_fixture(cls.base)

    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()

    def verify(self):
        # Exercise actual existing JSON/source/numerical validators without
        # fabricating an unavailable Poppler result. The separate full-wrapper
        # test below owns the real PDF-page gate.
        records={};outcomes={};legacy={p:self.identities[p] for p in owner.child_source_inventory(ROOT)}
        for name in owner.names()|owner.names(True):
            p=self.out/name
            owner.require(p.is_file() and not p.is_symlink() and p.stat().st_mtime_ns>=1,'original_missing_or_stale')
            records[name]=dict(basename=name,bytes=p.stat().st_size,sha256=owner.sha256(p),mtime_ns=p.stat().st_mtime_ns)
        for setting,eps in ev.SETTINGS:
            stem=f'{setting}_kaon_Q4p4W2p74_{eps}_pion-t-template-comparison'
            parent=owner.strict_json((self.out/(stem+'.json')).read_text())
            child=owner.strict_json((self.out/(stem.replace('pion-t-','pion-t-phi-')+'.json')).read_text())
            high=eps=='highe'
            outcome=owner.verify_payload(parent,setting,'a'*40,legacy,records[stem+'.pdf'],high=high)
            outcome['child']=owner.verify_child_payload(child,setting,'a'*40,legacy,parent,self.runtime,high=high)
            outcomes.setdefault(eps,{})[setting]=outcome
        return owner._verify_check1_packets(ev,self.out,'a'*40,self.identities,1,
            dict.fromkeys(owner.names(check1=True)),self.runtime,self.reference,self.authority,records,outcomes)

    def test_whole_five_setting_owner_with_real_poppler_page_gate(self):
        if not shutil.which('pdfinfo'):self.skipTest('NOT VERIFIED: full PDF-inclusive five-setting owner requires real Poppler pdfinfo')
        records,_=owner.verify_check1_artifacts(ev,ROOT,self.out,'a'*40,self.identities,1,
            dict.fromkeys(owner.names(check1=True)),self.runtime,self.reference,self.authority)
        self.assertEqual(len(records),35)

    def test_five_setting_numerical_admission_real_profile_and_collector_roundtrip(self):
        from testing import collect_pion_hgcer_validation_bundle as collector
        records,outcomes=self.verify()
        self.assertEqual(len(records),35);self.assertEqual(set(outcomes),{'lowe','highe'})
        profile=owner.resolved_profile(ROOT,'a'*40,check1=True)
        with tempfile.TemporaryDirectory() as tmp:
            effective=Path(tmp)/'profile.json';effective.write_bytes(ev.canonical_json(profile))
            output=Path(tmp)/'pending.zip'
            def runner(command,cwd):
                stdout=('a'*40+'\n') if list(command)==['git','rev-parse','HEAD'] else ''
                return dict(command=list(command),returncode=0,stdout=stdout,stderr='')
            result=collector.collect_validation_bundle(outdir=self.out,kinematic='Q4p4W2p74',output=output,
                profile_path=effective,repo_root=ROOT,command_runner=runner)
            self.assertEqual(result['returncode'],0,result['manifest']['errors'])
            owner.verify_zip(output,'a'*40,records,check1=True)
            import zipfile
            with zipfile.ZipFile(output) as archive:contents={n:archive.read(n) for n in archive.namelist()}
            for kind in ('member','manifest'):
                modified=dict(contents)
                if kind=='member':modified['global/Center_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json']+=b' '
                else:
                    manifest=json.loads(modified['manifest.json']);manifest['global_artifacts'][next(iter(manifest['global_artifacts']))]['sha256']='0'*64
                    modified['manifest.json']=ev.canonical_json(manifest)
                bad=Path(tmp)/(kind+'.zip')
                with zipfile.ZipFile(bad,'w') as archive:
                    for name,raw in modified.items():archive.writestr(name,raw)
                with self.assertRaisesRegex(ValueError,'zip_artifact_identity_mismatch'):owner.verify_zip(bad,'a'*40,records,check1=True)
            with zipfile.ZipFile(output,'a') as archive:archive.writestr('global/unapproved',b'tampered')
            with self.assertRaisesRegex(ValueError,'zip_inventory'):owner.verify_zip(output,'a'*40,records,check1=True)

    def test_missing_low_high_sidecar_stale_and_tampered_fail_closed(self):
        for name in ('Left_kaon_Q4p4W2p74_lowe_pion-t-phi-template-comparison.json',
                     'Right_kaon_Q4p4W2p74_highe_pion-t-phi-template-comparison.json',
                     'Center_kaon_Q4p4W2p74_lowe_check1_data.npz'):
            path=self.out/name;raw=path.read_bytes();path.unlink()
            try:
                with self.assertRaises((ValueError,FileNotFoundError)):self.verify()
            finally:path.write_bytes(raw)
        path=self.out/'Center_kaon_Q4p4W2p74_lowe_check1_data.npz';raw=path.read_bytes()
        path.write_bytes(raw+b'tamper')
        try:
            with self.assertRaisesRegex(ValueError,'npz_identity'):self.verify()
        finally:path.write_bytes(raw)
        path=self.out/'Center_kaon_Q4p4W2p74_lowe_check1_data.json';stat=path.stat()
        os.utime(path,ns=(0,0))
        try:
            with self.assertRaisesRegex(ValueError,'stale'):self.verify()
        finally:os.utime(path,ns=(stat.st_atime_ns,stat.st_mtime_ns))

    def test_mixed_attempt_and_native_reference_mismatch_fail(self):
        p=self.out/'Center_kaon_Q4p4W2p74_lowe_check1_data.json';raw=p.read_bytes();doc=json.loads(raw)
        doc['attempt']='other-attempt';p.write_bytes(ev.canonical_json(doc))
        try:
            with self.assertRaisesRegex(ValueError,'packet_binding'):self.verify()
        finally:p.write_bytes(raw)
        doc=json.loads(raw);doc['lineage']['reference_source_sha']='0'*40;p.write_bytes(ev.canonical_json(doc))
        try:
            with self.assertRaisesRegex(ValueError,'packet_binding'):self.verify()
        finally:p.write_bytes(raw)
        doc=json.loads(raw);doc['axes']['canonical_t'][1]=.58;p.write_bytes(ev.canonical_json(doc))
        try:
            with self.assertRaisesRegex(ValueError,'packet_axes'):self.verify()
        finally:p.write_bytes(raw)
        a=self.reference['simc']['lowe']['simc_mm_center_values'];a[0,0,0]+=1.
        try:
            with self.assertRaisesRegex(ValueError,'reference_simc_mm'):self.verify()
        finally:a[0,0,0]-=1.

    def test_source_roles_original_runtime_and_input_admission(self):
        plan=[dict(source=str(self.out/n)) for n in self.runtime]
        admitted=owner.check1_admit_inputs(ev,self.reference,plan)
        self.assertEqual(len(admitted),30)
        n=next(iter(self.reference['summaries']['lowe']['original_inputs']))
        before=self.reference['summaries']['lowe']['original_inputs'][n]['original']['sha256']
        self.reference['summaries']['lowe']['original_inputs'][n]['original']['sha256']='0'*64
        try:
            with self.assertRaisesRegex(ValueError,'cohort_mismatch'):owner.check1_admit_inputs(ev,self.reference,plan)
        finally:self.reference['summaries']['lowe']['original_inputs'][n]['original']['sha256']=before

    def test_producer_archive_model_role_closure(self):
        self.verify()
        for setting,eps in ev.SETTINGS:
            roles=ev.simc_derivation(setting,eps)
            doc,_=ev.read_packet(self.out/f'{setting}_kaon_Q4p4W2p74_{eps}_check1_simc.json')
            binding=doc['source_bindings']['simc:simc']
            self.assertEqual(binding['basename'],f'{setting}_kaon_Simc_Q4p4W2p74_{eps}.root')
            self.assertEqual(binding['original_basename'],roles['archive'])
            self.assertEqual(binding['original']['sha256'],ev.file_hash(self.out/roles['archive']))
            self.assertEqual(binding['runtime']['sha256'],ev.file_hash(self.out/roles['model']))
            self.assertNotEqual(binding['original']['sha256'],binding['runtime']['sha256'])
            self.assertEqual(binding['runtime_hist']['runtime'],binding['archive_hist']['original'])
            self.assertEqual(binding['derivation'],self.authority['simc_derivations'][roles['model']])
            for source in roles['sources']:
                self.assertEqual(doc['executables'][source]['sha256'],ev.file_hash(ROOT/source))

    def test_simc_archive_and_hist_admission_remain_mandatory(self):
        roles=ev.simc_derivation('Left','lowe');plan=[dict(source=str(self.out/n)) for n in self.runtime]
        for role in ('archive','archive_hist'):
            name=roles[role]
            with self.subTest(missing=name),self.assertRaisesRegex(ValueError,'reference_input_missing'):
                owner.check1_admit_inputs(ev,self.reference,[p for p in plan if Path(p['source']).name!=name])
            expected=self.reference['summaries']['lowe']['original_inputs'][name]['original']
            saved=expected['sha256'];expected['sha256']=self.runtime[roles['model']]['sha256']
            try:
                with self.subTest(substituted=name),self.assertRaisesRegex(ValueError,'original_cohort_mismatch'):
                    owner.check1_admit_inputs(ev,self.reference,plan)
            finally:expected['sha256']=saved

    def test_missing_mismatched_substituted_archive_and_model_identities_rejected(self):
        roles=ev.simc_derivation('Left','lowe')
        for role in ('archive','archive_hist','model','model_hist'):
            name=roles[role];saved=self.runtime.pop(name)
            try:
                with self.subTest(missing=role),self.assertRaisesRegex(ValueError,'archive_binding|simulation_identity_missing'):
                    self.verify()
            finally:self.runtime[name]=saved
            saved=self.runtime[name];self.runtime[name]=dict(saved,sha256='0'*64)
            try:
                with self.subTest(mismatch=role),self.assertRaisesRegex(ValueError,'archive_binding|model_binding|archive_copy_identity'):
                    self.verify()
            finally:self.runtime[name]=saved
        model=roles['model'];saved=self.runtime[model]
        self.runtime[model]=self.runtime[roles['archive']]
        try:
            with self.assertRaisesRegex(ValueError,'simc_model_binding'):self.verify()
        finally:self.runtime[model]=saved
        for name in (roles['archive'],roles['archive_hist']):
            saved=self.authority['original_inputs'].pop(name)
            try:
                with self.assertRaisesRegex(ValueError,'archive_binding'):self.verify()
            finally:self.authority['original_inputs'][name]=saved
        saved=self.authority['simc_derivations'][model]
        self.authority['simc_derivations'][model]=dict(saved,archive=roles['model'])
        try:
            with self.assertRaisesRegex(ValueError,'derivation_authority'):self.verify()
        finally:self.authority['simc_derivations'][model]=saved
        # A same-basename replacement is detected by actual owner readback,
        # independently of the producer's retained consumed-model hash.
        path=self.out/model;raw=path.read_bytes();saved=self.runtime
        try:
            path.write_bytes(b'substituted model')
            self.runtime=owner.check1_runtime_models(ev,self.out,self.runtime)
            with self.assertRaisesRegex(ValueError,'simc_model_binding'):self.verify()
        finally:path.write_bytes(raw);self.runtime=saved

    def test_export_rejects_archive_path_missing_model_and_substituted_hist(self):
        setting,eps='Left','lowe';roles=ev.simc_derivation(setting,eps)
        sink,hists=self.histograms[setting+':'+eps]
        ap=self.base/'authority.json'
        def export(name):
            return ev.export_simc(sink,types.SimpleNamespace(GetName=lambda:str(self.out/name)),3.,inp(eps),
                setting,self.reference['children'][setting+':'+eps]['t_edges'],
                self.reference['children'][setting+':'+eps]['phi_edges'],hists,self.out)
        with patch.dict(os.environ,{ev.AUTHORITY_ENV:str(ap),ev.AUTHORITY_HASH_ENV:ev.file_hash(ap)}):
            with self.assertRaisesRegex(ValueError,'simulation_runtime_role'):export(roles['archive'])
            for role in ('model','model_hist'):
                path=self.out/roles[role];raw=path.read_bytes();path.unlink()
                try:
                    with self.assertRaisesRegex(ValueError,'simulation_runtime_file'):export(roles['model'])
                    with self.assertRaisesRegex(ValueError,'runtime_model_missing_or_alias'):
                        owner.check1_runtime_models(ev,self.out,self.runtime)
                finally:path.write_bytes(raw)
            path=self.out/roles['model_hist'];raw=path.read_bytes();path.write_bytes(b'wrong hist')
            try:
                with self.assertRaisesRegex(ValueError,'archive_copy_identity'):export(roles['model'])
            finally:path.write_bytes(raw)

    def test_external_template_admission_is_bound_to_same_setting_and_component(self):
        with tempfile.TemporaryDirectory() as temp:
            path=Path(temp)/'template.root';path.write_bytes(b'controlled pion template')
            reference=deepcopy(self.reference)
            for p in reference['parents']['Left:lowe']['parents']:
                p['templates']['pi_n']['provenance']['sha256']=ev.file_hash(path)
            inventory=[dict(path=str(path),role='pion_template:neutron:Left:lowe',copied=False)]
            self.assertEqual(len(owner.check1_admit_external(reference,inventory)),1)
            path.write_bytes(b'changed model')
            with self.assertRaisesRegex(ValueError,'template_reference_mismatch'):owner.check1_admit_external(reference,inventory)

    def test_no_clobber_final_release_and_companion_failures(self):
        for mode in ('success','collision','companion_failure','verification_failure','preservation_failure','cleanup_failure',
                     'cleanup_exception','diagnostic_write_failure'):
            with self.subTest(mode=mode),tempfile.TemporaryDirectory() as temp:
                root=Path(temp);art=root/'art';transfer=root/'transfer';art.mkdir();transfer.mkdir()
                output=transfer/'result.zip';pending=transfer/'result.pending.zip';pending.write_bytes(b'validated pending bytes')
                info=pending.stat();pending_identity=(info.st_dev,info.st_ino)  # fixture-owned creation
                status=owner.GateStatus(art,output,'a'*40,high=True,check1=True)
                status.update(zip_identity=dict(sha256=owner.sha256(pending),bytes=pending.stat().st_size))
                status.update(artifact_verification_completed=True,zip_verification_completed=True,
                              worktree_cleanup_completed=True,check1_preservation_completed=True)
                paths=owner.destinations(art,output);paths['log'].write_bytes(b'child log');paths['summary'].write_bytes(b'bound run summary')
                if mode=='collision':output.write_bytes(b'competitor')
                if mode in ('verification_failure','preservation_failure','cleanup_failure'):
                    field={'verification_failure':'zip_verification_completed','preservation_failure':'check1_preservation_completed',
                           'cleanup_failure':'worktree_cleanup_completed'}[mode]
                    status.update(**{field:False})
                    with self.assertRaisesRegex(ValueError,'release_gate_incomplete'):owner.release_check1(status,paths,pending,output,pending_identity)
                    self.assertFalse(output.exists())
                elif mode in ('cleanup_exception','diagnostic_write_failure'):
                    temporary_directory=owner.tempfile.TemporaryDirectory
                    @contextmanager
                    def failing_cleanup(**kwargs):
                        with temporary_directory(**kwargs) as directory:yield directory
                        raise OSError('injected post-link cleanup failure')
                    real_update=status.update
                    def update(*args,**kwargs):
                        if mode=='diagnostic_write_failure' and kwargs.get('status')=='failed':
                            raise OSError('injected receipt write failure')
                        return real_update(*args,**kwargs)
                    with patch.object(owner.tempfile,'TemporaryDirectory',side_effect=failing_cleanup),patch.object(status,'update',side_effect=update):
                        with self.assertRaisesRegex(OSError,'post-link cleanup failure'):owner.release_check1(status,paths,pending,output,pending_identity)
                    self.assertFalse(output.exists())
                    self.assertIn('post-link cleanup failure',status.record['first_failed_invariant'])
                    if mode=='diagnostic_write_failure':self.assertIn('diagnostic_write_failure',status.record)
                elif mode=='companion_failure':
                    with patch.object(owner.os,'link',side_effect=OSError('injected companion failure')):
                        with self.assertRaisesRegex(OSError,'injected'):owner.release_check1(status,paths,pending,output,pending_identity)
                elif mode=='collision':
                    with self.assertRaises(FileExistsError):owner.release_check1(status,paths,pending,output,pending_identity)
                    self.assertEqual(output.read_bytes(),b'competitor')
                else:
                    owner.release_check1(status,paths,pending,output,pending_identity)
                    self.assertEqual(output.read_bytes(),b'validated pending bytes');self.assertFalse(pending.exists())
                    self.assertEqual(json.loads(paths['delivered_receipt'].read_text())['zip_identity']['sha256'],owner.sha256(output))
                if mode!='success':
                    self.assertEqual(status.record['status'],'failed')
                    self.assertEqual(status.record['failed_stage'],'check1_release')
                    self.assertFalse(paths['delivered_log'].exists());self.assertFalse(paths['delivered_summary'].exists())
                    self.assertFalse(paths['delivered_receipt'].exists())


class ExternalReferenceTests(unittest.TestCase):
    def test_failure_diagnostic_on_immutable_reference_structure_when_supplied(self):
        path=os.environ.get('KAONLT_CHECK1_REFERENCE_TEST_PATH')
        if not path:self.skipTest('supply read-only pinned ZIP through KAONLT_CHECK1_REFERENCE_TEST_PATH')
        reference=ev.load_reference(path)
        expected=reference['children']['Left:lowe'];actual=deepcopy(expected)
        parent=actual['parents'][0]['baseline_A']['pion_component_alignment']
        parent['components']['pi_delta']['immutable_source_hist_checksum']='e'*64
        parent['parent_alignment_hash']=ev.digest(ev.alignment_parent_preimage(parent))
        actual['cells'][0]['fit']['pion_component_alignment']['applied_components']['pi_delta']['parent_alignment_hash']=parent['parent_alignment_hash']
        diag=ev.scientific_failure_diagnostic(actual,expected,'a'*40,ev.REFERENCE_SHA256,'Left:lowe',ScientificFailureDiagnosticsTests.invariant)
        self.assertEqual(diag['first_difference']['path'],ScientificFailureDiagnosticsTests.invariant.split(':',2)[2])
        # Stored reference identity and a recomputed current-producer preimage
        # are distinct evidence. Verify reproduction, never force equality.
        function=next(n for n in ast.parse((ROOT/'src/cuts/pion_component_fits.py').read_bytes()).body
            if isinstance(n,ast.FunctionDef) and n.name=='_alignment_parent_hash')
        scope=dict(_alignment_hash=ev.digest,COMPONENT_NAMES=('pi_n','pi_delta','pi_sidis'))
        exec(compile(ast.Module(body=[function],type_ignores=[]),'<frozen-parent-hash>','exec'),scope)
        saved=expected['parents'][0]['baseline_A']['pion_component_alignment']
        recomputed=scope['_alignment_parent_hash'](saved['parent_setting_key'],saved['common_setting_shift_gev'],
            saved['components'],saved['resolved_configuration_hash'])
        self.assertEqual(diag['expected']['parent_hash_consistency']['recomputed_hash'],recomputed)
        self.assertEqual(diag['expected']['parent_hash_consistency']['status'],
            'matches' if recomputed==saved['parent_alignment_hash'] else 'mismatch')
        self.assertEqual(diag['expected']['child']['parent_alignment_hash']['value'],
            'c03775871b5033a2d9682e174293ecefe599ebefb1a2e9a10cbc102dd223c82a')
        self.assertTrue(diag['expected']['parent']['selection']['accepted']['value'])
        self.assertFalse(diag['expected']['child']['selection']['accepted']['value'])
        self.assertEqual(diag['first_parent_operand_difference']['path'],
            'parent_hash_operands.components.pi_delta.immutable_source_hist_checksum')

    def test_real_pinned_reference_when_supplied(self):
        path=os.environ.get('KAONLT_CHECK1_REFERENCE_TEST_PATH')
        if not path:self.skipTest('supply read-only pinned ZIP through KAONLT_CHECK1_REFERENCE_TEST_PATH')
        reference=ev.load_reference(path)
        self.assertEqual(set(reference['children']),{s+':'+e for s,e in ev.SETTINGS})
        self.assertEqual(sum(len(c['cells']) for c in reference['children'].values()),135)
        for d in reference['children'].values():
            ev.compare_science(ev.scientific_projection(d),ev.scientific_projection(deepcopy(d)))
        for eps in ('lowe','highe'):
            name='Center_kaon_Analysed_Data_Q4p4W2p74_'+eps+'.root'
            summary=reference['summaries'][eps]
            self.assertNotEqual(summary['original_inputs'][name]['original']['sha256'],summary['runtime_container_identities'][name]['sha256'])
            self.assertEqual(reference['children']['Center:'+eps]['runtime_inputs']['data']['sha256'],
                             summary['runtime_container_identities'][name]['sha256'])
        d=deepcopy(reference['children']['Center:lowe']);d['source_coefficients']['prompt']*=2
        with self.assertRaisesRegex(ValueError,'scientific_value'):
            ev.compare_science(ev.scientific_projection(d),ev.scientific_projection(reference['children']['Center:lowe']))
        self.assertEqual(ev.scientific_projection(dict(root_filename='/old/input.root',artifact_paths=['/old/alignment.json'])),
                         ev.scientific_projection(dict(root_filename='/private/input.root',artifact_paths=['/private/alignment.json'])))
        from data_coordinates import coordinate_fingerprint
        original=reference['children']['Center:lowe'];d=deepcopy(original)
        oldfp=d['coordinate']['coordinate_fingerprint']
        d['coordinate']['mm_shift_provenance']['plot_filename']='/private/'+Path(d['coordinate']['mm_shift_provenance']['plot_filename']).name
        newfp=coordinate_fingerprint(d['coordinate'])
        def replace(value):
            if isinstance(value,dict):return {k:replace(v) for k,v in value.items()}
            if isinstance(value,list):return [replace(v) for v in value]
            return newfp if value==oldfp else value
        d=replace(d);ev.validate_coordinate_contract(d)
        ev.compare_science(ev.scientific_projection(d),ev.scientific_projection(original))
        corrupt=deepcopy(d)
        corrupt['parents'][0]['baseline_A']['pion_component_alignment']['bin_key']['kaon_data_coordinate_fingerprint']='0'*64
        with self.assertRaisesRegex(ValueError,'coordinate_role_mismatch'):ev.validate_coordinate_contract(corrupt)
        with self.assertRaisesRegex(ValueError,'coordinate_role_mismatch'):ev.scientific_projection(corrupt)
        d['coordinate']['mm_shift']+=.001
        with self.assertRaisesRegex(ValueError,r'scientific_identity:reference.cells\[0\].fit.pion_component_alignment.bin_key.coordinate_fingerprint'):
            ev.compare_science(ev.scientific_projection(d),ev.scientific_projection(original))
        with tempfile.TemporaryDirectory() as temp:
            bad=Path(temp)/'unapproved.zip';bad.write_bytes(b'wrong packet')
            with self.assertRaisesRegex(ValueError,'reference_digest'):ev.load_reference(bad)


if __name__ == "__main__":
    unittest.main()
