"""Complete detached pion statistical/model-response propagation; no new fit."""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import zipfile
import numpy as np
if __package__:
    from . import bootstrap_pion_mm_delta as boot
else:
    import bootstrap_pion_mm_delta as boot
fit=boot.fit


def require(ok,reason):
    if not ok:raise ValueError('pion_propagation:'+reason)


def verified_bundle(path):
    with zipfile.ZipFile(path) as z:
        require(len(z.namelist())==len(set(z.namelist())),'duplicate_member')
        manifest=json.loads(z.read('members-sha256.json'))
        require(set(z.namelist())==set(manifest)|{'members-sha256.json'},'unbound_member')
        for name,info in manifest.items():
            b=z.read(name)
            require((len(b),fit.digest(b))==(info['bytes'],info['sha256']),'member_hash')
        result=json.loads(z.read('result.json'))
        interface=json.loads(z.read('covariance_interface.json'))
        raw=z.read('replicas_and_covariance_factors.npz')
    with np.load(io.BytesIO(raw),allow_pickle=False) as a:arrays={k:a[k] for k in a.files}
    return result,interface,arrays,raw


def supported_envelope(central,values,support):
    """Central-inclusive range, never a probability-weighted covariance."""
    central=np.asarray(central,dtype=float);values=np.asarray(values,dtype=float)
    support=np.asarray(support,dtype=bool)
    require(values.shape==(len(support),)+central.shape,'envelope_shape')
    require(np.isfinite(central).all() and np.isfinite(values[support]).all(),'nonfinite_envelope')
    if not support.any():return dict(status='unavailable_no_supported_variations',minimum=None,maximum=None,supported_cases=0,omitted_cases=len(support))
    x=np.concatenate((central[None],values[support]),axis=0)
    return dict(status='sensitivity_range_not_one_sigma',minimum=x.min(0).tolist(),maximum=x.max(0).tolist(),
                supported_cases=int(support.sum()),omitted_cases=int((~support).sum()))


def correlation(C):
    C=np.asarray(C,dtype=float);s=np.sqrt(np.maximum(0,np.diag(C)))
    return [[float(C[i,j]/(s[i]*s[j])) if s[i]>0 and s[j]>0 else None for j in range(len(s))] for i in range(len(s))]


def aggregate_factor(F,weights):
    F=np.asarray(F,dtype=float);weights=np.asarray(weights,dtype=float)
    require(F.ndim==2 and weights.ndim==2 and F.shape[1]==weights.shape[0],'factor_axis')
    require(np.isfinite(F).all() and np.isfinite(weights).all(),'nonfinite_factor')
    projected=F@weights
    return projected.T@projected


def response_cosine(vectors):
    """Co-movement geometry only: no scenario probabilities are assigned."""
    D=np.asarray(vectors,dtype=float);require(D.ndim==2 and np.isfinite(D).all(),'response_shape')
    G=D.T@D
    return correlation(G)


def varied_model(q,case):
    label=case['kind'];X=q['Xfull'].copy()
    if label.startswith('omit_delta_'):
        index=int(label[len('omit_delta_'):]);require(1<=index<=10,'delta_case')
    elif label.startswith('alignment_'):
        channels=[(ci,ch) for ci,ch in enumerate(fit.CHANNELS) if label.startswith('alignment_'+ch+'_')]
        require(len(channels)==1,'alignment_channel')
        ci,ch=channels[0];shift=float(label[len('alignment_'+ch+'_'):]);require(np.isfinite(shift),'alignment_shift')
        ref=deepcopy(q['parent_refs'][ci]);ref['shift']=shift
        L=boot.extraction.transform(ref);X[:,:,ci]=q['raw_parts'][ci]@L.T
        fit.close(X[:,:,ci][:,1:-1].sum(),1.,'propagated_parent_normalization')
    else:raise ValueError('pion_propagation:unapproved_model_family')
    beta=np.asarray(case['beta'],dtype=float)
    require(beta.shape==(3,) and np.isfinite(beta).all() and np.all(beta>=0),'case_amplitudes')
    return X,beta


def regional(q,X,beta):
    ed=q['edges'];N=len(ed)+1
    masks={'broad':np.r_[False,np.ones(N-2,bool),False],
           'lambda_complete':np.r_[False,(ed[:-1]>=1.10)&(ed[1:]<=1.16),False],
           'high_MM_complete':np.r_[False,ed[:-1]>=1.16,False],
           'low_MM_complete':np.r_[False,ed[1:]<=1.10,False]}
    model=X@beta;result={}
    for name,mask in masks.items():
        y=q['y'][:10,mask].sum(1);m=model[:10,mask].sum(1)
        result[name]={'data_delta':y.tolist(),'model_delta':m.tolist(),'residual_delta':(y-m).tolist(),
                      'data_total':float(y.sum()),'model_total':float(m.sum()),'residual_total':float((y-m).sum())}
    return result


def execute(args):
    out=Path(args.output);require(not out.exists(),'output_exists')
    old,interface,replicas,raw=verified_bundle(args.bootstrap_review)
    require(old['replicas']==1024 and old['full_scope'] and not old['production_applied'],'bootstrap_scope')
    require(old['source_provenance']['bootstrap_implementation_sha256']==fit.digest(Path(boot.__file__).read_bytes()),'bootstrap_source_identity')
    require(old['source_provenance']['current_implementation_sha256']==fit.digest(Path(fit.__file__).read_bytes()),'fit_source_identity')
    observed,models,roles,docs,accepted,rows,original,policies,A,provenance=boot.load_evidence(args.joint_review,args.fixed_comparison,args.projection,args.acquisition,args.tenbin)
    for k in ('fixed_archive_sha256','projection_sha256','historical_reference_sha256','joint_review_sha256','tenbin_sha256'):
        require(provenance[k]==old['source_provenance'][k],'input_binding_'+k)
    lookup={(p['setting'],p['parent']):p for p in accepted['parents']}
    order=interface['axes']['parent_order'];P=len(order)
    require(P==15 and {(p['setting'],p['parent']) for p in order}==set(lookup),'parent_axes')
    require(order==[{'setting':p['setting'],'parent':p['parent']} for p in old['parents']],'bootstrap_parent_order')
    native=interface['axes']['native_background_cells'];F=replicas['background_covariance_factor']
    indices=replicas['background_supported_replica_indices']
    bgreps=replicas['background_joint_replicas']
    require(F.shape[0]==len(indices)>=2 and np.all(np.diff(indices)>0) and np.all((indices>=0)&(indices<1024)),'replica_indices')
    fit.close(F,(bgreps[indices]-bgreps[indices].mean(0))/np.sqrt(len(indices)-1),'retained_covariance_factor')
    require(len(native)==15*10*102,'native_axes')
    columns=[v['factor_column'] for v in native if v['factor_column'] is not None]
    require(sorted(columns)==list(range(F.shape[1])),'factor_columns')
    stored_responses={(p['parent_slot'],p['kind']):p for p in old['model_responses']}
    summaries=[];cases=[];case_arrays=[];central_grids=[];central_phis=[];central_denoms=[];support_grids=[]
    mm_cov=[];phi_cov=[];phi_lambda_cov=[];denom_cov=[];denom_bg_cross=[];checks=[]
    response_arrays=replicas['model_background_response_vectors']
    for ip,item in enumerate(order):
        setting,parent=item['setting'],item['parent'];q=fit.build_problem(setting,parent-1,observed,models,roles,docs)
        fit.close(q['edges'],interface['axes']['MM_edges'],'native_MM_axis')
        source=lookup[setting,parent];beta=np.asarray(source['beta']);X=q['Xfull'];p=rows[setting]
        ones=np.ones(len(p['mm']));D=X.sum(0)@beta;num=X.sum(0)@A[setting,parent]
        grid,lp,supported,ratio=boot.background(p,ones,parent,q['edges'],num,D,policies[setting,parent])
        fit.close(beta,old['parents'][ip]['central']['beta'],'central_beta')
        fit.close(D,old['parents'][ip]['central']['denominator'],'central_denominator')
        fit.close(lp,old['parents'][ip]['central']['lambda_phi'],'central_Lambda_phi')
        require(supported['broad'] and supported['lambda_region'],'central_support')
        mapping=[v for v in native if v['parent_slot']==ip]
        require(len(mapping)==1020,'parent_cell_count')
        Wmm=np.zeros((F.shape[1],102));Wphi=np.zeros((F.shape[1],10));actual=np.zeros((10,102),int)
        use=(p['ti']==parent-1)&p['nommcuts'];mb=np.searchsorted(q['edges'],p['mm'],side='right');ph=np.where((p['pi']>=0)&(p['pi']<9),p['pi'],9)
        np.add.at(actual,(ph[use],mb[use]),1)
        active=[]
        for cell in mapping:
            ph0,mm=cell['phi_bin'],cell['mm_bin'];col=cell['factor_column']
            require((cell['setting'],cell['t_bin'])==(setting,parent),'parent_cell_identity')
            require(cell['application_event_support']==bool(actual[ph0,mm]),'application_cell_support')
            if col is not None:
                require(cell['application_event_support'],'factor_empty_cell');Wmm[col,mm]=1;Wphi[col,ph0]=1;active.append((col,ph0*102+mm))
        mm_cov.append(aggregate_factor(F,Wmm));phi_cov.append(aggregate_factor(F,Wphi))
        phi_lambda_cov.append(boot.covariance(replicas['joint/lambda_phi_background'][:,ip]))
        denom=replicas['denominator_joint_replicas'][:,ip];require(np.isfinite(denom).all(),'denominator_replicas')
        denom_cov.append(boot.covariance(denom))
        dd=denom[indices]-denom[indices].mean(0);bb=bgreps[indices]@Wmm
        denom_bg_cross.append(dd.T@(bb-bb.mean(0))/(len(indices)-1))
        native_complete=np.r_[False,(q['edges'][:-1]>=1.1)&(q['edges'][1:]<=1.16),False]
        central_region=regional(q,X,beta);local=[]
        for case in source['sensitivities']:
            if not (case['kind'].startswith('alignment_') or case['kind'].startswith('omit_delta_')):continue
            Xt,bv=varied_model(q,case);dv=Xt.sum(0)@bv;nv=Xt.sum(0)@A[setting,parent]
            gv,lv,sv,_=boot.background(p,ones,parent,q['edges'],nv,dv,policies[setting,parent])
            saved=stored_responses[ip,case['kind']]
            fit.close(dv-D,saved['denominator_response'],'model_denominator_replay')
            fit.close(lv.sum()-lp.sum(),saved['delta_background_lambda'],'model_Lambda_replay')
            require(sv['broad']==saved['supported'] and sv['lambda_region']==saved['lambda_supported'],'model_support_replay')
            rid=list(stored_responses).index((ip,case['kind']))
            cols=[v[0] for v in active];cells=[v[1] for v in active]
            fit.close((gv-grid).ravel()[cells],response_arrays[rid,cols],'model_native_response_replay')
            r=dict(parent_slot=ip,setting=setting,parent=parent,kind=case['kind'],family='alignment' if case['kind'].startswith('alignment_') else 'delta_influence',
                beta=bv.tolist(),delta_beta=(bv-beta).tolist(),boundary=case['boundary'],support=sv,
                background_broad=float(gv.sum()) if sv['broad'] else None,
                background_lambda_exact=float(lv.sum()) if sv['lambda_region'] else None,
                denominator_lambda_complete=float(dv[native_complete].sum()),regions=regional(q,Xt,bv))
            local.append(len(cases));cases.append(r)
            case_arrays.append((gv,lv,dv,nv,Xt@bv))
        central_grids.append(grid);central_phis.append(lp);central_denoms.append(D);support_grids.append(actual)
        envelopes={}
        for family in ('alignment','delta_influence','all_existing'):
            ids=[k for k in local if family=='all_existing' or cases[k]['family']==family]
            def env(c,ix,support_name=None):
                vals=np.asarray([case_arrays[k][ix] for k in ids]).reshape((len(ids),)+np.shape(c))
                valid=[cases[k]['support'][support_name] if support_name else True for k in ids]
                return supported_envelope(c,vals,valid)
            envelopes[family]=dict(beta=supported_envelope(beta,np.array([cases[k]['beta'] for k in ids]).reshape(len(ids),3),np.ones(len(ids),bool)),
                denominator_MM=env(D,2),numerator_MM=env(num,3),background_MM_phi=env(grid,0,'broad'),background_lambda_phi=env(lp,1,'lambda_region'),
                background_broad=supported_envelope(grid.sum(),np.asarray([case_arrays[k][0].sum() for k in ids]),[cases[k]['support']['broad'] for k in ids]),
                background_lambda_exact=supported_envelope(lp.sum(),np.asarray([case_arrays[k][1].sum() for k in ids]),[cases[k]['support']['lambda_region'] for k in ids]),
                denominator_lambda_complete=supported_envelope(D[native_complete].sum(),np.asarray([case_arrays[k][2][native_complete].sum() for k in ids]),np.ones(len(ids),bool)))
        model_vectors=np.array([[*(np.array(cases[k]['delta_beta'])),cases[k]['denominator_lambda_complete']-D[native_complete].sum(),cases[k]['background_lambda_exact']-lp.sum()] for k in local if cases[k]['support']['lambda_region']])
        summaries.append(dict(setting=setting,parent=parent,central=old['parents'][ip]['central'],bootstrap=old['parents'][ip]['modes'],
            amplitude_decomposition=old['parents'][ip]['amplitude_decomposition'],background_decomposition=old['parents'][ip]['background_decomposition'],
            original_and_fitted_regions={'original':regional(q,X,q['B']),'fitted':central_region},envelopes=envelopes,
            response_co_movement=dict(axes=['B_n','B_delta','B_sidis','control_Lambda_complete','background_Lambda_exact'],cosine=response_cosine(model_vectors) if len(model_vectors) else None,
                meaning='Uncentered normalized response-vector Gram geometry, not a probabilistic model correlation'),
            model_cases=local,prior_widths=q['widths'].tolist(),central_boundary=source['boundaries'],
            fit_support={'occupied_cells':source['occupied_fit_cells'],'rank':source['data_only_rank'],'condition':source['data_design_condition']},
            empty_cell_policy='Empty application cells retain zero conditional resampling support, not physical zero; native control count/MC masks retained in NPZ'))
        checks.append(dict(setting=setting,parent=parent,central_and_model_response_replays=len(local),Lambda_phi_closure=float(lp.sum()-old['parents'][ip]['central']['conditional_lambda_background'])))
        print('PROPAGATED',setting,parent,len(local),'existing cases',flush=True)
    require(len(cases)==len(stored_responses),'model_case_inventory')
    arrays=dict(central_background_MM_phi=np.array(central_grids),central_background_lambda_phi=np.array(central_phis),central_denominator_MM=np.array(central_denoms),
        application_event_counts=np.array(support_grids),model_background_MM_phi=np.array([a[0] for a in case_arrays]),model_background_lambda_phi=np.array([a[1] for a in case_arrays]),
        model_denominator_MM=np.array([a[2] for a in case_arrays]),model_numerator_MM=np.array([a[3] for a in case_arrays]),model_control_delta_MM=np.array([a[4] for a in case_arrays]),
        background_MM_covariance=np.array(mm_cov),background_broad_phi_covariance=np.array(phi_cov),background_lambda_phi_covariance=np.array(phi_lambda_cov),
        denominator_MM_covariance=np.array(denom_cov),denominator_background_MM_cross_covariance=np.array(denom_bg_cross))
    for ip,item in enumerate(order):
        q=fit.build_problem(item['setting'],item['parent']-1,observed,models,roles,docs)
        arrays[f'parent{ip}/observed_count_delta_MM']=q['count'];arrays[f'parent{ip}/observed_sumw2_delta_MM']=q['variance'];arrays[f'parent{ip}/observed_signed_delta_MM']=q['y'];arrays[f'parent{ip}/central_components_delta_MM']=q['Xfull']*np.asarray(lookup[item['setting'],item['parent']]['beta'])
        arrays[f'parent{ip}/finite_MC_template_support']=q['Xfull']>0
    provenance.update(bootstrap_review_sha256=fit.digest(Path(args.bootstrap_review).read_bytes()),propagation_implementation_sha256=fit.digest(Path(__file__).read_bytes()))
    result=dict(schema='kaonlt_pion_uncertainty_propagation/v1',status='DEVELOPMENT COMPLETE, FARM VALIDATION PENDING',parents=summaries,model_cases=cases,axes=interface['axes'],checks=checks,
        source_provenance=provenance,replicas=1024,new_fits=0,model_covariance=None,external_PID_covariance=None,full_physical_transfer_covariance=None,
        missing_terms=old['unavailable_terms'],statistical_scope=old['conditional_transfer'],production_applied=False,scientifically_accepted=False,
        model_scope='Completed existing scenario sensitivity ranges; no probability weights, one-sigma claim, quadrature sum or new calibration prerequisite',
        statistical_evidence='retained_bootstrap.npz is byte-identical to accepted bootstrap NPZ; retained_bootstrap_result/interface JSON preserve decomposition and axes',
        next_gate='Independent actual-diff and scientific review, then existing physical subtraction/final-yield gates with unchanged dependencies')
    encoded=json.dumps(result,indent=2,allow_nan=False)+'\n'
    out.mkdir(parents=True,exist_ok=False);(out/'result.json').write_text(encoded,encoding='utf-8')
    np.savez_compressed(out/'propagation_arrays.npz',**arrays)
    (out/'retained_bootstrap.npz').write_bytes(raw)
    (out/'retained_bootstrap_result.json').write_text(json.dumps(old,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (out/'covariance_interface.json').write_text(json.dumps(interface,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('bootstrap-review','joint-review','fixed-comparison','projection','acquisition','tenbin','output'):p.add_argument('--'+name,type=Path,required=True)
    execute(p.parse_args(argv))


if __name__=='__main__':main()
