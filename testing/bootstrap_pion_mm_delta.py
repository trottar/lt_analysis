"""Source-aware detached pion bootstrap; production and external PID stay frozen."""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import io
import json
from pathlib import Path
import platform
import time
import zipfile
import numpy as np
from scipy.linalg import cho_factor, cho_solve

if __package__:
    from . import fit_pion_mm_delta as fit
    from . import extract_h10_delta as extraction
else:
    import fit_pion_mm_delta as fit
    import extract_h10_delta as extraction

SOURCES = ('prompt', 'rand', 'dummy', 'dummy_rand')
MODES = ('data_only', 'MC_only', 'joint')


def require(ok, why):
    if not ok:
        raise ValueError('pion_bootstrap:'+why)


def identity_groups(source, tree, entry):
    """Retained cohort is supplied by caller; never identify a row by entry alone."""
    require(len(source)==len(tree)==len(entry), 'identity_lengths')
    keys=np.rec.fromarrays([source,np.asarray(tree,dtype='U128'),entry], names='source,tree,entry')
    unique,inverse=np.unique(keys,return_inverse=True)
    return len(unique), inverse


def multiplicities(rng, group):
    n,inverse=group
    return rng.poisson(1., n)[inverse]


def moments(index, weights, multiple, size):
    w=np.asarray(weights,dtype=float);k=np.asarray(multiple,dtype=float)
    require(len(index)==len(w)==len(k) and np.all(k>=0) and np.all(k==np.floor(k)), 'multiplicity')
    return (np.bincount(index,weights=k*w,minlength=size),
            np.bincount(index,weights=k*w*w,minlength=size))


def rebuilt_terms(S,V,L0):
    denominator=float((S*L0[1:-1].sum(0)).sum())
    require(denominator>0 and np.isfinite(denominator), 'empty_MC_parent')
    return fit.mc_terms(S,V,L0/denominator)


def transfer_ratio(numerator, denominator, floor, low, high):
    n=np.asarray(numerator,dtype=float);d=np.asarray(denominator,dtype=float)
    require(n.shape==d.shape and np.isfinite(n).all() and np.isfinite(d).all() and np.all(n>=0) and np.all(d>=0), 'transfer_model')
    unsupported=(d<=floor)&(n>floor)
    w=np.divide(n,d,out=np.zeros_like(n),where=d>floor)
    if low is not None:w=np.maximum(low,w)
    if high is not None:w=np.minimum(high,w)
    w[[0,-1]]=0.;unsupported[[0,-1]]=False
    return w,unsupported


def background(rows, multiple, parent, edges, numerator, denominator, policy):
    ratio,unsupported=transfer_ratio(numerator,denominator,**policy)
    N=len(edges)+1;mb=np.searchsorted(edges,rows['mm'],side='right')
    phi=np.where((rows['pi']>=0)&(rows['pi']<9),rows['pi'],9)
    use=(rows['ti']==parent-1)&rows['nommcuts']
    w=rows['coefficient']*multiple*ratio[mb]
    grid=np.bincount(phi[use]*N+mb[use],weights=w[use],minlength=10*N).reshape(10,N)
    lam=use&rows['allcuts']
    phi_lambda=np.bincount(phi[lam],weights=w[lam],minlength=10)
    valid=dict(broad=not np.any(use&(multiple>0)&unsupported[mb]),
               lambda_region=not np.any(lam&(multiple>0)&unsupported[mb]))
    return grid,phi_lambda,valid,ratio


def covariance(values):
    values=np.asarray(values,dtype=float)
    require(values.ndim==2 and len(values)>=2 and np.isfinite(values).all(), 'covariance_support')
    centered=values-values.mean(0)
    return centered.T@centered/(len(values)-1)


def covariance_summary(values):
    a=np.asarray(values,dtype=float)
    good=np.isfinite(a).all(1);x=a[good]
    result=dict(valid_replicas=int(good.sum()),invalid_replicas=int((~good).sum()),
                common_replica_indices=np.flatnonzero(good).tolist(),covariance=None)
    if len(x)<2:return result
    C=covariance(x);sigma=np.sqrt(np.maximum(0,np.diag(C)))
    corr=[[float(C[i,j]/(sigma[i]*sigma[j])) if sigma[i]>0 and sigma[j]>0 else None
           for j in range(len(sigma))] for i in range(len(sigma))]
    result.update(mean=x.mean(0).tolist(),standard_deviation=sigma.tolist(),covariance=C.tolist(),correlation=corr,
        empirical_percentiles=np.percentile(x,[2.5,16,50,84,97.5],axis=0).tolist(),
        percentile_interpretation='Conditional bootstrap empirical distribution; no physical coverage claim')
    n=len(a)//2;first=a[:n];fg=np.isfinite(first).all(1)
    result['first_half_standard_deviation']=np.sqrt(np.maximum(0,np.diag(covariance(first[fg])))).tolist() if fg.sum()>=2 else None
    # Eight batches expose simulation precision without imposing a convergence cut.
    batch=[]
    for b in np.array_split(x,8):
        if len(b)>=2:batch.append(np.diag(covariance(b)))
    result['batch_variance_MC_standard_error']=(np.std(batch,axis=0,ddof=1)/np.sqrt(len(batch))).tolist() if len(batch)>1 else None
    return result


def paired_decomposition(data,mc,joint,central):
    good=np.isfinite(data).all(1)&np.isfinite(mc).all(1)&np.isfinite(joint).all(1)
    if good.sum()<2:return dict(status='unavailable',common_replicas=int(good.sum()))
    D,M,J=data[good],mc[good],joint[good]
    interaction=J-D-M+np.asarray(central)
    factors=[x-x.mean(0) for x in (D,M,interaction)]
    n=len(D)-1
    components=[x.T@x/n for x in factors]
    cross=sum((factors[i].T@factors[j]+factors[j].T@factors[i])/n for i in range(3) for j in range(i+1,3))
    total=covariance(J)
    require(np.allclose(sum(components)+cross,total,rtol=1e-10,atol=1e-12),'paired_covariance_identity')
    return dict(status='available_conditional',common_replicas=int(good.sum()),
        covariance_data=components[0].tolist(),covariance_MC=components[1].tolist(),
        covariance_nonlinear_interaction=components[2].tolist(),covariance_paired_cross=cross.tolist(),
        covariance_joint=total.tolist(),closure_max_error=float(np.max(abs(sum(components)+cross-total))))


def policy_from_record(value):
    require('clip_max' in value and 'clip_min' in value and 'denom_floor' in value, 'missing_transfer_policy')
    maximum=value['clip_max']
    if isinstance(maximum,dict):
        # Owner always emits clip_max; snapshot(None) produces exactly this marker.
        require(maximum=={'status':'unavailable','reason':'not retained by producer'},'unknown_clip_policy')
        maximum=None
    return dict(floor=float(value['denom_floor']),low=value['clip_min'],high=maximum)


def load_evidence(joint_path,fixed_path,projection_path,acquisition_path,tenbin_path):
    observed,models,roles,docs,provenance=fit.load_inputs(fixed_path,projection_path)
    with zipfile.ZipFile(joint_path) as z:
        manifest=json.loads(z.read('members-sha256.json'))
        for n,m in manifest.items():
            b=z.read(n);require((len(b),fit.digest(b))==(m['bytes'],m['sha256']),'joint_member_hash')
        require(z.read('after/testing/fit_pion_mm_delta.py')==Path(fit.__file__).read_bytes(),'accepted_joint_source_changed')
        accepted=json.loads(z.read('result.json'))
    require(accepted['full_scope'] and not accepted['production_applied'],'joint_scope')
    require(fit.digest(Path(acquisition_path).read_bytes())==provenance['historical_reference_sha256'],'acquisition_hash')
    with zipfile.ZipFile(fixed_path) as z:
        expected=json.loads(z.read('verification.json'))['observed_tenbin_archive']['sha256']
    require(fit.digest(Path(tenbin_path).read_bytes())==expected,'tenbin_hash')
    with zipfile.ZipFile(tenbin_path) as z:
        allrows=fit.arrays(z.read('physical_pion_joint_rows.npz'))
    data={};original={};policies={};A={};identity={}
    with zipfile.ZipFile(acquisition_path) as z:
        for setting in fit.SETTINGS:
            name,eps=setting.split('_');stem=f'global/{name}_kaon_Q4p4W2p74_{eps}'
            meta=json.loads(z.read(stem+'_check1_data.json'));raw=z.read(stem+'_check1_data.npz')
            require(fit.digest(raw)==meta['npz']['sha256'],'event_npz_hash')
            table={k.split('/')[-1]:v for k,v in fit.arrays(raw).items() if k.startswith('rows/pion/')}
            p={k.split('/',1)[1]:v for k,v in allrows.items() if k.startswith(setting+'/')}
            index={(int(s),str(t),int(e)):i for i,(s,t,e) in enumerate(zip(table['source'],table['tree'],table['entry']))}
            order=np.array([index[(int(s),str(t),int(e))] for s,t,e in zip(p['source'],p['tree'],p['entry'])])
            require(len(np.unique(order))==len(order)==len(index),'event_binding')
            for field in ('source','entry','tree','mm','ti','pi','coefficient','weight','allcuts','nommcuts'):
                require(np.array_equal(p[field],table[field][order]),'event_field_binding')
            di=np.where((p['delta_index_exact']>=0)&(p['delta_index_exact']<10),p['delta_index_exact'],10)
            p['di']=di;data[setting]=p
            original[setting]=json.loads(z.read(stem+'_pion-t-template-comparison.json'))
            child=json.loads(z.read(stem+'_pion-t-phi-template-comparison.json'))
            identity[setting]=dict(source_bindings_sha256=meta['source_bindings_sha256'],row_count=len(order),
                group_key='retained cohort binding + source + PID tree + source entry',
                cross_PID_physical_event_map='unavailable; PID trees have local entry counters, no common run/event keys')
            for j,parent in enumerate(original[setting]['parents']):
                f=parent['candidate_B'];dx=f['diagnostics']['kaon'];protected=dx['pi_delta_signal_protected_fit']
                amps=protected['protected_applied_amplitudes'] if protected['enabled'] else dx['refined_amplitudes'] or dx['staged_amplitudes_scaled']
                A[setting,j+1]=np.array([amps[c] for c in fit.CHANNELS])
                available=[c['applications']['Local_t'] for c in child['cells'] if c['t_index']==j and c['applications']['Local_t']['status']=='available']
                require(bool(available),'missing_transfer')
                require(all(a['weights']==available[0]['weights'] for a in available),'inconsistent_parent_transfer')
                policies[setting,j+1]=policy_from_record(available[0]['policy'])
                original[setting]['parents'][j]['_weights']=available[0]['weights']
    provenance.update(joint_review_sha256=fit.digest(Path(joint_path).read_bytes()),tenbin_sha256=expected,
        bootstrap_implementation_sha256=fit.digest(Path(__file__).read_bytes()),experimental_identity=identity)
    return observed,models,roles,docs,accepted,data,original,policies,A,provenance


def fit_replica(problem,y,variance,terms):
    N=y.shape[1]
    domain=np.zeros(N,bool);domain[np.unique(problem['take']%N)]=True
    take=np.flatnonzero(((variance>0)&domain)[0:10])
    require(len(take)>3,'no_replica_fit_support')
    Xfull=np.stack([t[2] for t in terms],axis=-1)
    C=np.diag(variance.ravel()[take])+fit.mc_covariance(terms,problem['B'],take)
    r=fit.solve(Xfull.reshape(-1,3)[take],y.ravel()[take],C,problem['target'],problem['widths'])
    return r,Xfull


def execute(args):
    output=Path(args.output);require(not output.exists(),'output_exists')
    require(args.replicas>=2,'replicas')
    evidence=load_evidence(args.joint_review,args.fixed_comparison,args.projection,args.acquisition,args.tenbin)
    observed,models,roles,docs,accepted,data,original,policies,A,provenance=evidence
    require(set(args.settings)<=set(fit.SETTINGS) and len(set(args.settings))==len(args.settings),'settings')
    ns=len(args.settings);R=args.replicas;P=ns*3
    params={m:np.full((R,P,3),np.nan) for m in MODES}
    control={m:np.full((R,P),np.nan) for m in MODES}
    bg_lambda={m:np.full((R,P),np.nan) for m in MODES}
    phi_lambda={m:np.full((R,P,10),np.nan) for m in MODES}
    bounds={m:np.zeros((R,P,3),bool) for m in MODES}
    denom_joint=np.full((R,P,102),np.nan)
    control_delta_joint=np.full((R,P,10),np.nan)
    residual_delta_joint=np.full((R,P,10),np.nan)
    residual_delta_fixed=np.full((R,P,10),np.nan)
    failures=[];closures=[];axes=[];central=[];background_vectors=[];model_responses=[];discrepancy=[]
    order=[]
    for setting in args.settings:
        for j in (3,1,2):order.append((setting,j))
    accepted_lookup={(x['setting'],x['parent']):x for x in accepted['parents']}
    # Store one factor column per actual application cell; empty cells remain explicit in axes.
    native_rows=[];active_map={};offset=0
    problems={}
    for ip,(setting,parent) in enumerate(order):
        q=fit.build_problem(setting,parent-1,observed,models,roles,docs);problems[setting,parent]=q;N=q['y'].shape[1]
        require(N==102,'native_MM_scope')
        b=np.array(accepted_lookup[setting,parent]['beta']);p=data[setting];ones=np.ones(len(p['mm']))
        X=q['Xfull'];D=X.sum(0)@b;num=X.sum(0)@A[setting,parent]
        grid,lphi,valid,ratio=background(p,ones,parent,q['edges'],num,D,policies[setting,parent])
        require(valid['broad'] and valid['lambda_region'],'central_transfer_unsupported')
        replay,_=fit_replica(q,q['y'],q['variance'],q['terms'])
        fit.close(replay['beta'],b,'central_joint_refit_replay')
        oldD=X.sum(0)@q['B'];_,_,oldvalid,oldratio=background(p,ones,parent,q['edges'],num,oldD,policies[setting,parent])
        fit.close(oldratio,original[setting]['parents'][parent-1]['_weights'],'original_transfer_replay')
        mb=np.searchsorted(q['edges'],p['mm'],side='right')
        use=(p['ti']==parent-1)&p['nommcuts'];fit.close(p['weight'][use],p['coefficient'][use]*oldratio[mb[use]],'original_event_transfer')
        phi=np.where((p['pi']>=0)&(p['pi']<9),p['pi'],9);active=np.unique(phi[use]*N+mb[use])
        active_map[setting,parent]=(active,offset);offset+=len(active)
        for index in range(10*N):
            ph,mm=divmod(index,N);native_rows.append(dict(parent_slot=ip,setting=setting,t_bin=parent,phi_bin=ph,
                mm_bin=mm,application_event_support=bool(index in active),factor_column=int(np.searchsorted(active,index)+active_map[setting,parent][1]) if index in active else None))
        normal=np.r_[False,np.ones(N-2,bool),False];lam=np.r_[False,(q['edges'][:-1]>=1.1)&(q['edges'][1:]<=1.16),False]
        central.append(dict(setting=setting,parent=parent,beta=b.tolist(),lambda_control=float(D[lam].sum()),
            conditional_lambda_background=float(lphi.sum()),lambda_phi=lphi.tolist(),denominator=D.tolist(),
            original_lambda_background=float(np.sum(p['weight'][(p['ti']==parent-1)&p['allcuts']]))))
        closures.append(dict(setting=setting,parent=parent,original_transfer_max_error=float(np.max(abs(oldratio-original[setting]['parents'][parent-1]['_weights']))),
            accepted_beta_replay_error=q['original_replay_max_error'],rows=int(use.sum())))
        # Correlated model responses from accepted, source-defined variants; no invented law.
        for variant in accepted_lookup[setting,parent]['sensitivities']:
            if not (variant['kind'].startswith('alignment_') or variant['kind'].startswith('omit_delta_')):continue
            Xt=X.copy()
            if variant['kind'].startswith('alignment_'):
                for ci,ch in enumerate(fit.CHANNELS):
                    lead='alignment_'+ch+'_'
                    if variant['kind'].startswith(lead):
                        ref=deepcopy(q['parent_refs'][ci]);ref['shift']=float(variant['kind'][len(lead):]);L=extraction.transform(ref)
                        Xt[:,:,ci]=q['raw_parts'][ci]@L.T
            bv=np.asarray(variant['beta']);Dv=Xt.sum(0)@bv;nv=Xt.sum(0)@A[setting,parent]
            gv,lv,vvalid,_=background(p,ones,parent,q['edges'],nv,Dv,policies[setting,parent])
            model_responses.append(dict(parent_slot=ip,kind=variant['kind'],supported=vvalid['broad'],lambda_supported=vvalid['lambda_region'],delta_beta=(bv-b).tolist(),
                delta_control_lambda=float(Dv[lam].sum()-D[lam].sum()),delta_background_lambda=float(lv.sum()-lphi.sum()),
                background_response_active=(gv-grid).ravel()[active].tolist(),denominator_response=(Dv-D).tolist()))
        empirical=q['y'][:10,normal].sum(1);model=(X@b)[:10,normal].sum(1)
        discrepancy.append(dict(setting=setting,parent=parent,observed_delta=empirical.tolist(),central_delta=model.tolist(),
            residual_bias_direction=(empirical-model).tolist(),nuisance_covariance=None,
            interpretation='Measured correlated bias direction; not calibrated model variance and not removed from central residual'))
    background_joint=np.full((R,offset),np.nan)
    for iset,setting in enumerate(args.settings):
        p=data[setting];N=102;ps=[problems[setting,j] for j in (3,1,2)]
        group=identity_groups(p['source'],p['tree'],p['entry'])
        require(group[0]==len(p['entry']),'duplicate_population_rows')
        seed=lambda role:int.from_bytes(hashlib.sha256(f'{args.seed}|{setting}|{role}'.encode()).digest()[:8],'big')
        rngD=np.random.default_rng(seed('experimental'))
        rngM=[np.random.default_rng(seed('SIMC_'+ch)) for ch in fit.CHANNELS]
        # One index/multiplicity for all canonical parents and selections in a cohort.
        pmb=np.searchsorted(ps[0]['edges'],p['mm'],side='right')
        accepted_p=(p['ti']>=0)&(p['ti']<3)&p['nommcuts']
        pi_index=(p['ti'][accepted_p]*11+p['di'][accepted_p])*N+pmb[accepted_p]
        mc=[]
        for ch in fit.CHANNELS:
            a=models[setting+'_'+ch];valid=(a['t_bin']>=0)&(a['t_bin']<3)&(a['phi_bin']>=0)
            di=np.where(a['delta_bin']<0,10,a['delta_bin'])
            index=(a['t_bin'][valid]*11+di[valid])*N+a['mm_bin'][valid]
            require(len(np.unique(a['entry']))==len(a['entry']),'duplicate_SIMC_entry')
            phi_index=(a['t_bin'][valid]*9+a['phi_bin'][valid])*N+a['mm_bin'][valid]
            mc.append((a,valid,index,phi_index))
        start=time.monotonic()
        for rep in range(R):
            kd=multiplicities(rngD,group)
            yd,vd=moments(pi_index,p['coefficient'][accepted_p],kd[accepted_p],3*11*N)
            yd=yd.reshape(3,11,N);vd=vd.reshape(3,11,N)
            mc_arrays=[]
            for ci,(a,valid,index,phi_index) in enumerate(mc):
                km=rngM[ci].poisson(1.,len(a['entry']))
                S,V=moments(index,a['weight'][valid],km[valid],3*11*N)
                SP,VP=moments(phi_index,a['weight'][valid],km[valid],3*9*N)
                S=S.reshape(3,11,N);V=V.reshape(3,11,N)
                fit.close(S.sum(1),SP.reshape(3,9,N).sum(1),'bootstrap_phi_parent_content')
                fit.close(V.sum(1),VP.reshape(3,9,N).sum(1),'bootstrap_phi_parent_Sumw2')
                mc_arrays.append((S,V))
            for loc,q in enumerate(ps):
                parent=q['parent'];j=parent-1;ip=iset*3+loc
                try:
                    tmc=[rebuilt_terms(mc_arrays[c][0][j],mc_arrays[c][1][j],q['terms'][c][0]) for c in range(3)]
                except (ValueError,np.linalg.LinAlgError) as e:
                    failures.append(dict(setting=setting,parent=parent,replica=rep,mode='MC_only_and_joint',reason=str(e)));tmc=None
                for mode in MODES:
                    try:
                        terms=q['terms'] if mode=='data_only' else tmc
                        require(terms is not None,'failed_MC_parent')
                        Y=q['y'] if mode=='MC_only' else yd[j];V=q['variance'] if mode=='MC_only' else vd[j]
                        ans,X=fit_replica(q,Y,V,terms);b=ans['beta'];D=X.sum(0)@b;num=X.sum(0)@A[setting,parent]
                        lam=np.r_[False,(q['edges'][:-1]>=1.1)&(q['edges'][1:]<=1.16),False]
                        params[mode][rep,ip]=b;bounds[mode][rep,ip]=ans['active'];control[mode][rep,ip]=D[lam].sum()
                        k= np.ones(len(kd)) if mode=='MC_only' else kd
                        grid,lphi,supported,_=background(p,k,parent,q['edges'],num,D,policies[setting,parent])
                        if supported['lambda_region']:
                            bg_lambda[mode][rep,ip]=lphi.sum();phi_lambda[mode][rep,ip]=lphi
                        else:failures.append(dict(setting=setting,parent=parent,replica=rep,mode=mode,reason='unsupported_Lambda_transfer_on_resampled_event'))
                        if supported['broad']:
                            if mode=='joint':
                                active,offset0=active_map[setting,parent]
                                background_joint[rep,offset0:offset0+len(active)]=grid.ravel()[active]
                        else:failures.append(dict(setting=setting,parent=parent,replica=rep,mode=mode,reason='unsupported_broad_transfer_on_resampled_event'))
                        if mode=='joint':
                            denom_joint[rep,ip]=D
                            control_delta_joint[rep,ip]=(X@b)[:10,1:-1].sum(1)
                            residual_delta_joint[rep,ip]=(Y-X@b)[:10,1:-1].sum(1)
                            residual_delta_fixed[rep,ip]=(Y-X@q['B'])[:10,1:-1].sum(1)
                    except (ValueError,np.linalg.LinAlgError) as e:
                        failures.append(dict(setting=setting,parent=parent,replica=rep,mode=mode,reason=str(e)))
            if (rep+1)%64==0:print('BOOTSTRAP',setting,rep+1,'/',R,'elapsed_s',round(time.monotonic()-start,1),flush=True)
    summaries=[]
    for ip,(setting,parent) in enumerate(order):
        modes={}
        for mode in MODES:
            c=covariance_summary(params[mode][:,ip]);good=np.isfinite(params[mode][:,ip]).all(1)
            c['boundary_hit_frequency']=bounds[mode][good,ip].mean(0).tolist() if good.any() else None
            c['lambda_control']=covariance_summary(control[mode][:,ip:ip+1])
            c['conditional_lambda_background']=covariance_summary(bg_lambda[mode][:,ip:ip+1])
            modes[mode]=c
        summaries.append(dict(setting=setting,parent=parent,central=central[ip],modes=modes,
            amplitude_decomposition=paired_decomposition(params['data_only'][:,ip],params['MC_only'][:,ip],params['joint'][:,ip],central[ip]['beta']),
            background_decomposition=paired_decomposition(bg_lambda['data_only'][:,ip:ip+1],bg_lambda['MC_only'][:,ip:ip+1],bg_lambda['joint'][:,ip:ip+1],[central[ip]['conditional_lambda_background']])))
    good=np.isfinite(background_joint).all(1)
    factor=(background_joint[good]-background_joint[good].mean(0))/np.sqrt(good.sum()-1) if good.sum()>=2 else np.empty((0,offset))
    arrays=dict(background_joint_replicas=background_joint,background_covariance_factor=factor,background_supported_replica_indices=np.flatnonzero(good),denominator_joint_replicas=denom_joint,
                control_delta_joint_replicas=control_delta_joint,residual_delta_joint_replicas=residual_delta_joint,residual_delta_fixed_replicas=residual_delta_fixed)
    for mode in MODES:
        arrays.update({mode+'/amplitudes':params[mode],mode+'/lambda_control':control[mode],mode+'/lambda_background':bg_lambda[mode],mode+'/lambda_phi_background':phi_lambda[mode],mode+'/boundaries':bounds[mode]})
    model_factors=np.zeros((len(model_responses),offset))
    for i,row in enumerate(model_responses):
        ip=row['parent_slot'];setting,parent=order[ip];active,offset0=active_map[setting,parent]
        model_factors[i,offset0:offset0+len(active)]=row.pop('background_response_active')
    arrays['model_background_response_vectors']=model_factors
    axes=dict(parent_order=[dict(setting=s,parent=j) for s,j in order],native_background_cells=native_rows,
              background_factor_order='All MM/phi cells with retained application-event support; unrepresented cells have conditional resampling zero, not known physical zero',
              phi_bins='0..8 canonical; 9 explicit out-of-range flow',MM_edges=ps[0]['edges'].tolist(),
              covariance_reconstruction='C_stat = factor.T @ factor; model response D requires supplied nuisance covariance Sigma_eta: D.T @ Sigma_eta @ D')
    payload=dict(schema='kaonlt_pion_source_bootstrap/v1',status='DEVELOPMENT COMPLETE, FARM VALIDATION PENDING',source_provenance=provenance,
        versions=dict(python=platform.python_version(),numpy=np.__version__),replicas=R,seed=args.seed,
        full_scope=set(args.settings)==set(fit.SETTINGS),resampling='Poisson(1) per retained source/tree/entry and SIMC file/entry; same multiplier across linked uses',
        precision_policy='Recomputed per replica: counting and correlated finite-MC at original B; staged targets/widths and alignment conditional/fixed',
        parents=summaries,closure=closures,failures=failures,model_responses=model_responses,model_discrepancy=discrepancy,
        bootstrap_parent_phi_closure_checks=R*P*3,
        global_conditional_lambda_covariance=covariance_summary(bg_lambda['joint']),global_conditional_phi_lambda_covariance=covariance_summary(phi_lambda['joint'].reshape(R,-1)),
        global_conditional_amplitude_covariance=covariance_summary(params['joint'].reshape(R,-1)),
        full_transfer_covariance=None,model_covariance=None,external_PID_covariance=None,
        unavailable_terms=['Protected kaon-side numerator coefficient refit and its learned/training covariance','Common raw-event map across PID trees and raw upstream PID/proton-learning populations',
                           'Re-estimated staged targets and widths/alignment learning','Physical model-nuisance probability law and covariance','External PID systematic and cross covariance',
                           'Unknown shared normalization/cross-setting systematic covariance; block-independent source draws do not estimate it'],
        conditional_transfer='Fixed protected A, shared resampled SIMC numerator/denominator, shared pion data fit/application; not total physical uncertainty',
        no_independence_assumed_for_missing_terms=True,production_applied=False,scientifically_accepted=False)
    interface=dict(schema='kaonlt_covariance_combination_interface/v1',axes=axes,units='existing signed-yield units; covariance squared yield',
        statistical_covariance={'status':'available_conditional' if good.sum()>=2 else 'unavailable_insufficient_common_replicas','factor_npz_key':'background_covariance_factor','replica_indices_key':'background_supported_replica_indices'},
        model={'status':'response_vectors_only_uncalibrated','response_npz_key':'model_background_response_vectors','row_order_and_validity':'result.json model_responses; supported/lambda_supported flags required','nuisance_covariance':None,'covariance':None},
        PID={'status':'unavailable_external_owner','covariance':None,'cross_covariance_statistical_model':None,'assume_independent':False},
        full_transfer={'status':'unavailable_missing_numerator_and_learning_covariance','covariance':None},
        combination_rule='Verify axes, physical populations, PID definition/overlap and cross-covariances before summing; never substitute null with zero')
    encoded=json.dumps(payload,indent=2,allow_nan=False)+'\n';output.mkdir(parents=True,exist_ok=False)
    (output/'result.json').write_text(encoded,encoding='utf-8');(output/'covariance_interface.json').write_text(json.dumps(interface,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    np.savez_compressed(output/'replicas_and_covariance_factors.npz',**arrays)
    return payload


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('joint-review','fixed-comparison','projection','acquisition','tenbin','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--replicas',type=int,default=1024);p.add_argument('--seed',type=int,default=20261111)
    p.add_argument('--settings',nargs='+',choices=fit.SETTINGS,default=list(fit.SETTINGS))
    execute(p.parse_args(argv))


if __name__=='__main__':main()
