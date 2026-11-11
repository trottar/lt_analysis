"""Deterministic tests for detached joint MM/delta numerical ownership."""
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import zipfile

import numpy as np
from scipy.optimize import nnls

if __package__:
    from . import fit_pion_mm_delta as fit
else:
    import fit_pion_mm_delta as fit


class JointNumerics(unittest.TestCase):
    def setUp(self):
        self.X = np.array([[1.,0,0],[0,1,0],[0,0,1],[1,2,3],[3,2,1]])
        self.B = np.array([2.,3.,4.])
        self.C = np.eye(5)

    def solve(self, y=None, X=None, C=None, target=None, widths=None):
        return fit.solve(self.X if X is None else X, self.X@self.B if y is None else y,
                         self.C if C is None else C, self.B if target is None else target,
                         np.ones(3) if widths is None else widths)

    def test_exact_shared_coefficients(self):
        np.testing.assert_allclose(self.solve()['beta'], self.B, atol=1e-12)

    def test_same_three_coefficients_across_replicated_delta(self):
        X = np.tile(self.X, (10,1))/10
        r = self.solve(X=X, y=X@self.B, C=np.eye(50))
        self.assertEqual(r['beta'].shape, (3,))
        np.testing.assert_allclose(r['beta'], self.B, atol=1e-12)

    def test_signed_negative_observation_retained(self):
        y = self.X@self.B
        y[0] = -2
        a = self.solve(y=y)
        self.assertGreater(a['data_quadratic'], 0)
        self.assertLess(a['beta'][0], self.B[0])

    def test_nonnegative_boundary(self):
        r = fit.solve(np.eye(3), np.array([-3.,2.,1.]), np.eye(3), np.zeros(3), np.full(3,np.inf))
        self.assertTrue(r['active'][0])
        self.assertEqual(r['beta'][0], 0.)
        p = fit.profiles(r)
        self.assertEqual(p[0]['lower'], 0.)
        self.assertGreater(p[0]['upper'], 0)

    def test_all_bound_is_reportable(self):
        X = np.eye(3)
        r = fit.solve(X, -np.ones(3), np.eye(3), np.zeros(3), np.full(3,np.inf))
        c = fit.sampling_covariance(X, np.eye(3), np.ones(3), r)
        np.testing.assert_array_equal(c['total'], np.zeros((3,3)))
        self.assertTrue(all(p['lower_is_boundary'] for p in fit.profiles(r)))

    def test_matches_independent_nnls(self):
        y = self.X@self.B+np.array([1,-1,.5,-.1,.3])
        r = self.solve(y=y)
        expected = nnls(np.vstack([self.X,np.eye(3)]),np.r_[y,self.B])[0]
        np.testing.assert_allclose(r['beta'], expected, rtol=1e-10)

    def test_correlated_covariance_whitening(self):
        C = np.eye(5)+.2*np.ones((5,5))
        y = self.X@self.B+np.array([.3,-.7,.5,1,-1])
        r = self.solve(C=C,y=y)
        expected = np.linalg.solve(self.X.T@np.linalg.solve(C,self.X)+np.eye(3), self.X.T@np.linalg.solve(C,y)+self.B)
        np.testing.assert_allclose(r['beta'], expected, rtol=1e-10)

    def test_scale_equivariance(self):
        y = self.X@self.B+np.array([1,-1,.5,-.1,.3])
        r = self.solve(y=y)
        f = 1000
        s = self.solve(y=y*f,C=self.C*f*f,target=self.B*f,widths=np.ones(3)*f)
        np.testing.assert_allclose(s['beta'], r['beta']*f, rtol=1e-10)

    def test_no_new_regularization(self):
        r = self.solve(widths=np.full(3,np.inf))
        self.assertEqual(len(r['prior_target']),0)

    def test_nan_fails(self):
        with self.assertRaisesRegex(ValueError,'nonfinite'):
            self.solve(y=np.full(5,np.nan))

    def test_shape_fails(self):
        with self.assertRaisesRegex(ValueError,'shape'):
            self.solve(C=np.eye(4))

    def test_nonpositive_prior_fails(self):
        with self.assertRaisesRegex(ValueError,'prior_widths'):
            self.solve(widths=np.array([1.,0,1]))

    def test_nonsymmetric_covariance_fails(self):
        C = self.C.copy();C[0,1]=.5
        with self.assertRaisesRegex(ValueError,'symmetry'):
            self.solve(C=C)

    def test_indefinite_covariance_fails(self):
        C = self.C.copy();C[0,0]=-1
        with self.assertRaises(np.linalg.LinAlgError):
            self.solve(C=C)

    def test_empty_fit_fails(self):
        with self.assertRaises(ValueError):
            fit.solve(np.zeros((0,3)),np.zeros(0),np.zeros((0,0)),self.B,np.ones(3))

    def test_nonconverged_solver_fails(self):
        with mock.patch.object(fit, 'lsq_linear', return_value=mock.Mock(success=False)):
            with self.assertRaisesRegex(ValueError,'nonconvergence'):
                self.solve()

    def test_data_MC_covariance_decomposes(self):
        C = self.C+.2*np.ones((5,5))
        r = self.solve(C=C)
        c = fit.sampling_covariance(self.X,C,np.ones(5),r)
        np.testing.assert_allclose(c['total'],c['data']+c['mc'],atol=1e-14)
        self.assertGreater(np.linalg.eigvalsh(c['total']).min(),0)
        np.testing.assert_allclose(np.diag(c['correlation']),1)

    def test_no_prior_sampling_equals_curvature(self):
        r = self.solve(widths=np.full(3,np.inf))
        c = fit.sampling_covariance(self.X,self.C,np.ones(5),r)
        np.testing.assert_allclose(c['total'],c['curvature'],atol=1e-14)

    def test_same_data_prior_is_not_random_observation(self):
        r = self.solve(widths=np.ones(3)*.1)
        c = fit.sampling_covariance(self.X,self.C,np.ones(5),r)
        self.assertTrue(np.all(np.diag(c['total']) < np.diag(c['curvature'])))

    def test_profile_exact_quadratic(self):
        r = self.solve()
        c = fit.sampling_covariance(self.X,self.C,np.ones(5),r)
        p = fit.profiles(r)
        for j in range(3):
            self.assertAlmostEqual((p[j]['upper']-self.B[j])**2,c['curvature'][j,j],places=10)


class NormalizationCovariance(unittest.TestCase):
    def setUp(self):
        self.S = np.array([[0.,3.,7.,0.],[0.,2.,8.,0.]])
        self.V = self.S*2  # Deliberately unequal raw content and Sumw2.
        self.L = np.eye(4)/20
        self.terms = [fit.mc_terms(self.S,self.V,self.L)]

    def test_dense_J_matches_block_formula(self):
        take=np.arange(8);L,V,Q,u,a,g=self.terms[0]
        J=np.kron(np.eye(2),L)-np.outer(Q.ravel(),np.tile(g,2))
        expected=(J*V.ravel())@J.T
        C=fit.mc_covariance(self.terms,[1.],take)
        np.testing.assert_allclose(C,expected,atol=1e-15)
        self.assertGreater(np.linalg.eigvalsh(C).min(),-1e-14)

    def test_normalized_total_variance_zero(self):
        mask=np.ones((2,4));mask[:,[0,-1]]=False
        self.assertAlmostEqual(fit.integrated_mc_variance(self.terms,[1.],mask),0,places=14)

    def test_delta_slices_correlated(self):
        C=fit.mc_covariance(self.terms,[1.],np.arange(8))
        self.assertLess(C[1:3,5:7].sum(),0)

    def test_scalar_functional_matches_dense(self):
        C=fit.mc_covariance(self.terms,[2.],np.arange(8))
        mask=np.zeros((2,4));mask[0,1:3]=1
        self.assertAlmostEqual(fit.integrated_mc_variance(self.terms,[2.],mask),float(mask.ravel()@C@mask.ravel()),places=14)

    def test_bad_normalization_fails(self):
        with self.assertRaisesRegex(ValueError,'normalization_identity'):
            fit.mc_terms(self.S,self.V,np.eye(4))

    def test_variance_independent_of_content(self):
        other=[fit.mc_terms(self.S,self.V*4,self.L)]
        np.testing.assert_allclose(fit.mc_covariance(other,[1.],np.arange(8)),4*fit.mc_covariance(self.terms,[1.],np.arange(8)),atol=1e-15)

    def test_template_parameter_derivative_with_model_residual(self):
        rng=np.random.default_rng(92)
        raw=[];terms=[]
        for j in range(3):
            S=np.zeros((2,6));S[:,1:-1]=rng.uniform(5,30,(2,4));raw.append(S)
            terms.append(fit.mc_terms(S,S,np.eye(6)/S[:,1:-1].sum()))
        Xfull=np.stack([t[2] for t in terms],axis=-1)
        take=np.array([1,2,3,4,7,8,9,10]);X=Xfull.reshape(-1,3)[take]
        target=np.array([2.,3.,4.]);y=X@target+np.array([.03,-.05,.04,.02,-.01,.02,-.03,.01])
        C=np.eye(8)*.01+fit.mc_covariance(terms,target,take)
        a=fit.solve(X,y,C,target,np.ones(3));c=fit.sampling_covariance(X,C,np.full(8,.01),a)
        q=dict(take=take,y=np.zeros((2,6)),terms=terms)
        derivatives,mc=fit.template_parameter_derivatives(q,X,y,C,a,c)
        self.assertGreater(np.linalg.eigvalsh(mc).min(),-1e-14)
        for ci in range(3):
            h=1e-3;predictions=[]
            for sign in (-1,1):
                S=raw[ci].copy();S[0,2]+=sign*h
                Q=S/S[:,1:-1].sum()
                xx=X.copy();xx[:,ci]=Q.ravel()[take]
                predictions.append(fit.solve(xx,y,C,target,np.ones(3))['beta'])
            np.testing.assert_allclose((predictions[1]-predictions[0])/(2*h),derivatives[ci][:,2],rtol=1e-5,atol=1e-9)


class DetachedOwnership(unittest.TestCase):
    def problem_fixture(self):
        setting='Center_highe';prefix=setting+'/t1';N=6
        edges=np.linspace(.7,1.3,5);B=np.array([2.,3.,4.])
        shapes=np.array([[0,.6,.2,.1,.1,0],[0,.1,.6,.2,.1,0],[0,.1,.1,.3,.5,0]]).T
        y=np.zeros((11,N));y[0]=shapes@B;count=(y>0).astype(float);variance=count.copy()
        observed={prefix+'/B':B,prefix+'/mm_edges':edges}
        for field,a in [('data',y),('sumw2',variance),('count',count)]:
            observed[prefix+'/nommcuts/'+field]=a
            for source in ('prompt','rand','dummy','dummy_rand'):
                observed[f'{prefix}/nommcuts/{source}/{field}']=a if source=='prompt' else np.zeros_like(a)
        models={};roles={}
        for ci,ch in enumerate(fit.CHANNELS):
            S=np.zeros((11,N));S[0]=shapes[:,ci]*1000
            if ci==0:S[0]=np.array([0,400,400,100,100,0])  # Legitimately regenerated shape.
            L=np.eye(N)/1000;fields=dict(raw_sumw=S,raw_sumw2=S,count=S,alignment_operator=L,aligned=S@L.T)
            models[setting+'_'+ch]={'t1_'+k:v for k,v in fields.items()}
            for k,v in fields.items():observed[f'{prefix}/{ch}/{k}']=v
            roles[setting+'_'+ch]={'closure_reference':{'parents':[{'renormalize':True}]}}
        joint=dict(coefficients=dict(zip(fit.CHANNELS,B)),regularization_widths=dict.fromkeys(fit.CHANNELS,1.),regularization_enabled=True,fit_bin_indices=[1,2,3,4],data_chi2_contribution=0.)
        doc={'pion_diagnostics':dict(joint_refinement=joint,staged_component_amplitudes=dict(zip(fit.CHANNELS,B))),
             'aligned_shapes':{ch:dict(values=shapes[:,ci],edges=edges) for ci,ch in enumerate(fit.CHANNELS)},
             'frozen_inputs':{'pion':dict(values=y.sum(0),stored_sumw2=variance.sum(0))},'pion_component_alignment':{}}
        return setting,observed,models,roles,{setting:{'parents':[doc]}}

    def test_regenerated_shapes_need_not_reproduce_historical_fit(self):
        s,o,m,r,d=self.problem_fixture()
        q=fit.build_problem(s,0,o,m,r,d)
        np.testing.assert_allclose(q['Xfull'][0,:,0],[0,.4,.4,.1,.1,0])
        self.assertLess(q['original_replay_max_error'],1e-10)

    def test_changed_fixed_projection_binding_fails(self):
        s,o,m,r,d=self.problem_fixture();o[s+'/t1/pi_n/aligned']=o[s+'/t1/pi_n/aligned']+.1
        with self.assertRaisesRegex(ValueError,'fixed_projection_binding'):
            fit.build_problem(s,0,o,m,r,d)

    def test_negative_counting_variance_fails(self):
        s,o,m,r,d=self.problem_fixture();o[s+'/t1/nommcuts/sumw2'][0,1]=-1
        with self.assertRaisesRegex(ValueError,'negative_array'):
            fit.build_problem(s,0,o,m,r,d)

    def test_empty_data_not_assigned_variance(self):
        s,o,m,r,d=self.problem_fixture();o[s+'/t1/nommcuts/count'][0,1]=0
        with self.assertRaisesRegex(ValueError,'empty_variance_classification'):
            fit.build_problem(s,0,o,m,r,d)

    def test_no_ROOT_or_production_imports(self):
        import ast
        tree=ast.parse(Path(fit.__file__).read_text(encoding='utf-8'))
        names=[n.name for n in ast.walk(tree) if isinstance(n,ast.alias)]
        self.assertFalse(any(n.startswith(('ROOT','main','pion_component_fits')) for n in names))

    def test_existing_output_refused_before_loading(self):
        with tempfile.TemporaryDirectory() as d:
            with mock.patch.object(fit,'load_inputs') as load:
                with self.assertRaisesRegex(ValueError,'output_exists'):
                    fit.run('missing','missing',d)
                load.assert_not_called()

    def archive_pair(self, folder, *, wrong_projection=False, bad_status=False, bad_reference=False, bad_member=False):
        p=Path(folder)/'projection.zip';f=Path(folder)/'fixed.zip'
        receipt=dict(status='bad' if bad_status else 'verified_projection', parent_channel_closures=45,
                     primary_reference_sha256='wrong' if bad_reference else 'original',runtime_head='head',delta_edges=fit.DELTA.tolist(),roles=[])
        with zipfile.ZipFile(p,'w') as z:z.writestr('receipt.json',json.dumps(receipt))
        verification=dict(projection={'sha256':'wrong' if wrong_projection else fit.digest(p.read_bytes())},historical_acquisition={'sha256':'original'},head='head')
        b=io.BytesIO();np.savez_compressed(b)
        members={'verification.json':json.dumps(verification).encode(),'numerical_inputs.npz':b.getvalue(),'evidence/frozen_fit_input_extracts.json':b'[]'}
        manifest={k:dict(bytes=len(v),sha256=fit.digest(v)) for k,v in members.items()}
        if bad_member:manifest['verification.json']['sha256']='wrong'
        with zipfile.ZipFile(f,'w') as z:
            for n,b in members.items():z.writestr(n,b)
            z.writestr('members-sha256.json',json.dumps(manifest))
        return f,p

    def test_projection_hash_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'projection_hash'):
                fit.load_inputs(*self.archive_pair(d,wrong_projection=True))

    def test_fixed_member_hash_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'fixed_member_hash'):
                fit.load_inputs(*self.archive_pair(d,bad_member=True))

    def test_failed_projection_status_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'projection_status'):
                fit.load_inputs(*self.archive_pair(d,bad_status=True))

    def test_reference_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'reference_binding'):
                fit.load_inputs(*self.archive_pair(d,bad_reference=True))

    def test_missing_roles_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'five_setting_roles'):
                fit.load_inputs(*self.archive_pair(d))


if __name__=='__main__':
    unittest.main()
