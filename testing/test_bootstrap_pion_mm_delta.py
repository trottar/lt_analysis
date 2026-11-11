"""Deterministic shared-event, normalization and covariance bootstrap tests."""
import unittest
from unittest import mock
import numpy as np
if __package__:
    from . import bootstrap_pion_mm_delta as boot
else:
    import bootstrap_pion_mm_delta as boot


class EventResampling(unittest.TestCase):
    def test_entry_alone_is_not_identity(self):
        n,ix=boot.identity_groups([0,1,0],['treeA','treeA','treeB'],[4,4,4])
        self.assertEqual(n,3)
        self.assertEqual(len(set(ix)),3)

    def test_shared_record_uses_one_draw_across_views(self):
        g=boot.identity_groups([0,0,1],['tree','tree','tree'],[7,7,9])
        rng=mock.Mock();rng.poisson.return_value=np.array([2,5])
        np.testing.assert_array_equal(boot.multiplicities(rng,g),[2,2,5])
        rng.poisson.assert_called_once_with(1.,2)

    def test_signed_moments_use_repeated_event_Sumw2(self):
        y,v=boot.moments([0,0],[2.,-3.],[4,2],1)
        self.assertEqual(y[0],2)
        self.assertEqual(v[0],34)  # 4*2^2 + 2*(-3)^2, never k^2.

    def test_all_ones_replays_original(self):
        y,v=boot.moments([0,1,0],[1.,-2.,3.],[1,1,1],2)
        np.testing.assert_array_equal(y,[4,-2]);np.testing.assert_array_equal(v,[10,4])

    def test_zero_draw_stays_zero_without_variance_floor(self):
        y,v=boot.moments([0],[7.],[0],1)
        self.assertEqual(y[0],0);self.assertEqual(v[0],0)

    def test_fractional_draw_rejected(self):
        with self.assertRaisesRegex(ValueError,'multiplicity'):
            boot.moments([0],[2.],[.5],1)

    def test_reproducible_draws(self):
        g=boot.identity_groups([0,0,0],['t']*3,[1,2,3])
        np.testing.assert_array_equal(boot.multiplicities(np.random.default_rng(7),g),boot.multiplicities(np.random.default_rng(7),g))

    def test_MC_and_data_generators_are_separate(self):
        a=np.random.default_rng(8).poisson(1,200);b=np.random.default_rng(9).poisson(1,200)
        self.assertFalse(np.array_equal(a,b))


class TemplateNormalization(unittest.TestCase):
    def test_shared_parent_denominator_rebuilt_once(self):
        S=np.array([[0.,4,6,0],[0.,8,2,0]]);V=S*3;L=np.eye(4)/10
        terms=boot.rebuilt_terms(S,V,L)
        Q=terms[2]
        self.assertAlmostEqual(Q[:,1:-1].sum(),1)
        self.assertAlmostEqual(Q[0,1:-1].sum(),.5)
        self.assertAlmostEqual(Q[1,1:-1].sum(),.5)

    def test_changed_parent_population_changes_slice_fractions(self):
        S=np.array([[0.,2,2,0],[0.,6,10,0]])
        Q=boot.rebuilt_terms(S,S,np.eye(4)/50)[2]
        self.assertAlmostEqual(Q[0].sum(),.2)
        self.assertAlmostEqual(Q[1].sum(),.8)

    def test_zero_MC_parent_is_failure(self):
        with self.assertRaisesRegex(ValueError,'empty_MC_parent'):
            boot.rebuilt_terms(np.zeros((2,4)),np.zeros((2,4)),np.eye(4))

    def test_correlated_normalization_covariance_PSD(self):
        S=np.array([[0.,4,6,0],[0.,8,2,0]])
        t=boot.rebuilt_terms(S,S,np.eye(4)/10)
        C=boot.fit.mc_covariance([t],[2.],np.arange(8))
        self.assertGreater(np.linalg.eigvalsh(C).min(),-1e-13)
        self.assertLess(C[1:3,5:7].sum(),0)


class ExistingTransfer(unittest.TestCase):
    def test_ratio_and_zero_flows(self):
        w,u=boot.transfer_ratio([9,2,0,8],[3,4,0,2],1e-12,0,None)
        np.testing.assert_array_equal(w,[0,.5,0,0]);self.assertFalse(u.any())

    def test_shared_numerator_denominator_cancels_MC(self):
        for D in ([0.,.2,.8,0.],[0.,.9,.1,0.]):
            w,_=boot.transfer_ratio(D,D,1e-12,0,None)
            np.testing.assert_array_equal(w,[0,1,1,0])

    def test_unsupported_denominator_reported_not_floored(self):
        w,u=boot.transfer_ratio([0.,2,0,0],[0.,0,0,0],1e-12,0,None)
        self.assertEqual(w[1],0);self.assertTrue(u[1])

    def test_recorded_clipping_used(self):
        w,_=boot.transfer_ratio([0.,4,2,0],[0.,1,4,0],1e-12,.1,2.)
        np.testing.assert_array_equal(w,[0,2,.5,0])

    def test_snapshot_None_is_explicit_source_semantics(self):
        p=boot.policy_from_record({'clip_max':{'status':'unavailable','reason':'not retained by producer'},'clip_min':0.,'denom_floor':1e-12})
        self.assertIsNone(p['high'])

    def test_missing_upper_policy_rejected(self):
        with self.assertRaisesRegex(ValueError,'missing_transfer_policy'):
            boot.policy_from_record({'clip_min':0.,'denom_floor':1e-12})

    def test_unknown_upper_marker_rejected(self):
        with self.assertRaisesRegex(ValueError,'unknown_clip_policy'):
            boot.policy_from_record({'clip_min':0.,'denom_floor':1e-12,'clip_max':{'unknown':True}})

    def test_signed_background_reuses_multiplier(self):
        rows=dict(mm=np.array([.8,.8,1.1]),pi=np.array([0,0,1]),ti=np.array([0,0,0]),nommcuts=np.array([1,1,1],bool),allcuts=np.array([0,0,1],bool),coefficient=np.array([2.,-1.,3.]))
        grid,phi,valid,_=boot.background(rows,np.array([2,1,3]),1,np.array([.7,1.,1.3]),np.array([0.,2,4,0]),np.array([0.,4,2,0]),dict(floor=1e-12,low=0,high=None))
        self.assertAlmostEqual(grid[0,1],1.5);self.assertAlmostEqual(phi[1],18.)

    def test_broad_failure_does_not_discard_supported_Lambda(self):
        rows=dict(mm=np.array([.8,1.1]),pi=np.array([0,0]),ti=np.array([0,0]),nommcuts=np.array([1,1],bool),allcuts=np.array([0,1],bool),coefficient=np.ones(2))
        grid,phi,valid,_=boot.background(rows,np.ones(2),1,np.array([.7,1.,1.3]),np.array([0.,1,1,0]),np.array([0.,0,2,0]),dict(floor=1e-12,low=0,high=None))
        self.assertFalse(valid['broad']);self.assertTrue(valid['lambda_region']);self.assertEqual(phi[0],.5)


class CovarianceAlgebra(unittest.TestCase):
    def test_known_sample_covariance(self):
        C=boot.covariance(np.array([[1.,2.],[2.,4.],[3.,6.]]))
        np.testing.assert_array_equal(C,[[1,2],[2,4]])

    def test_covariance_requires_common_supported_replicas(self):
        a=np.array([[1.,2.],[2.,np.nan],[3.,6.]])
        r=boot.covariance_summary(a)
        self.assertEqual(r['common_replica_indices'],[0,2])
        np.testing.assert_array_equal(r['covariance'],[[2.,4.],[4.,8.]])

    def test_insufficient_support_not_zero_covariance(self):
        r=boot.covariance_summary([[1,2],[np.nan,np.nan]])
        self.assertIsNone(r['covariance'])

    def test_undefined_correlation_not_independence(self):
        r=boot.covariance_summary([[0.,1],[0.,2],[0.,3]])
        self.assertIsNone(r['correlation'][0][1])

    def test_nonlinear_paired_interaction_retained(self):
        D=np.array([[1.],[2.],[3.],[4.]]);M=np.array([[2.],[4.],[1.],[3.]])
        J=D+M+.1*D*M
        r=boot.paired_decomposition(D,M,J,[0.])
        self.assertEqual(r['status'],'available_conditional')
        self.assertGreater(r['covariance_nonlinear_interaction'][0][0],0)
        self.assertLess(r['closure_max_error'],1e-12)

    def test_factor_exactly_reconstructs_covariance(self):
        A=np.arange(30,dtype=float).reshape(10,3);F=(A-A.mean(0))/3
        np.testing.assert_allclose(F.T@F,boot.covariance(A))

    def test_covariance_rejects_nonfinite(self):
        with self.assertRaisesRegex(ValueError,'covariance_support'):
            boot.covariance([[1.,np.nan],[2.,3.]])


if __name__=='__main__':unittest.main()
