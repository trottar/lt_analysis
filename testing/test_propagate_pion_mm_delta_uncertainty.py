"""Deterministic uncertainty aggregation, range and source-identity guards."""
import argparse
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import zipfile
import numpy as np
if __package__:
    from . import propagate_pion_mm_delta_uncertainty as p
else:
    import propagate_pion_mm_delta_uncertainty as p


class SensitivityRanges(unittest.TestCase):
    def test_central_included_and_signed_variations_preserved(self):
        r=p.supported_envelope([1.,-2.],[[2.,-4.],[3.,-3.]],[True,True])
        self.assertEqual(r['minimum'],[1,-4]);self.assertEqual(r['maximum'],[3,-2])
        self.assertEqual(r['status'],'sensitivity_range_not_one_sigma')

    def test_unsupported_case_does_not_widen_range(self):
        r=p.supported_envelope(2.,np.array([1.,1000.]),[True,False])
        self.assertEqual((r['minimum'],r['maximum']),(1.,2.));self.assertEqual(r['omitted_cases'],1)

    def test_missing_family_is_unavailable_not_zero(self):
        r=p.supported_envelope([2.],np.empty((0,1)),[])
        self.assertIsNone(r['minimum']);self.assertIsNone(r['maximum'])

    def test_nonfinite_supported_case_rejected(self):
        with self.assertRaisesRegex(ValueError,'nonfinite_envelope'):p.supported_envelope(2.,[np.nan],[True])

    def test_nonfinite_unsupported_case_can_be_reported(self):
        r=p.supported_envelope(2.,[np.nan,3.],[False,True]);self.assertEqual(r['maximum'],3.)

    def test_axis_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError,'envelope_shape'):p.supported_envelope([1,2],[[3]],[True])


class CovarianceAggregation(unittest.TestCase):
    def test_cross_terms_retained_in_sum(self):
        F=np.array([[1.,2.],[-1.,-2.]])/np.sqrt(2)
        C=p.aggregate_factor(F,np.ones((2,1)))
        self.assertAlmostEqual(C[0,0],9.)  # Independent-bin substitution would give 5.

    def test_cancelling_correlated_bins(self):
        F=np.array([[1.,-1.],[-1.,1.]])/np.sqrt(2)
        self.assertAlmostEqual(p.aggregate_factor(F,np.ones((2,1)))[0,0],0.)

    def test_factor_matches_projected_replicas(self):
        A=np.random.default_rng(3).normal(size=(20,7));F=(A-A.mean(0))/np.sqrt(19)
        W=np.random.default_rng(4).normal(size=(7,3))
        np.testing.assert_allclose(p.aggregate_factor(F,W),np.cov((A@W).T))

    def test_bad_factor_axes_fail(self):
        with self.assertRaisesRegex(ValueError,'factor_axis'):p.aggregate_factor(np.ones((3,4)),np.ones((2,1)))

    def test_nonfinite_factor_fail(self):
        with self.assertRaisesRegex(ValueError,'nonfinite_factor'):p.aggregate_factor([[np.nan]],[[1.]])

    def test_undefined_correlation_is_null(self):
        self.assertIsNone(p.correlation([[0.,0.],[0.,1.]])[0][1])

    def test_single_response_sign_is_geometry_not_sampling(self):
        c=p.response_cosine([[1.,-2.,0.]])
        self.assertEqual(c[0][1],-1.);self.assertIsNone(c[0][2])


class ScientificBoundaries(unittest.TestCase):
    def test_unapproved_model_family_rejected(self):
        with self.assertRaisesRegex(ValueError,'unapproved_model_family'):p.varied_model({'Xfull':np.zeros((11,4,3))},{'kind':'new_acceptance_scale','beta':[1,2,3]})

    def test_delta_influence_does_not_change_shapes(self):
        X=np.ones((11,4,3));Y,b=p.varied_model({'Xfull':X},{'kind':'omit_delta_1','beta':[0.,2,3]})
        np.testing.assert_array_equal(X,Y);np.testing.assert_array_equal(b,[0,2,3]);self.assertIsNot(X,Y)

    def test_invalid_delta_index_rejected(self):
        with self.assertRaisesRegex(ValueError,'delta_case'):p.varied_model({'Xfull':np.zeros((11,4,3))},{'kind':'omit_delta_11','beta':[1,2,3]})

    def test_negative_case_amplitude_rejected(self):
        with self.assertRaisesRegex(ValueError,'case_amplitudes'):p.varied_model({'Xfull':np.zeros((11,4,3))},{'kind':'omit_delta_1','beta':[-1,2,3]})

    def test_native_boundaries_not_fractionally_filled(self):
        edges=np.array([.7,1.09,1.11,1.15,1.17,1.3]);q={'edges':edges,'y':np.ones((11,7))}
        X=np.ones((11,7,3));r=p.regional(q,X,np.array([1.,0,0]))
        self.assertEqual(r['lambda_complete']['data_total'],10.)
        self.assertEqual(r['broad']['data_total'],50.)
        self.assertEqual(r['high_MM_complete']['data_total'],10.)

    def test_refuses_existing_output_before_input_loading(self):
        with tempfile.TemporaryDirectory() as d:
            with mock.patch.object(p,'verified_bundle') as loader:
                with self.assertRaisesRegex(ValueError,'output_exists'):p.execute(argparse.Namespace(output=d))
                loader.assert_not_called()


class EvidenceIntegrity(unittest.TestCase):
    def test_member_mismatch_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'bad.zip'
            with zipfile.ZipFile(path,'w') as z:
                z.writestr('members-sha256.json',json.dumps({'result.json':{'bytes':2,'sha256':'wrong'}}));z.writestr('result.json','{}')
            with self.assertRaisesRegex(ValueError,'member_hash'):p.verified_bundle(path)

    def test_unbound_member_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'bad.zip'
            with zipfile.ZipFile(path,'w') as z:z.writestr('members-sha256.json','{}');z.writestr('extra','1')
            with self.assertRaisesRegex(ValueError,'unbound_member'):p.verified_bundle(path)

    def test_wrong_bootstrap_scope_rejected_before_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            with mock.patch.object(p,'verified_bundle',return_value=({'replicas':10,'full_scope':False},None,None,None)):
                with self.assertRaisesRegex(ValueError,'bootstrap_scope'):p.execute(argparse.Namespace(output=str(Path(d)/'new'),bootstrap_review='ignored'))

    def test_bootstrap_source_mismatch_rejected(self):
        old={'replicas':1024,'full_scope':True,'production_applied':False,'source_provenance':{'bootstrap_implementation_sha256':'different'}}
        with tempfile.TemporaryDirectory() as d:
            with mock.patch.object(p,'verified_bundle',return_value=(old,None,None,None)):
                with self.assertRaisesRegex(ValueError,'bootstrap_source_identity'):p.execute(argparse.Namespace(output=str(Path(d)/'new'),bootstrap_review='ignored'))


if __name__=='__main__':unittest.main()
