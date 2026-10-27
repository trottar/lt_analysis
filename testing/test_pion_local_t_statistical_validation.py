"""Deterministic conditional-statistics checks; no physical truth/farm acceptance."""
from copy import deepcopy
import json
import math
import types
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

from testing import pion_local_t_statistical_validation as stat

DOMAINS={'all':None,'signal':(1.1,1.2)}
EDGES=[.7,1.1,1.2,1.3]


class CovarianceTests(unittest.TestCase):
    def weights(self):return {'A':[0.,2.,2.,2.,0.],'B':[0.,3.,3.,3.,0.],'C':[0.,4.,4.,4.,0.]}

    def test_signed_cancellation_keeps_quadratic_source_covariance(self):
        # The four source signs affect the mean; their squares affect covariance.
        coeff=[2.,-1.,-.5,.25]
        s2=math.fsum(w*w for w in coeff)
        result=stat.conditional_covariance([0.,s2,0.,0.,0.],EDGES,self.weights(),DOMAINS)
        d=result['domains']['all']
        self.assertEqual(d['covariance'][0][1],s2*6)
        self.assertEqual(d['paired_difference_variance_from_B'][0],s2)
        self.assertEqual(d['naive_independent_variance_from_B'][0],s2*13)
        self.assertEqual(d['paired_difference_variance_from_B'][1],0.)
        cancelled=stat.conditional_covariance([0.,2.,0.,0.,0.],EDGES,self.weights(),DOMAINS)
        self.assertEqual(cancelled['domains']['all']['covariance'][0][0],8.)
        self.assertFalse(result['independent_model_validation'])

    def test_cut_domains_flows_and_linear_source_phi_closure(self):
        weights=self.weights()
        a=[10.,2.,3.,4.,20.]; b=[5.,5.,7.,1.,30.]
        left=stat.conditional_covariance(a,EDGES,weights,DOMAINS)
        right=stat.conditional_covariance(b,EDGES,weights,DOMAINS)
        joined=stat.conditional_covariance([x+y for x,y in zip(a,b)],EDGES,weights,DOMAINS)
        for domain in DOMAINS:
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual(joined['domains'][domain]['covariance'][i][j],
                        left['domains'][domain]['covariance'][i][j]+right['domains'][domain]['covariance'][i][j])
        self.assertEqual(left['domains']['signal']['covariance'][1][1],27.)
        self.assertEqual(left['domains']['all']['covariance'][1][1],81.)

    def test_PSD_and_same_event_identity(self):
        w={'A':[0.,1.,-2.,3.,0.],'B':[0.,2.,3.,1.,0.],'C':[0.,4.,1.,2.,0.]}
        s2=[0.,2.,4.,3.,0.]
        d=stat.conditional_covariance(s2,EDGES,w,DOMAINS)['domains']['all']
        names=list(w)
        for vector in ([1.,2.,3.],[-1.,4.,-2.],[0.,1.,-1.]):
            quad=math.fsum(vector[i]*d['covariance'][i][j]*vector[j] for i in range(3) for j in range(3))
            direct=math.fsum(s2[k]*math.fsum(vector[i]*w[n][k] for i,n in enumerate(names))**2 for k in range(1,4))
            self.assertGreaterEqual(quad,0.)
            self.assertAlmostEqual(quad,direct)

    def test_nearly_identical_models_do_not_lose_paired_variance(self):
        weights={'B':[0.,1.,1.,1.,0.],'C':[0.,1.+1e-12,1.,1.,0.]}
        d=stat.conditional_covariance([0.,2.,0.,0.,0.],EDGES,weights,DOMAINS)['domains']['all']
        self.assertGreater(d['paired_difference_variance_from_B'][1],0.)
        self.assertEqual(d['paired_difference_variance_from_B'][1],2*(weights['C'][1]-1)**2)

    def test_correlation_avoids_variance_product_overflow_and_underflow(self):
        for scale in (1e-300,1e300):
            weights={'B':[0.,1.,1.,1.,0.],'C':[0.,1.,1.,1.,0.]}
            d=stat.conditional_covariance([0.,scale,0.,0.,0.],EDGES,weights,DOMAINS)['domains']['all']
            self.assertAlmostEqual(d['correlation'][0][1],1.)

    def test_zero_support_is_explicit_and_inputs_immutable(self):
        weights=self.weights();before=deepcopy(weights)
        result=stat.conditional_covariance([0.]*5,EDGES,weights,DOMAINS)
        self.assertEqual(weights,before)
        self.assertEqual(result['domains']['all']['correlation'],[[None]*3]*3)
        self.assertEqual(result['domains']['all']['paired_difference_variance_from_B'],[0.]*3)

    def test_bad_geometry_nonfinite_variance_weights_and_flows_fail_closed(self):
        for variance,edges,weights in (
            ([0.,-1.,0.,0.,0.],EDGES,self.weights()),
            ([0.,float('nan'),0.,0.,0.],EDGES,self.weights()),
            ([0.]*5,[.7,.7,1.2,1.3],self.weights()),
            ([0.]*4,EDGES,self.weights()),
            ([0.]*5,EDGES,{'B':[0.,1.,None,1.,0.]}),
            ([0.]*5,EDGES,{'B':[1.,1.,1.,1.,0.]}),
            ([0.]*5,EDGES,{'A':[0.,1.,1.,1.,0.]})):
            with self.subTest(variance=variance,edges=edges,weights=weights),self.assertRaises(ValueError):
                stat.conditional_covariance(variance,edges,weights,DOMAINS)


class DependencyTests(unittest.TestCase):
    def test_full_sample_models_and_disjoint_folds_do_not_clear_dependency_gate(self):
        parent={'candidate_B':{'pion_component_alignment':{'training_only':True},'diagnostics':{'folds_disjoint':True},
                              'H_kaon_pion_bg_fit_total_refined':{'values':[1.]}}}
        before=deepcopy(parent);result=stat.training_dependency_audit(parent)
        self.assertEqual(parent,before)
        self.assertEqual(result['status'],'unavailable')
        self.assertFalse(result['independent_predictive_test_implemented'])
        self.assertFalse(result['disjoint_record_folds_are_sufficient'])
        self.assertFalse(result['established_fitter']['complete_training_chain_from_current_packet'])
        self.assertIn('raw learning populations',result['first_missing_step'])
        self.assertEqual(result['physical_bias'],'bias_truth_unavailable')
        self.assertEqual(result['MSE'],'MSE_selection_unavailable')
        self.assertIsNone(result['selected_lambda'])


class ExistingAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        global c
        from testing import pion_local_t_transfer_variance_comparison as c
        from testing.test_pion_local_t_transfer_variance_comparison import synthetic_documents
        cls.docs=synthetic_documents()
        cls.provenance={s:{'sha256':c.sha(c.encoded(cls.docs[s]))} for s in c.SETTINGS}
        cls.candidate=json.loads(c.encoded(c.compare(cls.docs,cls.provenance)))

    def test_existing_arrays_audit_never_refits_or_selects(self):
        before=c.encoded([self.docs,self.candidate])
        with patch.object(c,'least_squares',side_effect=AssertionError('no new fit')):
            result=stat.audit_existing(self.docs,self.provenance,self.candidate,c.SOURCE)
        self.assertEqual(c.encoded([self.docs,self.candidate]),before)
        self.assertEqual(len(result['parents']),6)
        self.assertIsNone(result['selected_lambda'])
        self.assertIn('unavailable',result['event_allcuts_validation'])

    def test_hash_saved_parent_selected_strength_policy_and_ratio_mismatch_rejected(self):
        for mutate in (
            lambda d:d.update(source_HEAD_of_input='0'*40),
            lambda d:d['provenance']['Center'].update(sha256='0'*64),
            lambda d:d['parents'][0]['A_B_unchanged'].update(t_bin_index=9),
            lambda d:d['parents'][0]['C'][0].update(selected=True),
            lambda d:d['policy'].update(selection='pick best'),
            lambda d:d.update(implementation_sha256='0'*64),
            lambda d:d['parents'][0]['C'][0]['transfer']['weights'].__setitem__(1,100.)):
            bad=deepcopy(self.candidate);mutate(bad)
            with self.subTest(mutate=mutate),self.assertRaises(ValueError):
                stat.audit_existing(self.docs,self.provenance,bad,c.SOURCE)

    def test_CLI_rejects_unreviewed_candidate_before_loading_fitting_or_output(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)
            inputs=[p/c.basename(s) for s in c.SETTINGS]+[p/'candidate.json']
            for item in inputs:item.write_text('{}')
            args=['--center-json',str(inputs[0]),'--left-json',str(inputs[1]),
                  '--candidate-json',str(inputs[2]),'--candidate-sha256','0'*64,
                  '--source-commit',c.SOURCE,'--output-json',str(p/'result.json'),
                  '--output-pdf',str(p/'result.pdf')]
            with patch.object(c,'load_inputs',side_effect=AssertionError('do not load unreviewed candidate')) as load:
                with self.assertRaisesRegex(ValueError,'reviewed_candidate_hash_mismatch'):stat.main(args)
                load.assert_not_called()
            self.assertFalse((p/'result.json').exists())
            self.assertFalse((p/'result.pdf').exists())


if __name__=='__main__':unittest.main()
