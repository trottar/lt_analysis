"""D1.a deterministic saved-source tests; no ROOT import or farm execution."""
import ast
import copy
import json
import math
from pathlib import Path
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FITS = ROOT / 'src/cuts/pion_component_fits.py'
PROTON = ROOT / 'src/cuts/proton_contamination_weights.py'


class Hist:
    def __init__(self, edges, values, errors):
        self.edges = list(edges)
        self.values = np.array(values, dtype=float)
        self.errors = np.array(errors, dtype=float)

    def GetNbinsX(self): return len(self.edges) - 1
    def GetXaxis(self): return self
    def GetXmin(self): return self.edges[0]
    def GetXmax(self): return self.edges[-1]
    def GetBinCenter(self, i): return (self.edges[i - 1] + self.edges[i]) / 2
    def GetBinLowEdge(self, i): return self.edges[i - 1]
    def GetBinContent(self, i): return float(self.values[i])
    def GetBinError(self, i): return float(self.errors[i])
    def GetName(self): return 'detached_test_hist'
    def Integral(self): return float(self.values[1:-1].sum())
    def FindBin(self, x): return int(np.searchsorted(self.edges, x, side='right'))
    def Add(self, other, scale=1):
        self.values += scale * other.values
        self.errors = np.sqrt(self.errors**2 + scale**2 * other.errors**2)
    def Scale(self, scale):
        self.values *= scale
        self.errors *= abs(scale)


def clone(h, name, reset=False):
    r = copy.deepcopy(h)
    if reset:
        r.values[:] = 0
        r.errors[:] = 0
    return r


def load_functions(path, names):
    nodes = [n for n in ast.parse(path.read_text(encoding='utf-8')).body
             if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in nodes} == set(names)
    ns = {'np': np, 'math': math, 'deepcopy': copy.deepcopy,
          '_clone_hist': clone, '_is_finite_number': lambda x: isinstance(x, (int, float)) and math.isfinite(x),
          '_component_plot_label': str}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), ns)
    return ns


FIT_NAMES = ['_validate_template_hist', '_build_single_template_fit_inputs',
             '_solve_nonnegative_template_amplitude', '_kaon_sidis_unsupported_bins',
             '_run_staged_component_pass']


class KaonSIDISOriginalSupportTests(unittest.TestCase):
    def setUp(self):
        self.ns = load_functions(FITS, FIT_NAMES)
        self.f = json.loads((ROOT / 'testing/fixtures/d1a_kaon_sidis_original_support.json').read_text())
        self.k, self.n, self.s = [Hist(self.f[x]['edges'], self.f[x]['values'], self.f[x]['errors'])
                                for x in ['kaon', 'neutron', 'sidis']]
        self.support = copy.deepcopy(self.f['support'])

    def run_stage(self, enabled=True, prefix='A', names=('pi_n', 'pi_sidis'), modes=None):
        windows = {x['component_name']: x['anchor_windows'] for x in self.f['stages']}
        windows['pi_delta'] = windows['pi_sidis']
        templates = {'pi_n': self.n, 'pi_sidis': self.s, 'pi_delta': self.s}
        return self.ns['_run_staged_component_pass'](
            self.k, templates, list(names), .7, 1.3, windows, {}, modes or {}, [], {},
            prefix, 0, kaon_sidis_experimental_support=self.support if enabled else None)

    def test_A_unchanged_eight_bin_boundary(self):
        amps, errors, _, steps = self.run_stage(False)
        self.assertAlmostEqual(amps['pi_n'], .182202244063977)
        self.assertEqual(amps['pi_sidis'], 0)
        d = steps[1]['fit_diagnostics']
        self.assertAlmostEqual(d['weighted_numerator'], -38.7723381448, places=8)
        self.assertEqual(d['fit_bin_indices'], list(range(61, 69)))

    def test_B_seven_bins_exact_weights_and_no_histogram_mutation(self):
        before = copy.deepcopy((self.k.values, self.k.errors, self.n.values, self.n.errors))
        amps, errors, _, steps = self.run_stage()
        d = steps[1]['fit_diagnostics']
        self.assertEqual(d['d1a_original_support']['excluded_bin_indices'], [62])
        self.assertEqual(d['fit_bin_indices'], [61, 63, 64, 65, 66, 67, 68])
        expected = self.run_stage(False)[3][1]['fit_diagnostics']['fit_sigma']
        self.assertTrue(np.array_equal(d['fit_sigma'], np.delete(expected, 1)))
        self.assertAlmostEqual(d['weighted_numerator'], 10.6094147226, places=8)
        self.assertAlmostEqual(d['weighted_denominator'], 48.4742198501456, places=10)
        self.assertAlmostEqual(amps['pi_sidis'], .218867157746124, places=12)
        for a, b in zip(before, (self.k.values, self.k.errors, self.n.values, self.n.errors)):
            self.assertTrue(np.array_equal(a, b))

    def test_C_signed_cancellation_is_supported(self):
        self.support['sources']['prompt']['entries'][62] = 1
        self.support['sources']['prompt']['sumw2'][62] = 1
        self.support['sources']['rand']['entries'][62] = 1
        self.support['sources']['rand']['sumw2'][62] = 1
        self.k.errors[62] = math.sqrt(2)
        _, _, _, steps = self.run_stage()
        self.assertIn(62, steps[1]['fit_diagnostics']['fit_bin_indices'])

    def test_D_supported_zero_and_negative_residuals_remain(self):
        _, _, _, steps = self.run_stage()
        d = steps[1]['fit_diagnostics']
        self.assertIn(61, d['fit_bin_indices'])
        self.assertLess(d['fit_y'][0], 0)
        self.k.values[61] = 0
        self.support['sources']['prompt']['entries'][61] = 1
        residual = copy.deepcopy(self.k)
        residual.errors[61] = 1e-10
        self.assertNotIn(61, self.ns['_kaon_sidis_unsupported_bins'](self.k, residual, self.support))

    def test_E_missing_and_zero_template_fail(self):
        for t in (None, Hist(self.k.edges, np.zeros(102), np.zeros(102))):
            r = self.ns['_solve_nonnegative_template_amplitude'](self.k, t, .7, 1.3)
            self.assertFalse(r['success'])
        self.s.values[:] = 0
        with self.assertRaisesRegex(RuntimeError, 'D1a_SIDIS_fit_failed'):
            self.run_stage()

    def test_E_template_without_anchor_support_fails(self):
        self.s.values[:] = 0
        self.s.values[2] = 1
        with self.assertRaisesRegex(RuntimeError, 'zero support inside anchor'):
            self.run_stage()

    def test_F_insufficient_supported_bins_fails(self):
        for i in range(61, 69):
            if i == 68: continue
            self.k.values[i] = self.k.errors[i] = 0
            for source in self.support['sources'].values():
                source['entries'][i] = 0
                source['sumw2'][i] = 0
        with self.assertRaisesRegex(RuntimeError, 'insufficient supported fit bins'):
            self.run_stage()

    def test_G_pion_neutron_delta_window_integral_unchanged(self):
        for prefix, names, modes in [('B', ('pi_n', 'pi_sidis'), {}),
                                      ('A', ('pi_n',), {}),
                                      ('A', ('pi_delta',), {}),
                                      ('A', ('pi_sidis',), {'pi_sidis': 'window_integral'})]:
            old = self.run_stage(False, prefix, names, modes)
            new = self.run_stage(True, prefix, names, modes)
            self.assertEqual(old[0], new[0])
            self.assertEqual(old[1], new[1])
            self.assertEqual(old[2], new[2])
            self.assertEqual([x['fit_diagnostics'] for x in old[3]], [x['fit_diagnostics'] for x in new[3]])

    def test_missing_malformed_or_axis_mismatched_support_fails(self):
        for support in ({}, {**self.support, 'complete': False}, {**self.support, 'edges': [0, 1]}):
            self.support = support
            with self.assertRaisesRegex(RuntimeError, 'D1a_original_experimental_support'):
                self.run_stage()

    def test_contradictory_empty_support_fails(self):
        self.k.errors[62] = 1
        with self.assertRaisesRegex(RuntimeError, 'conflicts_with_target'):
            self.run_stage()

    def test_producer_preserves_zero_cleaned_weight_observation_and_flow(self):
        p = load_functions(PROTON, ['_new_kaon_sidis_experimental_support', '_record_kaon_sidis_experimental_support'])
        specs = {'sources': {n: {'tree': [], 'coefficient': 1} for n in ['prompt', 'rand', 'dummy_prompt', 'dummy_rand']}}
        support = p['_new_kaon_sidis_experimental_support'](self.k, specs)
        # Metadata uses raw coefficient; even a zero cleaned factor cannot erase an entry.
        p['_record_kaon_sidis_experimental_support'](support, self.k, 'prompt', 1.069, 2)
        p['_record_kaon_sidis_experimental_support'](support, self.k, 'dummy_prompt', 1.5, -3)
        self.assertEqual(support['sources']['prompt']['entries'][62], 1)
        self.assertEqual(support['sources']['prompt']['sumw2'][62], 4)
        self.assertEqual(support['sources']['dummy']['entries'][-1], 1)
        specs['sources']['dummy_prompt']['tree'] = None
        self.assertFalse(p['_new_kaon_sidis_experimental_support'](self.k, specs)['complete'])

    def test_guarded_runtime_plumbing(self):
        tree = ast.parse(FITS.read_text(encoding='utf-8'))
        entry = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'fit_kaon_nosub_with_simc_pion_shapes')
        text = ast.get_source_segment(FITS.read_text(encoding='utf-8'), entry)
        self.assertIn('D1a_requires_staged_only_kaon_fit', text)
        self.assertIn('D1a_requires_fixed_alignment_no_residual_shift_selection', text)
        self.assertIn('and not _alignment_scoring_only', text)
        self.assertIn('kaon_experimental_support if d1a_enabled else None', text)
        producer = ast.parse(PROTON.read_text(encoding='utf-8'))
        fn = next(n for n in producer.body if isinstance(n, ast.FunctionDef) and n.name == 'apply_kaon_proton_cleaning_to_targets')
        rf = [n for n in ast.walk(fn) if isinstance(n, ast.If) and ast.unparse(n.test) == 'rf_accept']
        self.assertEqual(len(rf), 1)
        self.assertIn('collect_kaon_sidis_support and nommcuts', ast.unparse(rf[0]))
        self.assertIn('source_name, adj_mm, raw_weight', ast.unparse(rf[0]))

    def test_actual_kaon_entry_fails_closed_on_missing_support_or_unsupported_mode(self):
        entry = load_functions(FITS, ['fit_kaon_nosub_with_simc_pion_shapes'])
        entry.update({'BG_OPT_MM_PLOT_MIN': .7, 'BG_OPT_MM_PLOT_MAX': 1.3,
                      'KAON_SIGMA0_TEMPLATE_NAME': 'k_sigma0', 'KAON_SIGNAL_TEMPLATE_NAME': 'k_lambda',
                      'COMPONENT_NAMES': ('pi_n', 'pi_delta', 'pi_sidis')})
        for name in ['get_particle_subtraction_component_fit_window_config',
                     'resolve_particle_subtraction_component_fit_windows',
                     'resolve_particle_subtraction_component_stage_amplitude_windows',
                     'resolve_particle_subtraction_component_stage_amplitude_modes',
                     'resolve_particle_subtraction_component_prior_scales',
                     '_resolve_pi_delta_signal_protected_fit_config']:
            entry[name] = lambda *a, **k: {}
        entry['resolve_particle_subtraction_component_fit_excluded_windows'] = lambda *a, **k: []
        entry['resolve_particle_subtraction_component_fit_mode'] = lambda *a, **k: 'staged_only'
        entry['resolve_particle_subtraction_component_residual_shift_settings'] = lambda *a, **k: {'enabled': True}
        fn = entry['fit_kaon_nosub_with_simc_pion_shapes']
        inp = {'d1a_kaon_sidis_support_enabled': True}
        with self.assertRaisesRegex(RuntimeError, 'support_unavailable_or_invalid'):
            fn(self.k, self.n, self.s, self.s, None, None, inp)
        with self.assertRaisesRegex(RuntimeError, 'requires_fixed_alignment'):
            fn(self.k, self.n, self.s, self.s, None, None, inp, kaon_experimental_support=self.support)
        entry['resolve_particle_subtraction_component_fit_mode'] = lambda *a, **k: 'joint_refined'
        with self.assertRaisesRegex(RuntimeError, 'requires_staged_only'):
            fn(self.k, self.n, self.s, self.s, None, None, inp, kaon_experimental_support=self.support)


if __name__ == '__main__':
    unittest.main()
