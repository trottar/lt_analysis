"""Synthetic private-path and failure gates; never launch the analysis or farm."""
from copy import deepcopy
from contextlib import contextmanager, ExitStack
import os
from pathlib import Path
import tempfile
import subprocess
import sys
import json
import ast
from array import array
import unittest
from unittest.mock import patch
from unittest.mock import Mock
from types import SimpleNamespace

from testing import run_pion_local_t_matched_cohort_gate as owner
from testing import pion_local_t_event_validation_capture as cap

ROOT = Path(__file__).resolve().parents[1]


def external_fixture(worktree, private, volatile):
    """Read-only inputs for the real inventory, without ROOT or farm access."""
    worktree.mkdir(parents=True, exist_ok=True)
    private.mkdir(parents=True, exist_ok=True)
    (worktree / 'run_Prod_Analysis.sh').write_text('EffData="synthetic.csv"\n')
    paths = {'VOLATILEPATH': str(volatile), 'UTILPATH': str(worktree.parent / 'util'),
             'PARAMPATH': str(worktree.parent / 'params'), 'CUTPATH': str(worktree.parent / 'cuts')}
    def write(path, contents=b'synthetic input'):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(contents)
    write(Path(paths['UTILPATH']) / 'scripts/efficiency/OUTPUTS/synthetic.csv')
    write(Path(paths['PARAMPATH']) / 'Timing_Parameters.csv')
    for setting in owner.SETTINGS:
        write(private / f'{setting}_kaon_Analysed_Data_Q4p4W2p74_lowe.root')
        write(private / f'Prod_Coin_Q4p4W2p74{setting.lower()}_lowe.root')
        for dummy in ('', '_dummy'):
            write(Path(paths['UTILPATH']) / 'run_list/KaonLT' /
                  f'Q4p4W2p74{setting.lower()}_lowe{dummy}', b'1234\n')
        for component in ('neutron', 'delta', 'sidis'):
            write(volatile / 'OUTPUT/Analysis/KaonLT/background_samples' / component /
                  f'Prod_Coin_Q4p4W2p74{setting.lower()}_lowe.root')
    return paths


class PrivateBoundary(unittest.TestCase):
    def test_sigma0_unconfigured_inventory_preserves_fallback(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d).resolve()
            worktree, private = root / 'analysis', root / 'private'
            paths = external_fixture(worktree, private, root / 'volatile')
            with patch.object(owner, 'skim_paths', return_value={'effective_charge': str(root / 'skim')}):
                for value in (None, '', '  \t'):
                    with self.subTest(value=value):
                        env = {} if value is None else {
                            f'LT_BG_SIGMA0_{s.upper()}_LOW_ROOT': value for s in owner.SETTINGS}
                        records = owner.external_input_inventory(worktree, paths, private, env)
                        self.assertNotIn('sigma0_SIMC', [r['role'] for r in records.values()])
                        self.assertEqual(sum(r['role'].startswith('pion_SIMC_') for r in records.values()), 6)
                        owner.verify_external_inputs(records)

    def test_sigma0_configured_inventory_checks_identity_and_preservation(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d).resolve()
            worktree, private = root / 'analysis', root / 'private'
            paths = external_fixture(worktree, private, root / 'volatile')
            env = {}
            for s in owner.SETTINGS:
                p = root / (s + '-sigma.root')
                p.write_bytes((s + ' Sigma0').encode())
                env[f'LT_BG_SIGMA0_{s.upper()}_LOW_ROOT'] = '  ' + str(p) + '  '
            with patch.object(owner, 'skim_paths', return_value={'effective_charge': str(root / 'skim')}):
                records = owner.external_input_inventory(worktree, paths, private, env)
            sigma = {p: r for p, r in records.items() if r['role'] == 'sigma0_SIMC'}
            self.assertEqual(len(sigma), 2)
            for p, record in sigma.items():
                self.assertEqual(record['identity'], owner.identity(p))
            owner.verify_external_inputs(records)
            Path(next(iter(sigma))).write_bytes(b'changed configured source')
            with self.assertRaisesRegex(ValueError, 'external_input_changed'):
                owner.verify_external_inputs(records)

    def test_sigma0_configured_invalid_inventory_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d).resolve()
            worktree, private = root / 'analysis', root / 'private'
            paths = external_fixture(worktree, private, root / 'volatile')
            for setting in owner.SETTINGS:
                for bad in (root / 'missing.root', root):
                    with self.subTest(setting=setting, path=bad), patch.object(
                            owner, 'skim_paths', return_value={'effective_charge': str(root / 'skim')}):
                        with self.assertRaisesRegex(ValueError, 'required_external_input_missing:sigma0_SIMC'):
                            owner.external_input_inventory(worktree, paths, private, {
                                f'LT_BG_SIGMA0_{setting.upper()}_LOW_ROOT': str(bad)})

    def test_frozen_private_derived_branch_algorithm_preserves_raw_rows(self):
        source = ROOT / 'src/setup/shift_MM.py'
        tree = ast.parse(source.read_text())
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'build_multi_derived_tree')
        namespace = {'array': array}
        exec(compile(ast.Module(body=[function], type_ignores=[]), str(source), 'exec'), namespace)
        raw = [{'MM': 1.1, 'MandelT': -.2, 'PID': 0, 'MM_shift': 7., 't_shift': 8.},
               {'MM': 1.2, 'MandelT': -.4, 'PID': 3, 'MM_shift': 9., 't_shift': 10.}]
        before = deepcopy(raw)
        class NewTree:
            def __init__(self, parent):
                self.parent = parent
                self.kept = [k for k in raw[0] if parent.status[k]]
                self.branches, self.rows = {}, []
            def Branch(self, name, value, form): self.branches[name] = value
            def Fill(self):
                self.rows.append({**{k: self.parent.current[k] for k in self.kept},
                                  **{k: value[0] for k, value in self.branches.items()}})
        class InputTree:
            def __init__(self): self.status = {k: 1 for k in raw[0]}
            def GetBranch(self, name): return name in raw[0]
            def SetBranchStatus(self, name, value): self.status[name] = value
            def CloneTree(self, count): return NewTree(self)
            def __iter__(self):
                for row in raw:
                    self.current = row
                    yield SimpleNamespace(**row)
        original = InputTree()
        result = namespace['build_multi_derived_tree'](original, [
            ('MM_shift', lambda evt: evt.MM + .01), ('t_shift', lambda evt: -evt.MandelT + .02)])
        self.assertEqual(raw, before)
        self.assertEqual([{k: row[k] for k in ('MM', 'MandelT', 'PID')} for row in result.rows],
                         [{k: row[k] for k in ('MM', 'MandelT', 'PID')} for row in raw])
        self.assertAlmostEqual(result.rows[0]['MM_shift'], 1.11, places=6)
        self.assertAlmostEqual(result.rows[1]['t_shift'], .42, places=6)
        self.assertTrue(all(original.status.values()))

    def test_failed_child_cleans_owned_context_and_never_collects(self):
        self.exercise_exact_preflight_sigma0('absent')

    def test_exact_preflight_configured_sigma0_reaches_child(self):
        self.exercise_exact_preflight_sigma0('valid')

    def test_exact_preflight_invalid_sigma0_stops_before_child(self):
        self.exercise_exact_preflight_sigma0('invalid')

    def exercise_exact_preflight_sigma0(self, sigma_case):
        with tempfile.TemporaryDirectory() as d:
            parent = Path(d).resolve()
            repo, worktree, globus = parent / 'repo', parent / 'analysis', parent / 'globus'
            for p in (repo, worktree, globus, parent / 'runtime-python'):
                p.mkdir()
            (parent / 'runtime-python/ltsep').mkdir()
            (worktree / 'testing').mkdir()
            (worktree / owner.OWNER).write_text('synthetic owner')
            for name in ('src/cuts/shift_prep.py', 'src/setup/shift_MM.py'):
                p = worktree / name
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text('# frozen writer fixture')
            shared = parent / 'volatile/OUTPUT/Analysis/KaonLT'
            (shared / 'background_samples').mkdir(parents=True)
            for setting in owner.SETTINGS:
                for role in ('Data', 'Dummy'):
                    (shared / f'{setting}_kaon_Analysed_{role}_Q4p4W2p74_lowe.root').write_bytes((setting + role).encode())
                for suffix in ('.root', '.hist'):
                    (shared / f'Prod_Coin_Q4p4W2p74{setting.lower()}_lowe{suffix}').write_bytes(b'simc')
            # Exercise execute_gate -> external_input_inventory directly;
            # only the eventual analysis child is simulated and fails.
            external_paths = external_fixture(worktree, shared, parent / 'volatile')
            sigma_env = {}
            if sigma_case != 'absent':
                for s in owner.SETTINGS:
                    p = parent / (s + '-sigma.root')
                    if sigma_case == 'valid':
                        p.write_bytes((s + ' Sigma0').encode())
                    sigma_env[f'LT_BG_SIGMA0_{s.upper()}_LOW_ROOT'] = str(p)
            original = owner.inventory(shared)
            package = parent / 'installed'
            package.mkdir()
            for name in ('__init__.py', 'pathing.py'):
                (package / name).write_text('# synthetic only')
            head = 'c'*40
            before = {'head': head, 'origin_test': head, 'branch': 'test', 'porcelain': ''}
            paths = {'VOLATILEPATH': str(parent / 'volatile'), 'SKIMPATH': str(shared),
                     'LTANAPATH': str(worktree), 'ANATYPE': 'Kaon',
                     'OUTPATH': str(worktree / 'OUTPUT/Analysis/KaonLT'), **external_paths}
            cleaned = []
            @contextmanager
            def owned(*args):
                try:
                    yield worktree
                finally:
                    cleaned.append(args[-1])
            source_checks = Mock(return_value=('', [{'returncode': 0}]))
            collect = Mock(side_effect=AssertionError('failed child must not collect'))
            @contextmanager
            def collection(*args):
                yield SimpleNamespace(collect_source_checks=source_checks, collect_validation_bundle=collect)
            with ExitStack() as stack:
                def patched(obj, name, **kwargs):
                    return stack.enter_context(patch.object(obj, name, **kwargs))
                patched(owner.isolation.accepted, 'REPO', new=repo)
                patched(owner.isolation.accepted, 'ARTIFACTS', new=shared)
                patched(owner.isolation.accepted, 'GLOBUS', new=globus)
                patched(owner.isolation, 'snapshot', return_value=before)
                patched(owner.ab, 'git', return_value='')
                patched(owner.ab, 'GateStatus')
                patched(owner, 'toolchain', return_value={'synthetic': True})
                patched(owner, 'audit_default_output_writer', return_value={})
                patched(owner.isolation, 'owned_worktree', new=owned)
                patched(owner.isolation, 'collection_module', new=collection)
                patched(owner.isolation, 'import_identity', return_value={
                    'package_file': str(package / '__init__.py'), 'pathing_file': str(package / 'pathing.py')})
                patched(owner.isolation, 'prepare_runtime_overlay', return_value=(sigma_env, {
                    'baseline_paths': dict(paths, LTANAPATH=str(repo), OUTPATH=str(shared))}))
                patched(owner.isolation, 'probe_paths', return_value=[{'paths': paths}])
                patched(owner.isolation, 'external_symlink_preflight', return_value={'links': [], 'background_simc_path': ''})
                inventory_call = patched(owner, 'external_input_inventory', wraps=owner.external_input_inventory)
                patched(owner, 'skim_paths', return_value={'launcher': str(shared), 'effective_charge': str(shared)})
                patched(owner, 'probe_private_callers', return_value={})
                patched(owner.event_owner, 'committed_identities', return_value={})
                def failed_child(tree, log, env):
                    Path(log).write_text('synthetic child failure; no analysis executed')
                    return 1
                child = patched(owner.ab.debug, 'run_analysis', side_effect=failed_child)
                checkout = patched(owner.isolation, 'verify_preservation', return_value={})
                installed = patched(owner.isolation, 'verify_ltsep_preservation', return_value={})
                reason = 'required_external_input_missing:sigma0_SIMC' if sigma_case == 'invalid' else 'analysis_failed'
                with self.assertRaisesRegex(ValueError, reason):
                    owner.execute_gate(repo, shared, globus / 'new.zip', head)
                inventory_call.assert_called_once()
                self.assertEqual(inventory_call.call_args.args[3], sigma_env)
                if sigma_case == 'invalid':
                    child.assert_not_called()
                else:
                    child.assert_called_once()
                    child_env = child.call_args.args[-1]
                    self.assertEqual(child_env[owner.ab.ENABLE_ENV], child_env[cap.ENABLE_ENV])
                    self.assertEqual({k: child_env[k] for k in sigma_env}, sigma_env)
                collect.assert_not_called()
                self.assertGreaterEqual(checkout.call_count, 1)
                self.assertGreaterEqual(installed.call_count, 1)
            self.assertIn('analysis', cleaned)
            self.assertEqual({k: v for k, v in owner.inventory(shared).items() if k != 'new.log'}, original)
            self.assertEqual((shared / 'new.log').is_file(), sigma_case != 'invalid')

    def test_all_preservation_checks_run_after_first_failure(self):
        calls = []
        def fail():
            calls.append('shared')
            raise ValueError('changed')
        with self.assertRaisesRegex(ValueError, 'shared:changed'):
            owner.preservation_checks([('shared', fail), ('installed', lambda: calls.append('installed')),
                                       ('checkout', lambda: calls.append('checkout'))])
        self.assertEqual(calls, ['shared', 'installed', 'checkout'])

    def test_private_link_and_collision_preserve_external(self):
        with tempfile.TemporaryDirectory() as d:
            parent = Path(d).resolve()
            worktree = parent / 'analysis'
            worktree.mkdir()
            volatile = parent / 'external'
            volatile.mkdir()
            sentinel = volatile / 'ordinary.bin'
            sentinel.write_bytes(b'unchanged')
            paths = {'LTANAPATH': str(worktree), 'ANATYPE': 'Kaon',
                     'OUTPATH': str(worktree / 'OUTPUT/Analysis/KaonLT'), 'VOLATILEPATH': str(volatile)}
            output = owner.prepare_private_output(worktree, paths)
            self.assertEqual((worktree / 'OUTPUT/Analysis/KaonLT').resolve(), output)
            self.assertFalse(output.is_relative_to(volatile))
            with self.assertRaisesRegex(ValueError, 'collision'):
                owner.prepare_private_output(worktree, paths)
            self.assertEqual(sentinel.read_bytes(), b'unchanged')

    def test_private_copies_preserve_input_inode_and_mutation_is_local(self):
        with tempfile.TemporaryDirectory() as d:
            shared, private = Path(d) / 'shared', Path(d) / 'private'
            shared.mkdir(); private.mkdir()
            for setting in owner.SETTINGS:
                for role in ('Data', 'Dummy'):
                    (shared / f'{setting}_kaon_Analysed_{role}_Q4p4W2p74_lowe.root').write_bytes((setting + role).encode())
                for suffix in ('.root', '.hist'):
                    (shared / f'Prod_Coin_Q4p4W2p74{setting.lower()}_lowe{suffix}').write_bytes(b'simc')
            before = owner.inventory(shared)
            bindings = owner.stage_input_files(shared, private)
            self.assertEqual(len(bindings), 8)
            target = next(iter(bindings.values()))
            Path(target['private']).write_bytes(b'private mutation')
            self.assertEqual(owner.inventory(shared), before)
            with self.assertRaisesRegex(ValueError, 'collision'):
                owner.stage_input_files(shared, private)

    def test_incomplete_high_bundle_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            shared, private = Path(d) / 'shared', Path(d) / 'private'
            shared.mkdir(); private.mkdir()
            (shared / 'Center_kaon_Analysed_Data_Q4p4W2p74_lowe.root').write_bytes(b'data')
            with self.assertRaisesRegex(ValueError, 'incomplete_data_dummy_pair'):
                owner.stage_input_files(shared, private)

    def test_diamond_mid_and_fallback_input_namespace_is_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            shared, private = Path(d) / 'shared', Path(d) / 'private'
            shared.mkdir(); private.mkdir()
            for setting in owner.SETTINGS:
                for role in ('Data', 'Dummy'):
                    (shared / f'{setting}_kaon_Analysed_{role}_Q4p4W2p74_lowe.root').write_bytes((setting+role).encode())
                for suffix in ('.root', '.hist'):
                    (shared / f'Prod_Coin_Q4p4W2p74{setting.lower()}_lowe{suffix}').write_bytes(b'simc')
            extra = shared / 'Center_kaon_Data_Q4p4W2p74_mide.root'
            extra.write_bytes(b'independent diamond candidate')
            records = owner.stage_input_files(shared, private)
            self.assertIn(extra.name, records)
            self.assertEqual(owner.diamond_selections(shared), owner.diamond_selections(private))
            self.assertEqual((private / extra.name).read_bytes(), extra.read_bytes())

    def test_two_observers_only_in_owned_environment(self):
        with tempfile.TemporaryDirectory() as d:
            manifest = Path(d) / 'manifest.json'
            manifest.write_bytes(b'{}')
            source = {'PATH': 'preserved', 'LT_ANALYSIS_DEBUG_LEFT_LOW': '1'}
            env = owner.child_environment(source, manifest)
            self.assertEqual(env[owner.ab.ENABLE_ENV], '1')
            self.assertEqual(env[cap.ENABLE_ENV], '1')
            self.assertNotIn('LT_ANALYSIS_DEBUG_LEFT_LOW', env)
            self.assertEqual(source['LT_ANALYSIS_DEBUG_LEFT_LOW'], '1')
            with self.assertRaisesRegex(ValueError, 'untrusted_analysis_environment'):
                owner.child_environment({'LT_ANALYSIS_PION_LOCAL_T_EVENT_VALIDATION': '1'}, manifest)

    def test_default_writer_source_is_reached_without_plot_mode(self):
        report = owner.audit_default_output_writer(ROOT)
        self.assertEqual(report['mode'], [])
        self.assertIn('table_kaon_', report['destination'])

    def test_no_child_before_path_gates_and_no_retry(self):
        import inspect
        source = inspect.getsource(owner.execute_gate)
        self.assertLess(source.index('probe_private_callers('), source.index('run_owned_child('))
        self.assertLess(source.index('stage_input_files('), source.index('run_owned_child('))
        self.assertLess(source.index('verify_joint('), source.index('module.collect_validation_bundle('))
        self.assertEqual(source.count('run_owned_child('), 1)

    def test_wrong_mode_destination_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            tree = Path(d).resolve()
            base = {'VOLATILEPATH': str(tree / 'external'), 'OUTPATH': str(tree / 'external'),
                    'ANATYPE': 'Kaon', 'LTANAPATH': str(tree)}
            record = {'paths': base, 'package_file': str(tree / 'runtime-python/ltsep/__init__.py')}
            with patch.object(owner, 'reached_resolver_callers', return_value=[]), patch.object(owner.isolation, 'probe_paths', return_value=[]), patch.object(owner.resolver, 'probe', return_value=record):
                with self.assertRaisesRegex(ValueError, 'default_root_private_outpath_mismatch'):
                    owner.probe_private_callers(tree, tree, {}, base, tree / 'private')

    def test_reached_writers_cover_main_and_setup_contexts(self):
        callers = owner.reached_resolver_callers(ROOT)
        self.assertIn(('src/normalize/get_eff_charge.py', ('Plot_LTSep',), 'src'), callers)
        self.assertIn(('src/binning/calculate_yield.py', ('Plot_LTSep',), 'src'), callers)
        self.assertIn(('src/setup/findEffectiveCharge.py', (), 'src/setup'), callers)
        self.assertIn(('src/utility/utility.py', ('Plot_LTSep',), 'src/setup'), callers)
        self.assertIn(('src/utility/utility.py', ('Plot_LTSep',), 'src'), callers)
        self.assertFalse(any('main_iter.py' in row[0] or '/HeeP/' in row[0] for row in callers))

    def test_observer_open_audit_detects_duplicate_and_escape(self):
        with tempfile.TemporaryDirectory() as d:
            directory = Path(d)
            manifest = directory / 'matched-cohort-authority.json'
            manifest.write_bytes(cap.encoded({'private_outpath': str(directory)}))
            env = owner.child_environment(os.environ, manifest)
            name = owner.comparison_name('Center')
            script = owner.OBSERVER_AUDIT + '\np = output / ' + repr(name) + '\np.write_text("first")\np.write_text("final")\n'
            run = subprocess.run([sys.executable, '-B', '-c', script], env=env, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            run = subprocess.run([sys.executable, '-B', '-c', script + 'p.write_text("overwrite")'], env=env, capture_output=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn(b'observer_output_overwrite', run.stderr)
            escaped = owner.OBSERVER_AUDIT + '\n(pathlib.Path(output.parent) / ' + repr(name) + ').write_text("escape")'
            run = subprocess.run([sys.executable, '-B', '-c', escaped], env=env, capture_output=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn(b'observer_output_escape', run.stderr)

    def test_matched_authority_requires_fresh_AB_and_exact_parents(self):
        inputs = {'data': {'sha256': 'a'*64, 'bytes': 10}, 'dummy': {'sha256': 'd'*64, 'bytes': 12}}
        parents = [dict(t_bin_index=i, t_edges=[i, i+1], pion_parent_id='p'+str(i),
                        canonical_interval_pair_id='pair', canonical_interval_pair_hash='f'*64,
                        coordinate_fingerprint='c'*64, baseline_parent_fingerprint='b'*64,
                        pion_input={'synthetic': 1}, kaon_input={'synthetic': 2}) for i in range(3)]
        packet = {'source_HEAD': 'a'*40, 'input_identity': inputs, 'parents': parents}
        rows = [{'A_B_unchanged': dict(p, frozen_inputs={'pion': p['pion_input'], 'kaon': p['kaon_input']})}
                for p in parents]
        with tempfile.TemporaryDirectory() as d:
            root = Path(d).resolve() / 'analysis'
            (root / 'testing').mkdir(parents=True)
            (root / owner.OWNER).write_text('synthetic owner')
            output = root.parent / 'matched-output/Analysis/KaonLT'
            output.mkdir(parents=True)
            manifest = {'schema_version': 'pion_matched_cohort_authority/v1', 'mode': 'matched_cohort',
                'source_HEAD': packet['source_HEAD'], 'inputs': {s: packet['input_identity'] for s in owner.SETTINGS},
                'private_outpath': str(output), 'owner_source': cap.file_identity(root / owner.OWNER),
                'started_ns': 0, 'attempt_id': 'synthetic'}
            document = {'source_HEAD': packet['source_HEAD'], 'status': 'available', 'diagnostic_only': True,
                        'production_applied': False, 'phi_setting': 'Center',
                        'original_analysis_input_identity': packet['input_identity'],
                        'parents': [row['A_B_unchanged'] for row in rows]}
            abpath = output / owner.comparison_name('Center')
            abpath.write_bytes(cap.encoded(document))
            manifest_path = output / 'matched-cohort-authority.json'
            manifest_path.write_bytes(cap.encoded(manifest))
            env = {cap.MATCHED_ENV: str(manifest_path), cap.MATCHED_HASH_ENV: cap.digest(manifest)}
            with patch.dict(os.environ, env):
                auth = cap.matched_authority(root, output, 'Center', packet['input_identity'], packet['parents'], packet['source_HEAD'])
                self.assertEqual(auth['AB'], cap.file_identity(abpath))
                for key in ('pion_parent_id', 'coordinate_fingerprint', 'canonical_interval_pair_hash'):
                    bad = deepcopy(packet['parents'])
                    bad[0][key] = 'wrong'
                    with self.assertRaisesRegex(ValueError, 'matched_AB_parent_identity'):
                        cap.matched_authority(root, output, 'Center', packet['input_identity'], bad, packet['source_HEAD'])
                wrong = deepcopy(packet['input_identity'])
                wrong['dummy']['bytes'] += 1
                with self.assertRaisesRegex(ValueError, 'matched_input_identity'):
                    cap.matched_authority(root, output, 'Center', wrong, packet['parents'], packet['source_HEAD'])
                abpath.unlink()
                with self.assertRaisesRegex(ValueError, 'matched_AB_freshness'):
                    cap.matched_authority(root, output, 'Center', packet['input_identity'], packet['parents'], packet['source_HEAD'])

    def test_original_and_rewritten_runtime_container_identities_are_distinct(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d).resolve() / 'analysis'
            output = root.parent / 'matched-output/Analysis/KaonLT'
            (root / 'testing').mkdir(parents=True)
            output.mkdir(parents=True)
            (root / owner.OWNER).write_text('synthetic owner')
            writers = {}
            for name in ('src/cuts/shift_prep.py', 'src/setup/shift_MM.py'):
                p = root / name
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text('# frozen writer fixture')
                writers[name] = cap.file_identity(p)
            original, runtime, bindings = {}, {}, {}
            for role, title in (('data', 'Data'), ('dummy', 'Dummy')):
                source = root.parent / (role + '.bin')
                source.write_bytes(('original ' + role).encode())
                name = f'Center_kaon_Analysed_{title}_Q4p4W2p74_lowe.root'
                private = output / name
                private.write_bytes(('runtime with derived branches ' + role).encode())
                original[role], runtime[role] = cap.file_identity(source), cap.file_identity(private)
                bindings[name] = {'source': str(source), 'original': owner.identity(source),
                                  'private': str(private), 'private_mutation': 'frozen_shiftprep_derived_branches'}
            parents = [dict(t_bin_index=i, t_edges=[i, i+1], pion_parent_id='p'+str(i),
                            canonical_interval_pair_id='pair', canonical_interval_pair_hash='f'*64,
                            coordinate_fingerprint='c'*64, baseline_parent_fingerprint='b'*64,
                            pion_input={}, kaon_input={}) for i in range(3)]
            authority = {'schema_version': 'pion_matched_cohort_authority/v1', 'mode': 'matched_cohort',
                'source_HEAD': 'a'*40, 'started_ns': 0, 'private_outpath': str(output), 'attempt_id': 'synthetic',
                'inputs': {'Center': original, 'Left': original}, 'owner_source': cap.file_identity(root / owner.OWNER),
                'input_bindings': bindings, 'private_input_policy': 'frozen_shiftprep_private_copy',
                'private_writer_source_files': writers}
            manifest = output / 'matched-cohort-authority.json'
            manifest.write_bytes(cap.encoded(authority))
            document = {'source_HEAD': 'a'*40, 'phi_setting': 'Center', 'status': 'available',
                'diagnostic_only': True, 'production_applied': False, 'original_analysis_input_identity': runtime,
                'parents': [dict(p, frozen_inputs={'pion': {}, 'kaon': {}}) for p in parents]}
            (output / owner.comparison_name('Center')).write_bytes(cap.encoded(document))
            with patch.dict(os.environ, {cap.MATCHED_ENV: str(manifest), cap.MATCHED_HASH_ENV: cap.digest(authority)}):
                result = cap.matched_authority(root, output, 'Center', runtime, parents, 'a'*40)
                self.assertEqual(result['original_input_identity'], original)
                self.assertEqual(result['runtime_input_identity'], runtime)
                self.assertNotEqual(original, runtime)
                Path(bindings[next(iter(bindings))]['source']).write_bytes(b'changed shared input')
                with self.assertRaisesRegex(ValueError, 'matched_original_input_changed'):
                    cap.matched_authority(root, output, 'Center', runtime, parents, 'a'*40)


if __name__ == '__main__':
    unittest.main()
