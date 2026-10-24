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
import errno
import shutil
import venv
from array import array
import unittest
from unittest.mock import patch
from unittest.mock import Mock
from types import SimpleNamespace

from testing import run_pion_local_t_matched_cohort_gate as owner
from testing import pion_local_t_event_validation_capture as cap

ROOT = Path(__file__).resolve().parents[1]


class ResolverImports(unittest.TestCase):
    """Real subprocess imports in an isolated interpreter, never fake Root results."""
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix='kaonlt-import-interpreter-')
        cls.venv = Path(cls.temp.name).resolve()
        venv.EnvBuilder(with_pip=False).create(cls.venv)
        cls.python = cls.venv / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
        cls.site = Path(subprocess.check_output([str(cls.python), '-c',
            'import sysconfig; print(sysconfig.get_path("purelib"))'], text=True).strip())

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    @contextmanager
    def fixture(self, behavior='', dependency_behavior='', smoke_cleanup=False, real_dependencies=False, reached_callers=False):
        with tempfile.TemporaryDirectory(prefix='kaonlt-import-fixture-') as d:
            root = Path(d).resolve()
            repo, worktree = root / 'ordinary', root / 'attempt/analysis'
            installed, copied = root / 'installed/ltsep', root / 'attempt/runtime-python/ltsep'
            fields = {k: str(root / 'external' / k) for k in owner.resolver.FIELDS}
            fields.update(USER='fixture', HOST='fixture', ANATYPE='Kaon', LTANAPATH=str(repo),
                          OUTPATH=str(repo / 'OUTPUT/Analysis/KaonLT'))
            # Replay the supplied farm geometry under a temporary root. No-mode
            # follows stable hallc_replay_lt/UTIL_KAONLT, Plot follows LTANAPATH.
            fields['UTILPATH'] = str(root / 'hallc_replay_lt/UTIL_KAONLT')
            reached = owner.reached_resolver_callers(ROOT) if smoke_cleanup or reached_callers else []
            for tree in (repo, worktree):
                for name in {s for _, _, _, _, s in owner.resolver.CALLERS} | {
                        'src/setup/calcTotalEffectiveCharge.py'} | {c for c, _, _ in reached}:
                    p = tree / name
                    p.parent.mkdir(parents=True, exist_ok=True)
                    p.write_text('# caller must never execute\nraise RuntimeError("caller_executed")\n')
                for name in ('findEffectiveCharge.py', 'getEfficiencyValue.py', 'getDataTable.py',
                             'calcTotalEffectiveCharge.py'):
                    shutil.copyfile(ROOT / 'src/setup' / name, tree / 'src/setup' / name)
                (tree / 'farm_env').mkdir()
                for name in ('ltsep_paths.py',):
                    shutil.copyfile(ROOT / 'farm_env' / name, tree / 'farm_env' / name)
                if tree == repo:
                    (tree / 'OUTPUT/Analysis/KaonLT').mkdir(parents=True)
            worktree.parent.chmod(0o700)
            private = worktree.parent / 'matched-output/Analysis/KaonLT'
            if not smoke_cleanup:
                private = owner.prepare_private_output(worktree, dict(fields, LTANAPATH=str(worktree),
                    OUTPATH=str(worktree / 'OUTPUT/Analysis/KaonLT')))
            util_output = Path(fields['UTILPATH']) / 'OUTPUT'
            util_output.parent.mkdir(parents=True)
            util_output.symlink_to(repo / 'OUTPUT', target_is_directory=True)
            for package, tree in ((installed, repo), (copied, worktree)):
                (package / 'PATH_TO_DIR').mkdir(parents=True)
                (package / '__init__.py').write_text('from .ltsep import Root\n')
                (package / 'ltsep.py').write_text(
                    'from pathlib import Path\n' + (dependency_behavior if real_dependencies else 'import uproot\n') + 'from .pathing import SetPath\n'
                    'class Root:\n'
                    ' def __init__(self,caller,*modes):\n'
                    '  p=SetPath(caller)\n'
                    '  for k in p.values: setattr(self,k,p.getPath(k))\n'
                    '  if modes not in ((), ("Plot_LTSep",)): raise ValueError("wrong mode")\n'
                    '  self.OUTPATH=str(Path(self.LTANAPATH if modes else self.UTILPATH)/"OUTPUT/Analysis"/("KaonLT" if modes else "None"))\n'
                    + behavior)
                (package / 'pathing.py').write_text(
                    'from pathlib import Path\nclass SetPath:\n'
                    ' def __init__(self,caller):\n'
                    '  self.values=dict(line.split("=",1) for line in (Path(__file__).parent/"PATH_TO_DIR/user.path").read_text().splitlines())\n'
                    ' def getPath(self,key): return self.values[key]\n')
                config = dict(fields, LTANAPATH=str(tree))
                (package / 'PATH_TO_DIR/user.path').write_text(''.join(k+'='+v+'\n' for k,v in config.items()))
            # These are actual interpreter site-packages, outside ltsep/stdlib.
            for name, text in (() if real_dependencies else (('uproot', 'import numpy\n' + dependency_behavior), ('numpy', 'fixture_dependency=True\n'))):
                p = self.site / name
                p.mkdir()
                (p / '__init__.py').write_text(text)
            (root / 'physics.root').write_bytes(b'never read by resolver')
            for suffix in ('root', 'csv', 'json', 'dat'):
                (self.site / ('physics.' + suffix)).write_bytes(b'never read even inside site-packages')
            baseline_env = dict(os.environ, PYTHONPATH=str(installed.parent), PYTHONDONTWRITEBYTECODE='1')
            private_env = dict(baseline_env, PYTHONPATH=str(copied.parent))
            before = {str(p): p.read_bytes() for p in root.rglob('*') if p.is_file()}
            try:
                with patch.object(owner.sys, 'executable', sys.executable if real_dependencies else str(self.python)), patch.dict(os.environ, baseline_env, clear=True):
                    yield SimpleNamespace(root=root, repo=repo, tree=worktree, package=installed,
                        copied=copied, baseline=fields, env=private_env,
                        outdir=private)
            finally:
                for name, raw in before.items():
                    if smoke_cleanup and Path(name).is_relative_to(worktree.parent):
                        self.assertFalse(Path(name).exists())
                        continue
                    self.assertEqual(Path(name).read_bytes(), raw, name)
                self.assertFalse((root / 'written').exists())
                for name in (() if real_dependencies else ('uproot', 'numpy')):
                    shutil.rmtree(self.site / name)
                for p in self.site.glob('physics.*'):
                    p.unlink()

    def test_smoke_real_resolver_and_no_child_success_failure(self):
        # Only isolation setup is a fixture. The shared preflight and every
        # ordinary baseline/private Root subprocess run unchanged.
        for behavior, error in (('', None),
                ('mutate_installed', 'installed_ltsep_package_changed'),
                ('mutate_copied', 'copied_ltsep_package_changed'),
                ('mutate_index', 'smoke_worktree_source_index_changed'),
                ('mutate_link', 'external_link_changed'),
                ('collision', 'private_output_collision'),
                ('  self.OUTPATH="/unexpected/output"\n', 'outpath_mismatch'),
                ('  raise RuntimeError("dependency_import_failed")\n', 'dependency_import_failed')):
            with self.subTest(error=error), self.fixture(
                    behavior='' if behavior.startswith('mutate_') or behavior == 'collision' else behavior,
                    smoke_cleanup=True) as f:
                f.tree.parent.chmod(0o700)
                shared = Path(f.baseline['VOLATILEPATH']) / 'OUTPUT/Analysis/KaonLT'
                (shared / 'background_samples').mkdir(parents=True)
                f.baseline['OUTPATH'] = str(shared)
                # Baseline lexical OUTPATH is the ordinary OUTPUT alias.
                shutil.rmtree(f.repo / 'OUTPUT')
                (f.repo / 'OUTPUT').symlink_to(shared.parents[1], target_is_directory=True)
                private_paths = dict(f.baseline, LTANAPATH=str(f.tree),
                                     OUTPATH=str(f.tree / 'OUTPUT/Analysis/KaonLT'))
                head = 'c' * 40
                state = dict(head=head, origin_test=head, branch='test', porcelain='')
                cleaned = []
                @contextmanager
                def owned(*args):
                    try:
                        yield f.tree
                    finally:
                        shutil.rmtree(f.tree.parent)
                        cleaned.append(True)
                with ExitStack() as stack:
                    def patched(obj, name, **kwargs):
                        return stack.enter_context(patch.object(obj, name, **kwargs))
                    patched(owner.isolation.accepted, 'REPO', new=f.repo)
                    patched(owner.isolation.accepted, 'ARTIFACTS', new=shared)
                    patched(owner.isolation, 'snapshot', return_value=state)
                    patched(owner, 'smoke_source_state', return_value=state)
                    patched(owner.ab, 'git', return_value='')
                    patched(owner.ab, 'sha256', return_value='fixture-index')
                    patched(owner.shutil, 'which', return_value=str(self.python))
                    patched(owner, 'audit_default_output_writer', return_value={})
                    patched(owner.isolation, 'owned_worktree', new=owned)
                    patched(owner.isolation, 'import_identity', return_value={
                        'package_file': str(f.package / '__init__.py')})
                    patched(owner.isolation, 'prepare_runtime_overlay', return_value=(f.env, {
                        'baseline_paths': f.baseline, 'original_source_sha256': {}}))
                    patched(owner.isolation, 'probe_paths', return_value=[{'paths': private_paths}])
                    link = f.root / 'external-link'
                    link.symlink_to(shared, target_is_directory=True)
                    patched(owner.isolation, 'external_symlink_preflight', return_value={
                        'links': [{'path': str(link), 'literal_target': str(shared),
                                   'resolved_target': str(shared)}], 'background_simc_path': ''})
                    patched(owner, 'skim_paths', return_value={'launcher': str(shared)})
                    reached = owner.reached_resolver_callers(ROOT) if error is None else []
                    patched(owner, 'reached_resolver_callers', return_value=reached)
                    forbidden = [patched(owner, name, side_effect=AssertionError('forbidden:' + name))
                        for name in ('stage_input_files', 'external_input_inventory', 'inventory',
                                     'run_owned_child', 'verify_joint', 'verify_archive')]
                    forbidden += [patched(cap, 'file_identity',
                        side_effect=AssertionError('input artifact hashing forbidden')),
                        patched(owner.ab, 'GateStatus', side_effect=AssertionError('log/status forbidden')),
                        patched(owner.zipfile, 'ZipFile', side_effect=AssertionError('ZIP forbidden'))]
                    forbidden += [patched(owner.isolation, name,
                        side_effect=AssertionError('forbidden:' + name))
                        for name in ('collection_module', 'deliver_evidence')]
                    if behavior == 'collision':
                        (f.tree.parent / 'matched-output').mkdir()
                    real_preflight = owner.resolver_preflight
                    original_package = (f.package / '__init__.py').read_bytes()
                    def preflight(*args, **kwargs):
                        result = real_preflight(*args, **kwargs)
                        if behavior in ('mutate_installed', 'mutate_copied'):
                            target = f.package if behavior == 'mutate_installed' else f.copied
                            (target / '__init__.py').write_text('# mutated')
                        if behavior == 'mutate_link':
                            link.unlink()
                            link.symlink_to(f.root, target_is_directory=True)
                        if behavior == 'mutate_index':
                            owner.ab.sha256.return_value = 'mutated-index'
                        return result
                    patched(owner, 'resolver_preflight', side_effect=preflight)
                    real_probe = patched(owner, 'probe_resolver', wraps=owner.probe_resolver)
                    if error:
                        with self.assertRaisesRegex(ValueError, error):
                            owner.resolver_smoke(f.repo, head)
                    else:
                        result = owner.resolver_smoke(f.repo, head)
                        self.assertTrue(result['cleanup_completed'])
                        self.assertTrue(result['preservation_passed'])
                        self.assertFalse(result['analysis_started'])
                        self.assertFalse(result['archive_created'])
                        self.assertTrue(result['runtime_not_science'])
                        known = {(c, tuple(m)) for _, c, m, _, _ in owner.resolver.CALLERS}
                        expected = len(owner.resolver.CALLERS) + sum(
                            (c, tuple(m)) not in known for c, m, _ in reached)
                        self.assertEqual(real_probe.call_count, 2 * expected)
                        for record in result['probes']['all_callers'] + result['probes']['baseline_callers']:
                            self.assertIn('uproot', record['source_imports'])
                            self.assertIn('numpy', record['source_imports'])
                    (f.package / '__init__.py').write_bytes(original_package)
                    for mock in forbidden:
                        mock.assert_not_called()
                self.assertEqual(cleaned, [True])
                self.assertEqual(list(shared.iterdir()), [shared / 'background_samples'])

    def test_actual_resolver_rejects_wrong_identity_receipts(self):
        for key, reason in (('cwd', 'resolver_cwd_mismatch'),
                            ('caller', 'resolver_caller_mismatch'),
                            ('caller_realpath', 'resolver_caller_mismatch'),
                            ('mode', 'resolver_caller_mismatch'),
                            ('python', 'interpreter_mismatch'),
                            ('resolver_file', 'resolver_wrapper_import_mismatch'),
                            ('pathing_file', 'ltsep_module_identity_invalid'),
                            ('root_file', 'ltsep_module_identity_invalid')):
            with self.subTest(key=key), self.fixture() as f:
                run = subprocess.run
                def wrong_identity(*args, **kwargs):
                    result = run(*args, **kwargs)
                    if result.returncode == 0:
                        record = json.loads(result.stdout)
                        record[key] = ['Plot_LTSep'] if key == 'mode' else str(f.root)
                        result.stdout = json.dumps(record)
                    return result
                with patch.object(owner.subprocess, 'run', side_effect=wrong_identity):
                    with self.assertRaisesRegex(ValueError, reason):
                        owner.probe_resolver(f.tree, 'src/setup/findEffectiveCharge.py', (),
                                             f.env, 'src/setup', f.copied)

    def test_group_u_group_aliases_across_baseline_private_and_inline_probe(self):
        # Real POSIX root alias, actual subprocess/cwd/imports, both modes and
        # all reached callers. Python 3.12 path-shape evidence, not farm 3.9.
        with self.fixture(reached_callers=True) as f:
            physical = f.root / 'u/group'
            physical.parent.mkdir()
            physical.symlink_to(f.root, target_is_directory=True)
            group = f.root / 'group'
            group.symlink_to(physical, target_is_directory=True)
            repo, tree = group / 'ordinary', group / 'attempt/analysis'
            configs = [f.package / 'PATH_TO_DIR/user.path', f.copied / 'PATH_TO_DIR/user.path']
            originals = [p.read_bytes() for p in configs]
            try:
                for config, old, alias in zip(configs, (f.repo, f.tree), (repo, tree)):
                    config.write_bytes(config.read_bytes().replace(str(old).encode(), str(alias).encode()))
                env = dict(f.env, PYTHONPATH=str(group / 'attempt/runtime-python'))
                baseline_env = dict(os.environ, PYTHONPATH=str(group / 'installed'))
                with patch.dict(os.environ, baseline_env, clear=True), patch.object(
                        owner, 'reached_resolver_callers', return_value=owner.reached_resolver_callers(ROOT)):
                    result = owner.probe_private_callers(repo, tree, env, f.baseline,
                        f.outdir, group / 'installed/ltsep')
                for record in result['all_callers'] + result['baseline_callers']:
                    self.assertNotEqual(record['caller'], record['caller_realpath'])
                    self.assertEqual(Path(record['caller']).resolve(), Path(record['caller_realpath']))
                    self.assertEqual(Path(record['cwd']).resolve(), Path(record['cwd']))
            finally:
                for config, raw in zip(configs, originals):
                    config.write_bytes(raw)

    def test_missing_caller_and_redirected_caller_cwd_fail_before_import_or_hash(self):
        with self.fixture() as f:
            for tree in (f.repo, f.tree):
                path = tree / 'src/setup/findEffectiveCharge.py'
                raw = path.read_bytes()
                path.unlink()
                with patch.object(owner.subprocess, 'run') as run, patch.object(owner.resolver, 'tree_bytes') as hashes:
                    with self.assertRaisesRegex(ValueError, 'missing_or_nonregular:.*canonical.*lexical'):
                        owner.probe_resolver(tree, 'src/setup/findEffectiveCharge.py', (), f.env, 'src/setup', f.copied)
                    outside = f.root / 'wrong-checkout.py'
                    outside.write_bytes(raw)
                    path.symlink_to(outside)
                    with self.assertRaisesRegex(ValueError, 'alias_mismatch:.*authorized.*canonical.*lexical'):
                        owner.probe_resolver(tree, 'src/setup/findEffectiveCharge.py', (), f.env, 'src/setup', f.copied)
                    path.unlink()
                    path.write_bytes(raw)
                    cwd = tree / 'redirected-cwd'
                    cwd.symlink_to(f.root, target_is_directory=True)
                    with self.assertRaisesRegex(ValueError, 'cwd_alias_mismatch'):
                        owner.probe_resolver(tree, 'src/setup/findEffectiveCharge.py', (), f.env, 'redirected-cwd', f.copied)
                    run.assert_not_called()
                    hashes.assert_not_called()

    def test_ordinary_import_package_mutation_stops_on_success_and_exception(self):
        for installed in (False, True):
            for fails in (False, True):
                behavior = ('  Path(__file__).with_name("__init__.py").write_text("# changed")\n' +
                            ('  raise RuntimeError("import_failed")\n' if fails else ''))
                with self.subTest(installed=installed, fails=fails), self.fixture(behavior=behavior) as f:
                    package, tree, env = ((f.package, f.repo, dict(os.environ)) if installed else
                                          (f.copied, f.tree, f.env))
                    path = package / '__init__.py'
                    raw = path.read_bytes()
                    try:
                        with self.assertRaisesRegex(ValueError, 'ltsep_package_changed_during_probe'):
                            owner.probe_resolver(tree, 'src/setup/findEffectiveCharge.py', (),
                                                 env, 'src/setup', package)
                    finally:
                        path.write_bytes(raw)

    def test_ordinary_import_ctypes_and_legitimate_timezone_read(self):
        zone = Path('/usr/share/zoneinfo/UTC')
        behavior = 'import ctypes\nctypes.pythonapi.Py_GetVersion\n'
        if zone.is_file():
            behavior += 'from pathlib import Path\nassert Path("/usr/share/zoneinfo/UTC").read_bytes()\n'
        with self.fixture(dependency_behavior=behavior) as f:
            # The environment must not re-enable bytecode writes.
            f.env['PYTHONDONTWRITEBYTECODE'] = '0'
            self.run_private(f)
            self.assertFalse(list(f.copied.rglob('*.pyc')))

    def run_private(self, f):
        return owner.probe_private_callers(f.repo, f.tree, f.env, f.baseline, f.outdir, f.package)

    def test_real_import_closure_and_all_caller_modes(self):
        with self.fixture() as f:
            result = self.run_private(f)
            known = {(caller, tuple(modes)) for _, caller, modes, _, _ in owner.resolver.CALLERS}
            expected = len(owner.resolver.CALLERS) + sum(
                (caller, tuple(modes)) not in known for caller, modes, _ in owner.reached_resolver_callers(f.repo))
            self.assertEqual(len(result['all_callers']), expected)
            self.assertTrue(result['plot'])
            for record in result['all_callers']:
                self.assertEqual(record['cwd'], str(f.tree / ('src/setup' if not record['mode'] else
                    next(cwd for _, caller, modes, cwd, _ in owner.resolver.CALLERS
                         if str(f.tree / caller) == record['caller'] and list(modes) == record['mode']))))
                self.assertEqual(owner.resolver.module_identity(record), f.copied)
                self.assertEqual(Path(record['paths']['OUTPATH']).name, 'KaonLT' if record['mode'] else 'None')
                self.assertTrue({'ltsep', 'ltsep.ltsep', 'uproot', 'numpy'}.issubset(record['source_imports']))

    def test_unexpected_package_identity_fails(self):
        with self.fixture() as f:
            f.env['PYTHONPATH'] = str(f.package.parent)
            with self.assertRaisesRegex(ValueError, 'resolver_package_import_mismatch'):
                self.run_private(f)

    def test_ordinary_dependencies_and_timezone_reads(self):
        import importlib.util
        names = ('numpy', 'pandas', 'dateutil')
        missing = [name for name in names if importlib.util.find_spec(name) is None]
        if missing:
            self.skipTest('NOT VERIFIED: real dependencies unavailable: ' + ', '.join(missing))
        # Real host packages, never the isolated venv's dummy numpy/uproot.
        with self.fixture(real_dependencies=True, reached_callers=True, dependency_behavior=(
                'import numpy, pandas, dateutil, ctypes\n'
                'from dateutil import tz\nassert tz.gettz("UTC") is not None\n'
                'ctypes.pythonapi.Py_GetVersion\n')) as f:
            reached = owner.reached_resolver_callers(ROOT)
            with patch.object(owner, 'reached_resolver_callers', return_value=reached):
                result = self.run_private(f)
            known = {(c, tuple(m)) for _, c, m, _, _ in owner.resolver.CALLERS}
            expected = len(owner.resolver.CALLERS) + sum((c, tuple(m)) not in known for c, m, _ in reached)
            self.assertEqual(len(result['all_callers']), expected)
            self.assertEqual(len(result['baseline_callers']), expected)
            for record in result['all_callers'] + result['baseline_callers']:
                self.assertTrue(set(names + ('ctypes',)).issubset(record['source_imports']))
                self.assertNotIn(str(self.site), record['source_imports']['numpy']['path'])
            self.assertFalse(list(f.package.rglob('*.pyc')))
            self.assertFalse(list(f.copied.rglob('*.pyc')))

    def test_root_semantics_and_unique_config_fail_closed(self):
        cases = [
            ('  if not modes: self.OUTPATH=str(Path(self.LTANAPATH)/"wrong")\n', 'root_outpath_mismatch'),
            ('  if modes: self.OUTPATH=str(Path(self.LTANAPATH)/"OUTPUT/Analysis/None")\n', 'root_outpath_mismatch'),
            ('  self.UTILPATH="changed"\n', 'baseline_input_role_changed|stable_input_role_changed'),
        ]
        for behavior, message in cases:
            with self.subTest(behavior=behavior), self.fixture(behavior=behavior) as f:
                with self.assertRaisesRegex(ValueError, message):
                    self.run_private(f)
        for duplicate in (False, True):
            with self.fixture() as f:
                config = f.copied / 'PATH_TO_DIR/user.path'
                if duplicate:
                    second = config.with_name('second.path')
                    shutil.copyfile(config, second)
                    try:
                        with self.assertRaisesRegex(ValueError, 'matching_config_not_unique'):
                            self.run_private(f)
                    finally:
                        second.unlink()
                else:
                    # Valid parsed config with no matching LTANAPATH, not a missing file.
                    original = config.read_bytes()
                    config.write_bytes(original.replace(b'LTANAPATH=', b'OTHERPATH='))
                    try:
                        with self.assertRaisesRegex(ValueError, 'resolver_subprocess_failed|matching_config_not_unique'):
                            self.run_private(f)
                    finally:
                        config.write_bytes(original)

    def test_private_only_role_and_output_changes_fail(self):
        for statement, reason in (
                ('self.UTILPATH="changed"', 'stable_input_role_changed'),
                ('self.OUTPATH=str(Path(self.LTANAPATH)/"wrong")', 'private_outpath_mismatch')):
            with self.fixture(behavior='  if "runtime-python" in __file__: '+statement+'\n') as f:
                with self.assertRaisesRegex(ValueError, reason):
                    self.run_private(f)

    def test_package_and_config_realpath_changes_fail(self):
        for config in (False, True):
            with self.subTest(config=config), self.fixture() as f:
                path = f.copied / ('PATH_TO_DIR/user.path' if config else 'pathing.py')
                raw = path.read_bytes()
                target = f.root / 'escaped'
                target.write_bytes(raw)
                path.unlink()
                path.symlink_to(target)
                try:
                    with self.assertRaisesRegex(ValueError, 'package_symlink_not_supported|resolver_path_escape'):
                        self.run_private(f)
                finally:
                    path.unlink()
                    path.write_bytes(raw)

    def test_real_probe_success_and_failure_enter_owner_child_barrier(self):
        # Keep unrelated orchestration fixtures synthetic, but execute the full
        # subprocess import/caller gate rather than mocking probe return values.
        real_probe = owner.probe_private_callers
        for failure in (False, True):
            with self.subTest(failure=failure), self.fixture(behavior=
                    '  raise RuntimeError("resolver_fixture_failure")\n' if failure else '') as f:
                external = Path(f.baseline['UTILPATH'])
                external.mkdir(parents=True, exist_ok=True)
                (external / 'physics.root').write_bytes(b'never read')
                test = PrivateBoundary(methodName='runTest')
                test.exercise_exact_preflight_sigma0('absent',
                    probe_check=lambda *args: real_probe(f.repo, f.tree, f.env, f.baseline, f.outdir, f.package),
                    resolver_failure=failure)

    def test_farm_no_mode_alias_and_actual_csv_sink(self):
        import csv
        with self.fixture() as f:
            result = self.run_private(f)
            base = next(r for r in result['baseline_callers'] if r['caller'].endswith('findEffectiveCharge.py'))
            private = next(r for r in result['all_callers'] if r['caller'].endswith('findEffectiveCharge.py'))
            normalized = Path(base['paths']['OUTPATH'].replace('None', 'KaonLT'))
            self.assertIn('hallc_replay_lt/UTIL_KAONLT/OUTPUT/Analysis/None', base['paths']['OUTPATH'])
            self.assertNotEqual(normalized, Path(f.baseline['OUTPATH']))
            self.assertEqual(normalized.resolve(), Path(f.baseline['OUTPATH']).resolve())
            self.assertEqual(Path(private['paths']['OUTPATH']).parent, normalized.parent)
            sink = owner.effective_charge_outpath(f.tree, private['paths']) / 'table_kaon_Q4p4W2p74.csv'
            self.assertEqual(sink.parent.resolve(), f.outdir)
            # Execute the existing standard-library CSV sink, without importing
            # utility.py's analysis dependencies or calculating any physics.
            tree = ast.parse((ROOT / 'src/utility/utility.py').read_text())
            function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'data_to_csv')
            namespace = {'os': os, 'csv': csv}
            exec(compile(ast.Module(body=[function], type_ignores=[]), '<csv-sink>', 'exec'), namespace)
            namespace['data_to_csv'](str(sink), 'Total Efficiency', 0.75, '1234')
            self.assertTrue(sink.is_file())
            self.assertFalse((Path(f.baseline['OUTPATH']) / sink.name).exists())
            sink.unlink()

    def test_farm_wrong_baseline_and_shared_private_writer_fail(self):
        real_probe = owner.probe_private_callers
        for failure in ('baseline', 'private'):
            with self.subTest(failure=failure), self.fixture() as f:
                link = (Path(f.baseline['UTILPATH']) / 'OUTPUT' if failure == 'baseline' else f.tree / 'OUTPUT')
                target = os.readlink(link)
                wrong = f.root / 'wrong-output'
                (wrong / 'Analysis/KaonLT').mkdir(parents=True)
                link.unlink()
                link.symlink_to(wrong if failure == 'baseline' else f.repo / 'OUTPUT', target_is_directory=True)
                try:
                    reason = 'baseline_root_outpath_mismatch' if failure == 'baseline' else 'default_root_private_outpath_mismatch'
                    with self.assertRaisesRegex(ValueError, reason):
                        self.run_private(f)
                    if failure == 'private':
                        def probe(*args):
                            try:
                                return real_probe(f.repo, f.tree, f.env, f.baseline, f.outdir, f.package)
                            except ValueError as error:
                                raise ValueError('resolver_fixture_failure:' + str(error)) from error
                        PrivateBoundary(methodName='runTest').exercise_exact_preflight_sigma0(
                            'absent', probe_check=probe, resolver_failure=True, via_cli=True)
                finally:
                    link.unlink()
                    link.symlink_to(target, target_is_directory=True)

    def test_private_no_mode_can_already_point_private(self):
        with self.fixture(behavior='  if not modes and "runtime-python" in __file__: self.OUTPATH=str(Path(self.LTANAPATH)/"OUTPUT/Analysis/None")\n') as f:
            self.run_private(f)

    def test_private_backing_tree_alias_escape_fails(self):
        with self.fixture() as f:
            output = f.tree.parent / 'matched-output'
            saved = output.with_name('saved-output')
            output.rename(saved)
            output.symlink_to(f.repo / 'OUTPUT', target_is_directory=True)
            try:
                with self.assertRaisesRegex(ValueError, 'default_root_private_outpath_mismatch'):
                    self.run_private(f)
            finally:
                output.unlink()
                saved.rename(output)


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
    def test_old_inventory_reproduces_all_cycle_shapes_without_traversal(self):
        # Execute the actual audited pre-repair functions, not a conjectured mock.
        source = subprocess.check_output(['git', 'show',
            '842dab7a957a9823a51dcb484a765b4d11c6d9c1:' + owner.OWNER], cwd=ROOT).decode()
        functions = [n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef)
                     and n.name in ('inventory', 'link_inventory')]
        namespace = dict(vars(owner))
        exec(compile(ast.Module(body=functions, type_ignores=[]), '<old-owner>', 'exec'), namespace)
        for shape in ('self', 'two_link', 'relative_ancestor'):
            with self.subTest(shape=shape), tempfile.TemporaryDirectory() as d:
                root = Path(d).resolve()
                trial = root / 'cache_transfer/kaon/Q3p0W3p14/Trial_44'
                trial.mkdir(parents=True)
                link = trial / 'Trial_49'
                literal = {'self': 'Trial_49', 'two_link': 'peer',
                           'relative_ancestor': '../Trial_44/Trial_49'}[shape]
                link.symlink_to(literal, target_is_directory=True)
                if shape == 'two_link':
                    (trial / 'peer').symlink_to('Trial_49', target_is_directory=True)
                # A resolvable ancestor directory link is also never descended.
                (trial / 'ancestor').symlink_to('..', target_is_directory=True)
                walked = [Path(parent).relative_to(root).as_posix()
                          for parent, _, _ in os.walk(root, followlinks=False)]
                self.assertEqual(len(walked), 5)
                self.assertFalse(any('/Trial_49/' in p or '/ancestor/' in p for p in walked))
                with self.assertRaises((RuntimeError, OSError)):
                    namespace['inventory'](root)
                result = owner.inventory(root)
                key = link.relative_to(root).as_posix()
                self.assertEqual(result[key], {'link': literal, 'cyclic_target': True})
                self.assertEqual(owner.inventory(root), result)
                self.assertEqual(len(result), 6 + (shape == 'two_link'))

    def test_loop_exception_forms_are_metadata_without_retry(self):
        for error in (RuntimeError("Symlink loop from 'fixture'"), OSError(errno.ELOOP, 'cycle')):
            with self.subTest(error=error), tempfile.TemporaryDirectory() as d:
                link = Path(d) / 'link'
                link.symlink_to('link')
                with patch.object(Path, 'resolve', side_effect=error) as resolve, patch.object(
                        Path, 'is_file', side_effect=AssertionError('cyclic referent must not be read')):
                    self.assertEqual(owner.link_inventory(link), {'link': 'link', 'cyclic_target': True})
                resolve.assert_called_once_with(strict=True)

    def test_other_resolution_and_readlink_errors_propagate(self):
        errors = (PermissionError(errno.EACCES, 'denied'), OSError(errno.EIO, 'IO'),
                  OSError(errno.ENOTDIR, 'not directory'), RuntimeError('unrelated failure'))
        with tempfile.TemporaryDirectory() as d:
            link = Path(d) / 'link'
            link.symlink_to('link')
            for error in errors:
                with self.subTest(error=error), patch.object(Path, 'resolve', side_effect=error):
                    with self.assertRaises(type(error)) as caught:
                        owner.link_inventory(link)
                    self.assertIs(caught.exception, error)
            with patch.object(os, 'readlink', side_effect=PermissionError('readlink denied')):
                with self.assertRaisesRegex(PermissionError, 'readlink denied'):
                    owner.link_inventory(link)
            link.unlink()
            target = Path(d) / 'regular'
            target.write_bytes(b'live')
            link.symlink_to(target)
            with patch.object(owner, 'identity', side_effect=PermissionError('hash denied')):
                with self.assertRaisesRegex(PermissionError, 'hash denied'):
                    owner.link_inventory(link)

    def test_all_cycle_inventory_callers_reach_exactly_one_child_for_both_forms(self):
        for form in ('native', 'runtime', 'oserror'):
            with self.subTest(form=form):
                self.exercise_exact_preflight_sigma0('absent', cyclic=True, loop_form=form)

    def test_cycle_and_live_link_changes_fail_after_child_failure(self):
        for location in ('shared', 'CUTPATH', 'PARAMPATH'):
            for change in ('retarget', 'remove', 'replace', 'resolve', 'bytes', 'live_bytes', 'live_cycle'):
                with self.subTest(location=location, change=change):
                    self.exercise_exact_preflight_sigma0('absent', (location, change), cyclic=True)

    def test_required_cyclic_inputs_stop_before_child(self):
        for role in ('sigma0_SIMC', 'neutron', 'delta', 'sidis'):
            with self.subTest(role=role):
                self.exercise_exact_preflight_sigma0('cycle', required_role=role)

    def test_preflight_traceback_and_independent_preservation_failures(self):
        self.exercise_exact_preflight_sigma0('absent', preflight_failure=True,
                                            preservation_failure=True)

    def test_completed_child_failed_artifact_gate_never_collects(self):
        self.exercise_exact_preflight_sigma0('absent', cyclic=True, artifact_failure=True)

    def test_completed_child_changed_cycle_never_collects(self):
        for location in ('shared', 'CUTPATH', 'PARAMPATH'):
            with self.subTest(location=location):
                self.exercise_exact_preflight_sigma0('absent', (location, 'resolve'),
                                                    cyclic=True, completed_child=True)

    def test_shared_dangling_link_metadata_and_preservation(self):
        for change in ('retarget', 'remove', 'replace', 'resolve', 'bytes'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as d:
                root = Path(d).resolve()
                shared = root / 'KaonLT'
                shared.mkdir()
                target = root / 'OUTPUT/Analysis/HeeP'
                target.parent.mkdir(parents=True)
                link = shared / 'unrelated'
                link.symlink_to(target, target_is_directory=True)
                regular = shared / 'regular.bin'
                regular.write_bytes(b'original')
                before = owner.inventory(shared)
                self.assertEqual(before['unrelated'], {'link': str(target), 'missing_target': True})
                self.assertEqual(owner.inventory(shared), before)
                self.assertFalse(target.exists())
                if change in ('retarget', 'remove', 'replace'):
                    link.unlink()
                if change == 'retarget':
                    link.symlink_to(root / 'other-missing')
                elif change == 'replace':
                    link.write_bytes(b'replacement')
                elif change == 'resolve':
                    target.mkdir(parents=True)
                elif change == 'bytes':
                    regular.write_bytes(b'mutated')
                self.assertNotEqual(owner.inventory(shared), before)

    def test_external_dangling_metadata_is_preserved(self):
        for role in ('CUTPATH', 'PARAMPATH'):
            for change in ('retarget', 'remove', 'replace', 'resolve'):
                with self.subTest(role=role, change=change), tempfile.TemporaryDirectory() as d:
                    root = Path(d).resolve()
                    worktree, private = root / 'analysis', root / 'private'
                    paths = external_fixture(worktree, private, root / 'volatile')
                    directory = Path(paths[role])
                    directory.mkdir(exist_ok=True)
                    target = root / 'OUTPUT/Analysis/HeeP'
                    target.parent.mkdir(parents=True)
                    link = directory / 'unrelated'
                    link.symlink_to(target, target_is_directory=True)
                    with patch.object(owner, 'skim_paths', return_value={'effective_charge': str(root / 'skim')}):
                        records = owner.external_input_inventory(worktree, paths, private, {})
                    self.assertEqual(records[str(link)], {'role': role, 'unresolved_link': {
                        'link': str(target), 'missing_target': True}})
                    owner.verify_external_inputs(records)
                    self.assertFalse(target.exists())
                    if change in ('retarget', 'remove', 'replace'):
                        link.unlink()
                    if change == 'retarget':
                        link.symlink_to(root / 'other-missing')
                    elif change == 'replace':
                        link.write_bytes(b'replacement')
                    elif change == 'resolve':
                        target.mkdir(parents=True)
                    with self.assertRaisesRegex(ValueError, 'external_input_changed'):
                        owner.verify_external_inputs(records)

    def test_required_dangling_inputs_still_fail_closed(self):
        for role in ('sigma0_SIMC', 'pion_SIMC_neutron'):
            with self.subTest(role=role), tempfile.TemporaryDirectory() as d:
                root = Path(d).resolve()
                worktree, private = root / 'analysis', root / 'private'
                paths = external_fixture(worktree, private, root / 'volatile')
                env = {}
                if role == 'sigma0_SIMC':
                    link = root / 'sigma.root'
                    env['LT_BG_SIGMA0_CENTER_LOW_ROOT'] = str(link)
                else:
                    link = root / 'volatile/OUTPUT/Analysis/KaonLT/background_samples/neutron/Prod_Coin_Q4p4W2p74center_lowe.root'
                    link.unlink()
                link.symlink_to(root / 'missing.root')
                with self.assertRaises(FileNotFoundError):
                    owner.identity(link)
                with patch.object(owner, 'skim_paths', return_value={'effective_charge': str(root / 'skim')}):
                    with self.assertRaisesRegex(ValueError, 'required_external_input_missing:' + role):
                        owner.external_input_inventory(worktree, paths, private, env)

    def test_intact_link_keeps_referent_identity(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d).resolve()
            target = root / 'regular.bin'
            target.write_bytes(b'original')
            link = root / 'link'
            link.symlink_to(target)
            before = owner.inventory(root)
            self.assertEqual(before['link'], {'link': str(target),
                'resolved_path': str(target), 'referent': owner.identity(link)})
            target.write_bytes(b'mutated')
            self.assertNotEqual(owner.inventory(root), before)

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

    def test_exact_preflight_dangling_sigma0_stops_before_child(self):
        self.exercise_exact_preflight_sigma0('dangling')

    def test_child_dangling_link_mutations_fail_preservation(self):
        for location in ('shared', 'CUTPATH', 'PARAMPATH'):
            for change in ('retarget', 'remove', 'replace', 'resolve', 'bytes'):
                with self.subTest(location=location, change=change):
                    self.exercise_exact_preflight_sigma0('absent', (location, change))

    def test_cheap_full_acquisition_failure_does_not_read_large_shared_tree(self):
        self.exercise_exact_preflight_sigma0('absent', cheap_bulk_failure=True)

    def test_missing_full_acquisition_caller_stops_before_audit_and_bulk_reads(self):
        self.exercise_exact_preflight_sigma0('absent', cheap_bulk_failure='caller')

    def test_success_has_distinct_precollection_and_final_byte_checks_without_duplicate_final_scan(self):
        self.exercise_exact_preflight_sigma0('absent', complete_owner=True)

    def test_cleanup_mutation_still_blocks_delivery_after_collection(self):
        self.exercise_exact_preflight_sigma0('absent', complete_owner=True, cleanup_mutation=True)

    def test_status_io_failure_cannot_suppress_final_preservation(self):
        self.exercise_exact_preflight_sigma0('absent', complete_owner=True, final_status_failure=True)

    def test_full_acquisition_accepts_authorized_cli_root_aliases(self):
        self.exercise_exact_preflight_sigma0('absent', complete_owner=True, root_alias=True)

    def exercise_exact_preflight_sigma0(self, sigma_case, mutation=None, *, cyclic=False,
                                       loop_form='native', required_role='sigma0_SIMC',
                                       preflight_failure=False, preservation_failure=False,
                                       artifact_failure=False, completed_child=False,
                                       probe_check=None, resolver_failure=False,
                                       cheap_bulk_failure=False, complete_owner=False, cleanup_mutation=False,
                                       final_status_failure=False, root_alias=False, via_cli=False):
        completed_child = completed_child or complete_owner
        with tempfile.TemporaryDirectory() as d:
            parent = Path(d).resolve()
            repo, worktree, globus = parent / 'repo', parent / 'analysis', parent / 'globus'
            for p in (repo, worktree, globus, parent / 'runtime-python'):
                p.mkdir()
            entry_repo, entry_shared, entry_output = repo, None, globus / 'new.zip'
            if root_alias:
                alias = parent / 'group'
                alias.symlink_to(parent, target_is_directory=True)
                entry_repo, entry_output = alias / 'repo', alias / 'globus/new.zip'
            # The cheap source/caller gate now precedes all scientific I/O.
            for tree in (repo, worktree):
                names = {caller for _, caller, _, _, _ in owner.resolver.CALLERS}
                names.update(('src/main.py', 'src/setup/calcTotalEffectiveCharge.py',
                              'src/setup/set_sig_fortran.py'))
                for name in names - {'.'}:
                    p = tree / name
                    p.parent.mkdir(parents=True, exist_ok=True)
                    p.write_text('# source-only caller fixture\n')
            (parent / 'runtime-python/ltsep').mkdir()
            (worktree / 'testing').mkdir()
            (worktree / owner.OWNER).write_text('synthetic owner')
            for name in ('src/cuts/shift_prep.py', 'src/setup/shift_MM.py'):
                p = worktree / name
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text('# frozen writer fixture')
            shared = parent / 'volatile/OUTPUT/Analysis/KaonLT'
            entry_shared = alias / 'volatile/OUTPUT/Analysis/KaonLT' if root_alias else shared
            (shared / 'background_samples').mkdir(parents=True)
            for setting in owner.SETTINGS:
                for role in ('Data', 'Dummy'):
                    (shared / f'{setting}_kaon_Analysed_{role}_Q4p4W2p74_lowe.root').write_bytes((setting + role).encode())
                for suffix in ('.root', '.hist'):
                    (shared / f'Prod_Coin_Q4p4W2p74{setting.lower()}_lowe{suffix}').write_bytes(b'simc')
            # Exercise execute_gate -> external_input_inventory directly;
            # only the eventual analysis child is simulated and fails.
            external_paths = external_fixture(worktree, shared, parent / 'volatile')
            # Every exact preflight exercises both broad inventory callers.
            target = parent / 'volatile/OUTPUT/Analysis/HeeP'
            links = {}
            peers, live_targets = {}, {}
            for location, directory in dict(shared=str(shared), **external_paths).items():
                if location not in ('shared', 'CUTPATH', 'PARAMPATH'):
                    continue
                directory = Path(directory)
                directory.mkdir(exist_ok=True)
                (directory / 'regular.bin').write_bytes(b'original')
                links[location] = directory / 'unrelated'
                if cyclic:
                    peers[location] = parent / (location + '-peer')
                    peers[location].symlink_to(links[location])
                    links[location].symlink_to(peers[location], target_is_directory=True)
                else:
                    links[location].symlink_to(target, target_is_directory=True)
                live_targets[location] = parent / (location + '-live.bin')
                live_targets[location].write_bytes(b'live original')
                (directory / 'live-link').symlink_to(live_targets[location])
            sigma_env = {}
            if sigma_case == 'cycle':
                p = (parent / 'sigma-cycle.root' if required_role == 'sigma0_SIMC' else
                     shared / 'background_samples' / required_role / 'Prod_Coin_Q4p4W2p74center_lowe.root')
                if p.exists():
                    p.unlink()
                p.symlink_to(p.name)
                if required_role == 'sigma0_SIMC':
                    sigma_env['LT_BG_SIGMA0_CENTER_LOW_ROOT'] = str(p)
                with self.assertRaises((OSError, RuntimeError)):
                    owner.identity(p)
            elif sigma_case != 'absent':
                for s in owner.SETTINGS:
                    p = parent / (s + '-sigma.root')
                    if sigma_case == 'valid':
                        p.write_bytes((s + ' Sigma0').encode())
                    elif sigma_case == 'dangling':
                        p.symlink_to(parent / 'missing-sigma.root')
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
                    if cleanup_mutation and args[-1] == 'analysis':
                        (shared / 'regular.bin').write_bytes(b'cleanup mutated original')
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
                # Use the real atomic JSON status writer, including on failure.
                real_resolve = Path.resolve
                def resolve(path, *args, **kwargs):
                    if path in links.values() and loop_form != 'native':
                        if loop_form == 'runtime':
                            raise RuntimeError('Symlink loop from ' + repr(str(path)))
                        raise OSError(errno.ELOOP, 'cycle', str(path))
                    return real_resolve(path, *args, **kwargs)
                if cyclic:
                    stack.enter_context(patch.object(Path, 'resolve', new=resolve))
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
                scans = patched(owner, 'inventory', wraps=owner.inventory)
                bulk = shared / 'large-unrelated'
                if cheap_bulk_failure:
                    bulk.mkdir()
                    with (bulk / 'large.root').open('wb') as stream:
                        stream.truncate(8 * 1024**3)
                    for i in range(128):
                        (bulk / str(i)).mkdir()
                    scans.side_effect = AssertionError('cheap gate must not scan shared or external trees')
                    hashes = patched(owner, 'identity', side_effect=AssertionError('cheap gate must not hash science'))
                    if cheap_bulk_failure == 'caller':
                        (repo / 'src/setup/findEffectiveCharge.py').unlink()
                captured = {}
                real_inventory = owner.external_input_inventory
                def capture_inventory(*args, **kwargs):
                    records = real_inventory(*args, **kwargs)
                    if not kwargs.get('check_only'):
                        captured.update(records)
                    return records
                inventory_call = patched(owner, 'external_input_inventory', side_effect=capture_inventory)
                if preflight_failure:
                    def fail_source_caller(*args, **kwargs):
                        raise PermissionError('synthetic preflight denied')
                    inventory_call.side_effect = fail_source_caller
                patched(owner, 'skim_paths', return_value={'launcher': str(shared), 'effective_charge': str(shared)})
                if probe_check is None:
                    probe = patched(owner, 'probe_private_callers', return_value={})
                    if cheap_bulk_failure:
                        probe.side_effect = ValueError('synthetic_cheap_path_failure')
                else:
                    patched(owner, 'probe_private_callers', side_effect=probe_check)
                patched(owner.event_owner, 'committed_identities', return_value={})
                def failed_child(tree, log, env):
                    Path(log).write_text('All available low-epsilon settings completed in debug mode.\n'
                                         'Full high-epsilon processing is intentionally skipped in -d debug mode.'
                                         if artifact_failure or completed_child else
                                         'synthetic child failure; no analysis executed')
                    if cyclic:
                        for location in ('CUTPATH', 'PARAMPATH'):
                            record = captured[str(links[location])]
                            self.assertEqual(record, {'role': location, 'unresolved_link': {
                                'link': str(peers[location]), 'cyclic_target': True}})
                        owner.verify_external_inputs(captured)
                    if mutation:
                        location, change = mutation
                        link = links[location]
                        if change in ('retarget', 'remove', 'replace'):
                            link.unlink()
                        if change == 'retarget':
                            link.symlink_to(parent / 'other-missing')
                        elif change == 'replace':
                            link.write_bytes(b'replacement')
                        elif change == 'resolve':
                            if cyclic:
                                peers[location].unlink()
                                peers[location].write_bytes(b'repaired referent')
                            else:
                                target.mkdir()
                        elif change == 'live_bytes':
                            live_targets[location].write_bytes(b'mutated referent')
                        elif change == 'live_cycle':
                            live_targets[location].unlink()
                            live_targets[location].symlink_to(links[location])
                        elif change == 'bytes':
                            regular = (shared / 'Center_kaon_Analysed_Data_Q4p4W2p74_lowe.root'
                                       if location == 'shared' else Path(external_paths[location]) / 'regular.bin')
                            regular.write_bytes(b'mutated')
                    return 0 if artifact_failure or completed_child else 1
                child = patched(owner.ab.debug, 'run_analysis', side_effect=failed_child)
                if artifact_failure:
                    patched(owner.ab, 'committed_identities', return_value={})
                    patched(owner, 'verify_joint', side_effect=ValueError('synthetic joint gate failure'))
                elif completed_child:
                    patched(owner.ab, 'committed_identities', return_value={})
                    patched(owner, 'verify_joint', return_value=([], {}, {
                        s: {'input_identity': {}} for s in owner.SETTINGS}))
                checkout = patched(owner.isolation, 'verify_preservation', return_value={})
                installed = patched(owner.isolation, 'verify_ltsep_preservation', return_value={})
                delivery = patched(owner.isolation, 'deliver_evidence')
                if final_status_failure:
                    update = owner.ab.GateStatus.update
                    def fail_status(instance, stage=None, **fields):
                        if stage == 'final_preservation':
                            raise ValueError('synthetic_status_io_failure')
                        return update(instance, stage, **fields)
                    patched(owner.ab.GateStatus, 'update', new=fail_status)
                if complete_owner:
                    def collected(**kwargs):
                        Path(kwargs['output']).write_bytes(b'synthetic ZIP; archive validator tested separately')
                        return {'returncode': 0}
                    collect.side_effect = collected
                    patched(owner, 'verify_archive', return_value={})
                if preservation_failure:
                    checkout.side_effect = ValueError('checkout preservation failure')
                    installed.side_effect = ValueError('installed preservation failure')
                invalid = sigma_case in ('invalid', 'dangling', 'cycle') or preflight_failure or resolver_failure or cheap_bulk_failure
                role = 'sigma0_SIMC' if required_role == 'sigma0_SIMC' else 'pion_SIMC_' + required_role
                reason = ('missing_or_nonregular:.*canonical.*lexical' if cheap_bulk_failure == 'caller' else
                          'synthetic_status_io_failure' if final_status_failure else
                          'synthetic_cheap_path_failure' if cheap_bulk_failure else
                          'resolver_fixture_failure' if resolver_failure else 'preservation_failed' if mutation or cleanup_mutation else
                          'required_external_input_missing:' + role if invalid else
                          'synthetic joint gate failure' if artifact_failure else 'analysis_failed')
                expected = ValueError
                if sigma_case == 'cycle':
                    # is_file may propagate ELOOP on some Python versions.
                    expected = (ValueError, OSError, RuntimeError)
                    reason = 'required_external_input_missing:|[Ss]ymlink|[Ss]ymbolic'
                if preservation_failure:
                    reason = 'preservation_failed'
                if complete_owner and not cleanup_mutation and not final_status_failure:
                    owner.execute_gate(entry_repo, entry_shared, entry_output, head)
                elif via_cli:
                    import io
                    from contextlib import redirect_stderr
                    errors = io.StringIO()
                    with redirect_stderr(errors):
                        result = owner.main(['--repo', str(entry_repo), '--outdir', str(entry_shared),
                                             '--output', str(entry_output), '--source-commit', head])
                    self.assertEqual(result, 1)
                    self.assertIn('resolver_fixture_failure:default_root_private_outpath_mismatch', errors.getvalue())
                    self.assertIn('Owner gate status:', errors.getvalue())
                else:
                    with self.assertRaisesRegex(expected, reason):
                        owner.execute_gate(entry_repo, entry_shared, entry_output, head)
                shared_scans = [call for call in scans.call_args_list if Path(call.args[0]) == shared]
                self.assertEqual(len(shared_scans), 0 if cheap_bulk_failure or resolver_failure or invalid else
                                 3 if complete_owner else 2)
                if cheap_bulk_failure:
                    hashes.assert_not_called()
                    scans.assert_not_called()
                    shutil.rmtree(bulk)
                if resolver_failure or cheap_bulk_failure:
                    inventory_call.assert_not_called()
                else:
                    self.assertEqual(inventory_call.call_count, 1 if invalid else 2)
                    self.assertEqual(inventory_call.call_args.args[3], sigma_env)
                if invalid:
                    child.assert_not_called()
                else:
                    child.assert_called_once()
                    child_env = child.call_args.args[-1]
                    self.assertEqual(child_env[owner.ab.ENABLE_ENV], child_env[cap.ENABLE_ENV])
                    self.assertEqual({k: child_env[k] for k in sigma_env}, sigma_env)
                if complete_owner:
                    collect.assert_called_once()
                else:
                    collect.assert_not_called()
                    self.assertFalse((globus / 'new.zip').exists())
                if not complete_owner or cleanup_mutation or final_status_failure:
                    delivery.assert_not_called()
                else:
                    delivery.assert_called_once()
                status = json.loads(owner.isolation.accepted.gate_status_path(shared, globus / 'new.zip').read_text())
                if complete_owner and not cleanup_mutation and not final_status_failure:
                    self.assertEqual(status['status'], 'success')
                    if root_alias:
                        self.assertNotEqual(status['requested_paths'], status['canonical_paths'])
                        self.assertEqual(status['canonical_paths']['repo'], str(repo))
                    for name in ('source_preflight', 'overlay_path_validation', 'shared_output_snapshot',
                                 'private_ROOT_staging', 'external_input_verification', 'analysis',
                                 'post_child_preservation', 'collection', 'final_preservation'):
                        timing = status['stage_timings'][name]
                        self.assertTrue(timing['completed'])
                        self.assertGreaterEqual(timing['elapsed_seconds'], 0)
                    self.assertTrue(status['shared_inputs_outputs_preserved'])
                    return
                self.assertEqual(status['status'], 'failed')
                self.assertFalse(status['evidence_admissible'])
                if cheap_bulk_failure:
                    self.assertNotIn('shared_output_snapshot', status['stage_timings'])
                    self.assertFalse(status['stage_timings'][
                        'source_preflight' if cheap_bulk_failure == 'caller' else 'overlay_path_validation']['completed'])
                if cleanup_mutation or final_status_failure:
                    self.assertIn('synthetic_status_io_failure' if final_status_failure else 'shared_outputs_changed',
                                  status['preservation_failure_reason'])
                else:
                    self.assertTrue(status['failure_reason'])
                    self.assertIn('execute_gate', status['failure_traceback'])
                    self.assertLessEqual(len(status['failure_traceback']), 8192)
                if resolver_failure:
                    self.assertEqual(status['failed_stage'], 'overlay_path_validation')
                    self.assertFalse(status['analysis_started'])
                    self.assertIn('probe_private_callers', status['failure_traceback'])
                if preflight_failure:
                    self.assertEqual(status['failed_stage'], 'overlay_path_validation')
                    self.assertFalse(status['analysis_started'])
                    self.assertIn('fail_source_caller', status['failure_traceback'])
                    self.assertEqual(status['failure_reason'], 'synthetic preflight denied')
                    self.assertIn('ordinary_checkout:', status['preservation_failure_reason'])
                    self.assertIn('installed_ltsep:', status['preservation_failure_reason'])
                self.assertGreaterEqual(checkout.call_count, 1)
                if cheap_bulk_failure != 'caller':
                    self.assertGreaterEqual(installed.call_count, 1)
            if cheap_bulk_failure != 'caller':
                self.assertIn('analysis', cleaned)
            if not mutation and not cleanup_mutation:
                excluded = {p.name for p in owner.ab.destinations(shared, globus / 'new.zip').values()}
                self.assertEqual({k: v for k, v in owner.inventory(shared).items() if k not in excluded}, original)
                self.assertFalse(target.exists())
            self.assertEqual((shared / 'new.log').is_file(), not invalid)

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

    def test_effective_charge_change_is_path_only_and_old_sink_is_rejected(self):
        relative = 'src/setup/findEffectiveCharge.py'
        old = subprocess.check_output(['git', 'show', owner.BASE + ':' + relative], cwd=ROOT, text=True)
        new = (ROOT / relative).read_text()
        before, after = ast.parse(old), ast.parse(new)
        original_path = next(n for n in before.body if isinstance(n, ast.If) and
                             ast.unparse(n.test) == "'None' in OUTPATH")
        before.body.remove(original_path)
        new_path = [n for n in after.body if isinstance(n, ast.Assign) and
                    any(isinstance(t, ast.Name) and t.id == 'OUTPATH' for t in n.targets)][-1]
        after.body.remove(new_path)
        self.assertEqual(ast.dump(before), ast.dump(after))
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / relative
            target.parent.mkdir(parents=True)
            target.write_text(old)
            with self.assertRaisesRegex(ValueError, 'effective_charge_writer_path_changed'):
                owner.effective_charge_outpath(directory, {'LTANAPATH': directory, 'ANATYPE': 'Kaon'})

    def test_no_child_before_path_gates_and_no_retry(self):
        import inspect
        source = inspect.getsource(owner.execute_gate)
        self.assertLess(source.index('resolver_preflight('), source.index('run_owned_child('))
        self.assertLess(source.index('stage_input_files('), source.index('run_owned_child('))
        self.assertLess(source.index('verify_joint('), source.index('module.collect_validation_bundle('))
        self.assertEqual(source.count('run_owned_child('), 1)

    def test_smoke_cli_rejects_capture_and_offline_arguments(self):
        import contextlib, io
        for option in ('--outdir', '--output', '--output-dir', '--offline-bundle'):
            with self.subTest(option=option), patch.object(owner, 'resolver_smoke') as smoke:
                with contextlib.redirect_stderr(io.StringIO()):
                    self.assertEqual(owner.main(['--resolver-smoke', '--source-commit', 'c'*40,
                                                  option, '/unused']), 1)
                smoke.assert_not_called()

    def test_shared_source_authority_rejects_unknown_commit_and_dirty_index_source(self):
        good = dict(head='c'*40, origin_test='c'*40, branch='test', porcelain='')
        for change in ({'head': 'd'*40}, {'origin_test': 'd'*40}, {'branch': 'other'},
                       {'porcelain': 'M  testing/source.py'}, {'porcelain': ' M src/main.py'}):
            with self.subTest(change=change), patch.object(owner.isolation, 'snapshot',
                    return_value=dict(good, **change)), patch.object(owner.ab, 'git', return_value=''):
                with self.assertRaises(ValueError):
                    owner.published_source(ROOT, 'c'*40)
        with patch.object(owner.isolation, 'snapshot', return_value=good), patch.object(
                owner.ab, 'git', return_value='src/unauthorized.py'):
            with self.assertRaisesRegex(ValueError, 'unexpected_committed_source'):
                owner.published_source(ROOT, 'c'*40)

    def test_smoke_ordinary_import_failure_preserves_and_never_creates_worktree(self):
        state = dict(head='c'*40, origin_test='c'*40, branch='test', porcelain='')
        with ExitStack() as stack:
            def patched(obj, name, **kwargs):
                return stack.enter_context(patch.object(obj, name, **kwargs))
            patched(owner.isolation.accepted, 'REPO', new=ROOT)
            patched(owner.isolation, 'snapshot', return_value=state)
            patched(owner, 'smoke_source_state', return_value=state)
            patched(owner.ab, 'git', return_value='')
            patched(owner.shutil, 'which', return_value=owner.sys.executable)
            patched(owner, 'audit_default_output_writer', return_value={})
            patched(owner.isolation, 'import_identity', side_effect=ValueError('ordinary_import_failed'))
            worktree = patched(owner.isolation, 'owned_worktree')
            preservation = patched(owner.isolation, 'verify_preservation')
            with self.assertRaisesRegex(ValueError, 'ordinary_import_failed'):
                owner.resolver_smoke(ROOT, 'c'*40)
            worktree.assert_not_called()
            preservation.assert_called_once_with(ROOT, state)

    def test_smoke_and_acquisition_use_same_preflight(self):
        import inspect
        for function in (owner.resolver_smoke, owner.execute_gate):
            source = inspect.getsource(function)
            self.assertIn('published_source(repo, source_commit)', source)
            self.assertIn('isolation.import_identity(repo)', source)
            self.assertIn('isolation.prepare_runtime_overlay(repo, worktree)', source)
            self.assertIn('resolver_preflight(repo, worktree, env, paths,', source)
        preflight = inspect.getsource(owner.resolver_preflight)
        for call in ('external_symlink_preflight(', 'skim_paths(', 'probe_paths(',
                     'prepare_private_output(', 'probe_private_callers('):
            self.assertIn(call, preflight)

    def test_wrong_mode_destination_blocks(self):
        with tempfile.TemporaryDirectory() as d:
            tree = Path(d).resolve()
            base = {'VOLATILEPATH': str(tree / 'external'), 'OUTPATH': str(tree / 'external'),
                    'ANATYPE': 'Kaon', 'LTANAPATH': str(tree)}
            record = {'paths': base, 'package_file': str(tree / 'runtime-python/ltsep/__init__.py')}
            with patch.object(owner, 'reached_resolver_callers', return_value=[]), patch.object(owner, 'probe_resolver', return_value=record):
                with self.assertRaisesRegex(ValueError, 'baseline_root_outpath_mismatch'):
                    owner.probe_private_callers(tree, tree, {}, base, tree / 'private', tree)

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
