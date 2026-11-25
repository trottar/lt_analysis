"""Deterministic admission/runner/adapter tests; no actual SIMC/ROOT execution.

Optional KSIGMA0_TEST_ARCHIVE points to the original external ZIP for all real
member admissions. Synthetic rows/stub executables never establish farm physics.
"""
import copy
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

REPO = Path(__file__).resolve().parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


p = module('ksigma_prepare', REPO / 'background_samples/prepare_ksigma0.py')
v = module('ksigma_validate', REPO / 'src/setup/validate_ksigma0_output.py')


def executable(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(0o700)


def row(x=0., y=0., weight=1.):
    data = dict.fromkeys(v.IMMUTABLE, 0.)
    data.update(ssxfp=x, ssyfp=y, Weight=weight)
    return data


def hist():
    return """H(e,e'K) producing a SIGMA0
doing_kaon= T doing_pion= F doing_delta= F doing_semi= F doing_phsp= F
doing_hydkaon= T doing_decay= T using_rad= T mc_smear= T
using_Eloss= T correct_Eloss= T which_kaon= 1
Ngen (request) = 200000
Ncontribute = 200000
normfac = 0.463860D+08
Mtar_struc = 938.272 MeV Mrec_struc = 0.119264E+04 MeV
Ebeam = 10590.6631 MeV
"""


class AdmissionTests(unittest.TestCase):
    def test_exact_reviewed_lookup_and_all_25_stems(self):
        data = p.load_manifest()
        self.assertEqual(len({s['canonical_stem'] for s in data['settings']}), 25)
        for s in data['settings']:
            with self.subTest(stem=s['canonical_stem']):
                kin = s['kinematic'].replace('_', '')
                self.assertEqual(s['canonical_stem'], 'Prod_Coin_' + kin + s['setting'] + '_' + s['epsilon'] + 'e')
                self.assertEqual(p.select(data, [s['canonical_stem']]), [s])
                self.assertEqual(s['requested_successes'], 200000)
        for kin in {s['kinematic'] for s in data['settings']}:
            self.assertEqual({(s['setting'], s['epsilon']) for s in data['settings'] if s['kinematic'] == kin},
                             {('left', 'low'), ('center', 'low'), ('left', 'high'),
                              ('center', 'high'), ('right', 'high')})

    def test_unknown_right_low_duplicate_empty_and_hostile_setting(self):
        data = p.load_manifest()
        stem = data['settings'][0]['canonical_stem']
        for stems in ([], [stem, stem], ['Prod_Coin_Q4p4W2p74right_lowe'], ['x;touch nope'], ['other']):
            with self.subTest(stems=stems), self.assertRaises(p.AdmissionError):
                p.select(data, stems)

    def test_pinned_manifest_change_rejected(self):
        with patch.object(p, 'file_hash', return_value='0' * 64), self.assertRaises(p.AdmissionError):
            p.load_manifest()

    def archive(self, directory, hostile=None, alter=False):
        data = {'flags': p.load_manifest()['flags'], 'common_sources': {}, 'common_tables': {},
                'settings': [{'input_member': 's/input.inp', 'settings_member': 's/settings.json',
                              'historical_hist_member': 's/out.hist'}]}
        inp = '\n'.join(k + ' = ' + str(x) for k, x in data['flags'].items()).encode()
        cfg = json.dumps({'program_settings': data['flags'],
                          'recon_hcana_settings': {'mode': 'kaon', 'hadron_mass': 0.493677}}).encode()
        items = {'s/input.inp': inp, 's/settings.json': cfg, 's/out.hist': hist().encode()}
        for key in ('input', 'settings', 'historical_hist'):
            data['settings'][0][key + '_sha256'] = p.digest(items[data['settings'][0][key + '_member']])
        if alter:
            items['s/input.inp'] += b'changed'
        path = directory / 'archive.zip'
        with zipfile.ZipFile(path, 'w') as z:
            for name, b in items.items():
                z.writestr(name, b)
            if hostile:
                info = zipfile.ZipInfo(hostile)
                if hostile == 'link':
                    info.external_attr = (stat.S_IFLNK | 0o777) << 16
                z.writestr(info, b'outside')
        data['sha256'] = p.file_hash(path)
        return path, data

    def test_bounded_archive_positive_and_outer_member_hash_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, data = self.archive(Path(tmp))
            self.assertEqual(len(p.archive_payloads(path, data)), 3)
            data['sha256'] = '0' * 64
            with self.assertRaisesRegex(p.AdmissionError, 'archive SHA256'):
                p.archive_payloads(path, data)
            path, data = self.archive(Path(tmp), alter=True)
            with self.assertRaisesRegex(p.AdmissionError, 'member hash'):
                p.archive_payloads(path, data)

    def test_archive_traversal_absolute_symlink_and_duplicate(self):
        for name in ('../escape', '/absolute', 'a\\b', 'C:/x', 's/input.inp', 'link'):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmp:
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', UserWarning)
                    path, data = self.archive(Path(tmp), hostile=name)
                with self.assertRaises(p.AdmissionError):
                    p.archive_payloads(path, data)

    def test_archive_size_bound(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, data = self.archive(Path(tmp))
            with patch.object(p, 'MAX_MEMBER', 8), self.assertRaises(p.AdmissionError):
                p.archive_payloads(path, data)

    def test_wrong_reaction_and_duplicate_input_flag(self):
        data = p.load_manifest()
        raw = '\n'.join(k + ' = ' + str(x) for k, x in data['flags'].items()).encode()
        cfg = json.dumps({'program_settings': data['flags'],
                          'recon_hcana_settings': {'mode': 'kaon', 'hadron_mass': 0.493677}}).encode()
        for bad in (raw.replace(b'which_kaon = 1', b'which_kaon = 0'), raw + b'\nngen = 1'):
            with self.assertRaises(p.AdmissionError):
                p.validate_input(bad, cfg, data)

    def test_missing_and_wrong_approved_base_manifest(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'base.json'
            path.write_text('{}')
            for sha in ('', '0' * 64):
                with self.assertRaises(p.AdmissionError):
                    p.validate_base(Path(tmp), path, sha)

    def test_strict_private_paths_alias_missing_intermediate_and_overlap(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            target = root / 'base'
            target.mkdir()
            link = root / 'link'
            link.symlink_to(target, target_is_directory=True)
            for path in (link, link / 'missing' / 'x'):
                with self.assertRaisesRegex(p.AdmissionError, 'symlink'):
                    p.physical(path, exists=False)
            for a, b in ((root, target), (target, root), (root, root)):
                with self.assertRaises(p.AdmissionError):
                    p.disjoint(a, b)
            p.disjoint(target, root / 'private')

    def test_inventory_rejects_link_and_hardlink(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'a').write_text('a')
            os.link(root / 'a', root / 'b')
            with self.assertRaises(p.AdmissionError):
                p.inventory(root)

    def test_no_opt_in_and_bindings_fail_before_writes(self):
        with self.assertRaises(SystemExit), patch.object(p, 'prepare') as prepare:
            p.main(['--setting', 'anything'])
        prepare.assert_not_called()
        with patch.dict(os.environ, {'LT_BG_SIGMA0_CENTER_HIGH_ROOT': '/external'}), \
                patch.object(p, 'prepare') as prepare, self.assertRaises(p.AdmissionError):
            p.main(['--sigma-only', '--setting', 'x', '--archive', 'x', '--base', 'x',
                    '--base-manifest', 'x', '--base-manifest-sha256', '0' * 64, '--attempt-parent', 'x'])
        prepare.assert_not_called()

    def test_real_external_archive_all_members_and_historical_reports(self):
        archive = os.environ.get('KSIGMA0_TEST_ARCHIVE')
        if not archive:
            self.skipTest('original external archive not configured')
        data = p.load_manifest()
        payloads = p.archive_payloads(Path(archive), data)
        self.assertEqual(len(payloads), 900)
        for s in data['settings']:
            prefix = s['settings_member'].rsplit('/', 1)[0] + '/'
            raw = payloads[prefix + 'recon_hcana.C']
            expected = data['common_sources']['recon_hcana.C']['archive_sha256']
            derived = p.reconstruction_macro(raw, expected)
            self.assertNotEqual(derived, raw)
            self.assertEqual(derived.replace(p.FILE_OWNED_DESTRUCTOR, p.ARCHIVED_DESTRUCTOR), raw)
            result = v.parse_hist(payloads[s['historical_hist_member']].decode())
            self.assertEqual(result['normfac'], s['historical_normfac'])


    @unittest.skipUnless(shutil.which('git'), 'Git required for checkout bytes test')
    def test_published_manifest_survives_git_checkout(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            name = p.MANIFEST.relative_to(REPO)
            manifest = root / name
            manifest.parent.mkdir()
            shutil.copyfile(p.MANIFEST, manifest)
            shutil.copyfile(REPO / '.gitattributes', root / '.gitattributes')
            subprocess.run(['git', 'init', '-q', str(root)], check=True)
            git = ['git', '-C', str(root), '-c', 'core.autocrlf=input']
            blob = subprocess.check_output(git + ['hash-object', '-w', str(name)]).strip()
            manifest.write_bytes(subprocess.check_output(git + ['cat-file', 'blob', blob.decode()]))
            self.assertEqual(p.file_hash(manifest), p.MANIFEST_SHA)
            with patch.object(p, 'MANIFEST', manifest):
                self.assertEqual(len(p.load_manifest()['settings']), 25)

    def test_reconstruction_derives_only_the_file_ownership_change(self):
        raw = b'prefix\n' + p.ARCHIVED_DESTRUCTOR + b'suffix\n'
        self.assertEqual(p.reconstruction_macro(raw, p.digest(raw)),
                         b'prefix\n' + p.FILE_OWNED_DESTRUCTOR + b'suffix\n')
        with self.assertRaisesRegex(p.AdmissionError, 'hash'):
            p.reconstruction_macro(raw, '0' * 64)
        duplicate = raw + p.ARCHIVED_DESTRUCTOR
        with self.assertRaisesRegex(p.AdmissionError, 'destructor'):
            p.reconstruction_macro(duplicate, p.digest(duplicate))


class AdapterTests(unittest.TestCase):
    def test_pyroot_loaded_identity_and_validator_drift_fail_before_events(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            core, module_path = root / 'libCore.so', root / '__init__.py'
            core.write_bytes(b'core fixture')
            module_path.write_bytes(b'ROOT fixture')
            class ROOT:
                __file__ = str(module_path)
                class gROOT:
                    @staticmethod
                    def GetVersionInt(): return 62408
            real_open = open
            def opener(name, *args, **kwargs):
                if name == '/proc/self/maps':
                    return io.StringIO('000-111 r-xp 0 00:00 1 ' + str(core) + '\n')
                return real_open(name, *args, **kwargs)
            with patch('builtins.open', side_effect=opener):
                identity = v.root_identity(ROOT)
            self.assertEqual(identity['core_sha256'], p.file_hash(core))
            self.assertEqual(identity['version_int'], 62408)
            with patch.dict(sys.modules, {'ROOT': ROOT}), \
                    patch.object(v, 'root_identity', return_value=identity), \
                    patch.object(v, 'open_tree') as events, self.assertRaisesRegex(ValueError, 'identity changed'):
                v.validate_attempt(root, {'root_runtime': {'pyroot': {'wrong': True}}})
            events.assert_not_called()

    def test_current_sidecar_fields_and_sigma_recoil(self):
        parsed = v.parse_hist(hist())
        self.assertEqual(parsed['Ncontribute'], 200000)
        self.assertEqual(parsed['normfac'], 46386000)
        for old, new in (('SIGMA0', 'LAMBDA'), ('which_kaon= 1', 'which_kaon= 0'),
                         ('200000', '100'), ('normfac = 0.463860D+08', 'normfac = -1'),
                         ('0.119264E+04', '0.111568E+04'), ('doing_decay= T', 'doing_decay= F')):
            with self.subTest(old=old), self.assertRaises(ValueError):
                v.parse_hist(hist().replace(old, new))
        with self.assertRaises(ValueError):
            v.parse_hist(hist() + '\nnormfac = 3')

    def test_tray_boundaries_all_kinematics_and_float_rounding(self):
        for s in p.load_manifest()['settings']:
            stem = s['canonical_stem']
            tight = s['kinematic'] in ('Q4p4_W2p74', 'Q3p0_W3p14', 'Q5p5_W3p02')
            x, y = (45., 30.) if tight else (55., 50.)
            self.assertTrue(v.tray(row(), stem))
            for r in (row(x), row(-x), row(y=y), row(y=-y), row(x + 1)):
                self.assertFalse(v.tray(r, stem))
        self.assertFalse(v.tray(row(45. - 1e-8), 'Prod_Coin_Q4p4W2p74center_highe'))

    def test_surviving_order_signed_weight_and_sumw2(self):
        raw = [row(weight=2.), row(50., weight=90.), row(weight=-3.), row(weight=0.)]
        after = [raw[0], raw[2], raw[3]]
        result = v.check_rows(raw, after, 'Prod_Coin_Q4p4W2p74center_highe')
        self.assertEqual(result['surviving_entries'], 3)
        self.assertEqual(result['sum_Weight'], -1)
        self.assertEqual(result['sumw2_Weight'], 13)
        for bad in (after[::-1], after[:-1], raw, [dict(after[0], Weight=4), *after[1:]],
                    [dict(after[0], hsdelta=1), *after[1:]]):
            with self.assertRaises(ValueError):
                v.check_rows(raw, bad, 'Prod_Coin_Q4p4W2p74center_highe')

    def test_nonfinite_empty_and_zero_weight_rejected(self):
        for raw in ([], [row(weight=math.nan)], [row(weight=0)]):
            with self.assertRaises(ValueError):
                v.check_rows(raw, raw, 'Prod_Coin_Q4p4W2p74center_highe')

    def test_identity_clone_stub_preserves_all_branches_and_order(self):
        class Tree:
            Weight = 0.
            def GetListOfBranches(self):
                return []
            def CloneTree(self, count):
                return self
            def Branch(self, name, buffer, leaf):
                self.buffer = buffer
                self.values = []
                self.args = (name, leaf)
                return True
            def GetEntries(self):
                return 3
            def GetEntry(self, i):
                self.Weight = [2., -3., 0.][i]
                return 1
            def Fill(self):
                self.values.append(self.buffer[0])
                return 1
            def Write(self):
                return 1
        class File:
            def IsZombie(self): return False
            def cd(self): pass
            def Close(self): pass
        class ROOT:
            class TFile:
                @staticmethod
                def Open(*args): return File()
        with tempfile.TemporaryDirectory() as tmp:
            tree = Tree()
            v.clone_identity(ROOT, tree, Path(tmp) / 'candidate.root')
            self.assertEqual(tree.values, [2., -3., 0.])
            self.assertEqual(tree.args, ('iter_weight', 'iter_weight/D'))
            (Path(tmp) / 'candidate.root').touch()
            with self.assertRaises(ValueError):
                v.clone_identity(ROOT, tree, Path(tmp) / 'candidate.root')

    def test_full_root_validator_stub_and_incompatible_branches_or_weight_copy(self):
        class Named:
            def __init__(self, name): self.name = name
            def GetName(self): return self.name
            def GetLenStatic(self): return 1
            def GetLeafCount(self): return None
            def GetTypeName(self): return 'Float_t'
        class Tree:
            def __init__(self, data): self.data = data
            def GetListOfBranches(self): return [Named(n) for n in self.data[0]]
            def GetLeaf(self, name): return Named(name) if name in self.data[0] else None
            def GetEntries(self): return len(self.data)
            def GetEntry(self, i): self.current = self.data[i]; return 1
            def __getattr__(self, name): return self.current[name]
            def InheritsFrom(self, name): return name == 'TTree'
        class File:
            def __init__(self, tree): self.tree = tree
            def IsZombie(self): return False
            def Get(self, name): return self.tree if name == 'h10' else None
            def Close(self): pass
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stem = 'Prod_Coin_Q4p4W2p74center_highe'
            (root / 'OUTPUTS').mkdir()
            (root / 'OUTPUTS' / (stem + '.hist')).write_text(hist())
            raw = [row(weight=2.), row(x=50.), row(weight=-3.)]
            after = []
            for r in (raw[0], raw[2]):
                r = dict(r, **dict.fromkeys(set(v.REQUIRED) - set(r), 0.))
                r.update(Q2=4.4, W=2.74, t=-0.5, missmass=1.19264)
                after.append(r)
            for failure in (None, 'missing-axis', 'wrong-units', 'lossy-copy'):
                data = copy.deepcopy(after)
                if failure == 'missing-axis':
                    del data[0]['phgcer_x_det']
                if failure == 'wrong-units':
                    data[0]['phipq'] = 180.
                files = {str(root / 'recon_hcana' / (stem + '.root')): Tree(raw),
                         str(root / 'recon_hcana' / (stem + '_recon.root')): Tree(data)}
                class ROOT:
                    class TFile:
                        @staticmethod
                        def Open(path, mode): return File(files[path])
                def clone(ROOT, source, target):
                    copied = [dict(r, iter_weight=r['Weight']) for r in source.data]
                    if failure == 'lossy-copy':
                        copied[0]['iter_weight'] = 1.
                    files[str(target)] = Tree(copied)
                parsed = v.parse_hist(hist())
                parsed['Ncontribute'] = 3  # Explicit tiny ROOT stub; not generator evidence.
                with self.subTest(failure=failure), patch.dict(sys.modules, {'ROOT': ROOT}), \
                        patch.object(v, 'parse_hist', return_value=parsed), \
                        patch.object(v, 'root_identity', return_value={'stub': True}), \
                        patch.object(v, 'clone_identity', side_effect=clone):
                    if failure:
                        with self.assertRaises(ValueError):
                            v.validate_attempt(root, {'stem': stem, 'root_runtime': {'pyroot': {'stub': True}}})
                    else:
                        result = v.validate_attempt(root, {'stem': stem, 'root_runtime': {'pyroot': {'stub': True}}})
                        self.assertTrue(result['iter_weight_identity'])
                        self.assertFalse(result['scientific_acceptance'])
                        self.assertEqual(result['sumw2_Weight'], 13.)


@unittest.skipUnless(os.name == 'posix', 'Linux ownership/runtime fixtures')
class PrivateRunnerTests(unittest.TestCase):
    def fixture(self, directory):
        base = directory / 'base'
        base.mkdir()
        sources = {'Makefile': 'all: physics_iterate.o\n', 'modules.f': 'module structureModule\n',
                   'simulate.inc': 'include interface\n', 'simc.f': 'native writer\n',
                   'results_write.f': 'native branches\n', 'shms/mc_shms_hut.f': 'detector\n',
                   'util/root_tree/Makefile': 'converter\n', 'util/root_tree/make_root_tree.f': 'native\n',
                   'dbase.f': "which_kaon.eq.1 sigma0\n" + '\n'.join("open(1,file='" + n + "')"
                       for n in ('weightc.dat', 'weightd.dat', 'saghai_proton.dat', 'saghai_sigma0.dat'))}
        for n, value in sources.items():
            dest = base / n
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_text(value)
        tools = directory / 'tools'
        for name in p.TOOLS:
            executable(tools / name, '#!/bin/bash\nprintf "fixture-version\\n"\n')
        tool_records = {n: {'path': str(tools / n), 'sha256': p.file_hash(tools / n),
                            'version_sha256': p.digest(b'fixture-version\n')} for n in p.TOOLS}
        env = {n: '' for n in p.ENV_KEYS}
        env['PATH'] = str(tools) + ':' + os.environ['PATH']
        lib = directory / 'fixture-library'
        lib.write_text('runtime')
        spec = {'schema': 'kaonlt_ksigma0_base/v1', 'commit': '1' * 40,
                'files': p.inventory(base), 'tracked': sorted(sources), 'tools': tool_records,
                'environment': env, 'runtime_files': {str(lib): p.file_hash(lib)},
                'active_reads': ['par.pl', 'weightc.dat', 'weightd.dat', 'saghai_proton.dat', 'saghai_sigma0.dat']}
        approval = directory / 'approved.json'
        approval.write_text(json.dumps(spec))
        return base, spec, approval, env

    def base_validate(self, base, spec, approval, env):
        real_command = p.command
        def command(args, cwd, env=None, input_text=None):
            if str(args[0]) == 'git':
                if args[-1] == 'HEAD':
                    return spec['commit'] + '\n'
                return '\0'.join(spec['tracked']) + '\0'
            return real_command(args, cwd, env, input_text)
        with patch.dict(os.environ, env), patch.object(p, 'command', side_effect=command), \
                patch.object(p, 'root_preflight', return_value={'stub': True}), \
                patch.object(p.subprocess, 'check_output', return_value=b''.join(
                    b'100644 blob ' + p.git_blob((base / n).read_bytes()).encode() + b'\t' + n.encode() + b'\0'
                    for n in spec['tracked'] if (base / n).exists())):
            return p.validate_base(base, approval, p.file_hash(approval))

    def test_approved_base_identity_and_missing_changed_dependencies(self):
        with tempfile.TemporaryDirectory() as tmp:
            base, spec, approval, env = self.fixture(Path(tmp))
            admitted = self.base_validate(base, spec, approval, env)[0]
            self.assertEqual({k: value for k, value in admitted.items() if not k.startswith('_')}, spec)
            (base / 'simulate.inc').write_text('changed')
            with self.assertRaisesRegex(p.AdmissionError, 'inventory/hash'):
                self.base_validate(base, spec, approval, env)
            (base / 'simulate.inc').unlink()
            with self.assertRaises(p.AdmissionError):
                self.base_validate(base, spec, approval, env)

    def prepared(self, directory):
        base, spec, approval, env = self.fixture(directory)
        data = copy.deepcopy(p.load_manifest())
        setting = data['settings'][18]
        prefix = setting['settings_member'].rsplit('/', 1)[0] + '/'
        payload = {}
        for n in set(p.STAGE) | set(data['common_tables']):
            raw = ('fixture ' + n).encode()
            if n == 'recon_hcana.C':
                raw += b'\n' + p.ARCHIVED_DESTRUCTOR
            payload[prefix + n] = raw
            if n in data['common_sources']:
                data['common_sources'][n]['archive_sha256'] = p.digest(raw)
            if n in data['common_tables']:
                data['common_tables'][n] = p.digest(raw)
        payload[setting['input_member']] = b'input fixture'
        setting['input_sha256'] = p.digest(b'input fixture')
        archive = directory / 'archive.zip'
        archive.write_bytes(b'fixture archive')
        data['sha256'] = p.file_hash(archive)
        core = directory / 'libCore.so'
        core.write_bytes(b'fixture core')
        spec['_root_runtime'] = {'pyroot': {'core': str(core), 'core_sha256': p.file_hash(core)}}
        approval.write_text(json.dumps(spec))
        parent = directory / 'attempts'
        parent.mkdir()
        root = p.prepare(base, parent, setting, payload, data, spec, p.file_hash(approval))
        return root, base, spec, data, setting, payload, parent

    def test_prepare_independent_files_exact_bytes_and_exclusive_collision(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, base, spec, data, setting, payload, parent = self.prepared(Path(tmp))
            self.assertEqual(p.inventory(base), spec['files'])
            for name, expected in json.loads((root / 'receipt.json').read_text())['loaded'].items():
                self.assertEqual(p.file_hash(root / name), expected)
                self.assertEqual((root / name).stat().st_nlink, 1)
            self.assertFalse(list((root / 'OUTPUTS').iterdir()))
            self.assertFalse((root / 'run_simc.py').exists())
            archived = (root / p.ARCHIVED_MACRO).read_bytes()
            derived = (root / p.RECON_MACRO).read_bytes()
            self.assertEqual(derived.replace(p.FILE_OWNED_DESTRUCTOR, p.ARCHIVED_DESTRUCTOR), archived)
            suffix = root.name.rsplit('-', 1)[1]
            class UUID:
                hex = suffix
            with patch.object(p.uuid, 'uuid4', return_value=UUID()), self.assertRaises(FileExistsError):
                p.prepare(base, parent, setting, payload, data, spec, 'approval')

    def fake_run(self, root, calls, fail=None, missing=None):
        receipt = json.loads((root / 'receipt.json').read_text())
        stem = receipt['stem']
        identity = receipt['root_runtime']['pyroot']
        def run(argv, cwd, env, input, text, stdout, stderr, check):
            calls.append(([str(x) for x in argv], str(cwd), input, dict(env)))
            label = Path(argv[0]).name
            if fail and label == fail:
                raise subprocess.CalledProcessError(3, argv)
            if str(argv[0]) == sys.executable:
                if argv[4] == 'validate' and fail == 'validator':
                    raise subprocess.CalledProcessError(3, argv)
                if argv[4] == 'probe':
                    report = {'different': True} if fail == 'private-probe' else identity
                else:
                    report = {'stub_only': True}
                stdout.write('KSIGMA0_VALIDATOR=' + json.dumps(report) + '\n')
            if label == 'ldd':
                core = identity['core']
                if fail == 'converter-core':
                    other_directory = root / 'tmp'
                    alternate = other_directory / 'libCore.so'
                    alternate.write_bytes(Path(core).read_bytes())
                    core = str(alternate)
                stdout.write(f'\tlibCore.so => {core} (0x123abc)\n')
            if label == 'make' and 'clean' not in argv:
                for n in ('simc', 'physics_iterate.o', 'util/root_tree/make_root_tree'):
                    (root / n).write_text('built')
                    (root / n).chmod(0o700)
            if label == 'simc':
                for n in ['worksim/' + stem + '.bin'] + ['OUTPUTS/' + stem + ext for ext in ('.hist', '.gen', '.geni')]:
                    if n != missing:
                        (root / n).write_text('fresh output')
            if label == 'make_root_tree':
                (root / 'worksim' / (stem + '.root')).write_text('raw root stub')
            if label == 'root':
                (root / 'recon_hcana' / (stem + '_recon.root')).write_text('recon root stub')
            # Filesystem timestamp granularity must not make a deterministic
            # stub predate its start. Keep the real stale-artifact guard strict.
            started = json.loads((root / 'receipt.json').read_text())['started_ns']
            for path in root.rglob('*'):
                if path.is_file() and (path.name in ('simc', 'physics_iterate.o', 'make_root_tree') or
                                       path.suffix in ('.root', '.bin', '.hist', '.gen', '.geni')):
                    os.utime(path, ns=(started + 1000000000, started + 1000000000))
            return subprocess.CompletedProcess(argv, 0)
        return run

    def execute_args(self, root, setting):
        parent = root.parent
        base = parent.parent / 'base'
        approval = parent.parent / 'approved.json'
        archive = parent.parent / 'archive.zip'
        return root, setting['canonical_stem'], base, approval, p.file_hash(approval), parent, archive

    def execute(self, root, data, run):
        # Stub builds only after the independent admission boundary.
        args = self.execute_args(root, data['settings'][18])
        spec = json.loads(args[3].read_text())
        with patch.object(p, 'load_manifest', return_value=data), \
                patch.object(p, 'validate_base', return_value=(spec, dict(spec['environment']))) as admit, \
                patch.object(p.subprocess, 'run', side_effect=run):
            p.execute(*args)
        admit.assert_called_once_with(args[2], args[3], args[4])

    def test_direct_runner_cannot_substitute_receipt_for_approval(self):
        for failure in ('wrong_digest', 'modified_receipt', 'extra_makefile', 'wrong_parent'):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                root, _, spec, data, setting, *_ = self.prepared(Path(tmp))
                args = list(self.execute_args(root, setting))
                if failure == 'wrong_digest':
                    args[4] = '0' * 64
                if failure == 'modified_receipt':
                    path = root / 'receipt.json'
                    receipt = json.loads(path.read_text())
                    receipt['tools']['make']['path'] = '/unapproved/make'
                    p.write_json(path, receipt)
                if failure == 'extra_makefile':
                    (root / 'GNUmakefile').write_text('unapproved')
                if failure == 'wrong_parent':
                    args[5] = root.parent.parent
                validate_base = p.validate_base
                def admit(base, approval, digest):
                    if failure == 'wrong_digest':
                        return validate_base(base, approval, digest)
                    return spec, dict(spec['environment'])
                with patch.object(p, 'load_manifest', return_value=data), \
                        patch.object(p, 'validate_base', side_effect=admit) as validate, \
                        patch.object(p.subprocess, 'run') as launch, self.assertRaises(p.AdmissionError):
                    p.execute(*args)
                launch.assert_not_called()
                if failure == 'wrong_parent':
                    validate.assert_not_called()
                self.assertFalse((root / 'logs' / 'claimed').exists())

    def test_private_build_generate_convert_reconstruct_check_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root, base, spec, data, setting, *_ = self.prepared(Path(tmp))
            calls = []
            self.execute(root, data, self.fake_run(root, calls))
            receipt = json.loads((root / 'receipt.json').read_text())
            self.assertEqual(receipt['state'], 'checked')
            self.assertFalse(receipt['runtime_scientifically_accepted'])
            python = Path(sys.executable).name
            self.assertEqual([Path(a[0][0]).name for a in calls],
                             [python, 'make', 'make', 'simc', 'make', 'make', 'ldd', 'make_root_tree',
                              'root', python])
            self.assertEqual(calls[0][0][:4], calls[-1][0][:4])
            self.assertEqual(calls[0][0][4:], ['probe'])
            self.assertEqual(calls[-1][0][4:], ['validate', str(root)])
            self.assertEqual((calls[0][1], calls[0][3]), (calls[-1][1], calls[-1][3]))
            self.assertEqual((calls[6][1], calls[6][3]), (calls[7][1], calls[7][3]))
            self.assertEqual(calls[6][0][-1], calls[7][0][0])
            self.assertEqual(receipt['converter_runtime'], receipt['root_runtime']['pyroot'])
            self.assertEqual(calls[-2][0][-1], 'recon_hcana.C("' + setting['canonical_stem'] + '", "kaon", 0.493677)')
            self.assertEqual(calls[3][2], setting['canonical_stem'] + '\n')
            self.assertEqual(p.inventory(base), spec['files'])
            self.assertTrue(receipt['digests'])

    def test_midrun_failure_incomplete_output_and_failed_validator(self):
        for failure in ('simc', 'missing-sidecar', 'stale-sidecar', 'validator'):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                root, base, spec, data, setting, *_ = self.prepared(Path(tmp))
                calls = []
                missing = 'OUTPUTS/' + setting['canonical_stem'] + '.geni' if failure == 'missing-sidecar' else None
                stub = self.fake_run(root, calls, fail=failure, missing=missing)
                def run(*args, **kwargs):
                    result = stub(*args, **kwargs)
                    if failure == 'stale-sidecar' and Path(args[0][0]).name == 'simc':
                        path = root / 'OUTPUTS' / (setting['canonical_stem'] + '.geni')
                        os.utime(path, ns=(1, 1))
                    return result
                with self.assertRaises((ValueError, subprocess.CalledProcessError)):
                    self.execute(root, data, run)
                receipt = json.loads((root / 'receipt.json').read_text())
                self.assertEqual(receipt['state'], 'failed')
                self.assertFalse(receipt['accepted_candidate'])
                self.assertIn('error', receipt)
                self.assertEqual(p.inventory(base), spec['files'])

    def test_runtime_identity_failures_stop_at_the_affected_command(self):
        for failure, message, last in (
            ('private-probe', 'before generation', Path(sys.executable).name),
            ('converter-core', 'differs from approved ROOT/PyROOT', 'ldd'),
        ):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                root, _, _, data, *_ = self.prepared(Path(tmp))
                calls = []
                with self.assertRaisesRegex(p.AdmissionError, message):
                    self.execute(root, data, self.fake_run(root, calls, fail=failure))
                self.assertEqual(Path(calls[-1][0][0]).name, last)
                if failure == 'private-probe':
                    self.assertEqual(len(calls), 1)
                self.assertFalse((root / 'logs' / 'converter.log').exists())
                receipt = json.loads((root / 'receipt.json').read_text())
                self.assertEqual(receipt['state'], 'failed')
                self.assertFalse(receipt['accepted_candidate'])

    def test_converter_loader_requirements(self):
        with tempfile.TemporaryDirectory() as tmp:
            core = Path(tmp) / 'libCore.so.6.24'
            core.write_bytes(b'approved core')
            alias = core.with_name('libCore.so')
            alias.symlink_to(core.name)
            identity = {'core': str(core), 'core_sha256': p.file_hash(core)}
            line = f'\tlibCore.so => {alias} (0xabcdef)\n'
            self.assertEqual(p.converter_core(line, identity), identity)
            for output in ('', line + line, line + '\tlibImt.so => not found\n'):
                with self.assertRaises(p.AdmissionError):
                    p.converter_core(output, identity)
            core.write_bytes(b'changed core')
            with self.assertRaisesRegex(p.AdmissionError, 'differs'):
                p.converter_core(line, identity)

    @unittest.skipUnless(sys.platform == 'linux' and shutil.which('gcc') and shutil.which('ldd'),
                         'Linux GCC and ldd required for native loader fixture')
    def test_real_loader_distinguishes_link_configuration_from_resolved_core(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            library = directory / 'core.c'
            library.write_text('int fixture_core(void) { return 0; }\n')
            main = directory / 'main.c'
            main.write_text('int fixture_core(void); int main(void) { return fixture_core(); }\n')
            approved, alternate = directory / 'approved', directory / 'alternate'
            for location in (approved, alternate):
                location.mkdir()
                subprocess.run(['gcc', '-shared', '-fPIC', '-Wl,-soname,libCore.so',
                                library, '-o', location / 'libCore.so'], check=True)
            executable = directory / 'converter-fixture'
            subprocess.run(['gcc', main, '-L' + str(approved), '-Wl,-rpath,' + str(approved),
                            '-lCore', '-o', executable], check=True)
            library = approved / 'libCore.so'
            expected = {'core': str(library), 'core_sha256': p.file_hash(library)}
            loader = [shutil.which('ldd'), '--', executable]
            env = {'PATH': os.environ['PATH'], 'LC_ALL': 'C', 'LD_LIBRARY_PATH': str(approved)}
            self.assertEqual(p.converter_core(p.command(loader, directory, env), expected), expected)
            env['LD_LIBRARY_PATH'] = str(alternate)
            with self.assertRaisesRegex(p.AdmissionError, 'differs'):
                p.converter_core(p.command(loader, directory, env), expected)

    def test_existing_outputs_owner_mismatch_and_source_tamper_block_before_launch(self):
        for failure in ('outputs', 'owner', 'source', 'link'):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                root, _, _, data, setting, *_ = self.prepared(Path(tmp))
                if failure == 'outputs':
                    (root / 'worksim/stale.root').write_text('stale')
                if failure == 'source':
                    (root / 'par.pl').write_text('altered')
                if failure == 'link':
                    (root / 'tmp/escape').symlink_to(Path(tmp) / 'base', target_is_directory=True)
                if failure == 'owner':
                    receipt = json.loads((root / 'receipt.json').read_text())
                    receipt['owner_uid'] += 1
                    p.write_json(root / 'receipt.json', receipt)
                args = self.execute_args(root, setting)
                spec = json.loads(args[3].read_text())
                with patch.object(p, 'load_manifest', return_value=data), \
                        patch.object(p, 'validate_base', return_value=(spec, dict(spec['environment']))), \
                        patch.object(p.subprocess, 'run') as launch, self.assertRaises(p.AdmissionError):
                    p.execute(*args)
                launch.assert_not_called()


@unittest.skipUnless(os.name == 'posix', 'Linux reference compatibility fixtures')
class ReferenceCompatibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        reference = os.environ.get('KSIGMA0_TEST_REFERENCE')
        if not reference:
            raise unittest.SkipTest('actual reference Git repository not configured')
        cls.reference = Path(reference)
        cls.tree = p.git_tree(cls.reference, p.REFERENCE_COMMIT)
        cls.raw_tree = subprocess.check_output(['git', '-C', reference, 'ls-tree', '-rz',
                                               '--full-tree', p.REFERENCE_COMMIT])
        oids = list(dict.fromkeys(record[1] for record in cls.tree.values()))
        batch = subprocess.run(['git', '-C', reference, 'cat-file', '--batch'],
                               input=('\n'.join(oids) + '\n').encode(), check=True,
                               stdout=subprocess.PIPE).stdout
        cls.blobs = {}
        offset = 0
        for oid in oids:
            end = batch.index(b'\n', offset)
            actual, kind, size = batch[offset:end].decode().split()
            assert actual == oid and kind == 'blob'
            size = int(size)
            cls.blobs[oid] = batch[end + 1:end + 1 + size]
            assert p.git_blob(cls.blobs[oid]) == oid
            offset = end + size + 2

    def full_base(self, directory):
        base, spec, approval, env = PrivateRunnerTests.fixture(self, directory)
        for name, (mode, oid) in self.tree.items():
            path = base / name
            path.parent.mkdir(parents=True, exist_ok=True)
            if mode == '120000':
                path.symlink_to(self.blobs[oid].decode())
            else:
                path.write_bytes(self.blobs[oid])
                path.chmod(int(mode, 8) & 0o777)
        optics = p.dependency_policy(self.tree)
        spec.update(commit=p.REFERENCE_COMMIT, files=p.inventory(base, optics), tracked=sorted(self.tree))
        approval.write_text(json.dumps(spec))
        return base, spec, approval, env

    def admit(self, base, spec, approval, env):
        real_command = p.command
        def command(args, cwd, env=None, input_text=None):
            if str(args[0]) == 'git':
                return spec['commit'] + '\n' if args[-1] == 'HEAD' else '\0'.join(sorted(self.tree)) + '\0'
            return real_command(args, cwd, env, input_text)
        # Real Git object inventory/bytes, with only metadata location redirected
        # to a disposable export. No git init/add/commit or reference-tree writes.
        with patch.dict(os.environ, env), patch.object(p, 'command', side_effect=command), \
                patch.object(p.subprocess, 'check_output', return_value=self.raw_tree), \
                patch.object(p, 'root_preflight', return_value={'stub': True}):
            return p.validate_base(base, approval, p.file_hash(approval))[0]

    def test_actual_553_entry_tree_and_source_defined_dependencies_admit(self):
        self.assertEqual(len(self.tree), 553)
        self.assertEqual({n for n, r in self.tree.items() if r[0] == '120000'}, set(p.OPTICS))
        for name, expected in p.REFERENCE_DEPENDENCIES.items():
            self.assertEqual(tuple(self.tree[name]), expected)
        self.assertIn(b"file_name='shms/shms_forward.dat'", self.blobs[self.tree['shared/transp.f'][1]])
        self.assertIn(b"file='shms/shms_recon.dat'", self.blobs[self.tree['shms/mc_shms_recon.f'][1]])
        with tempfile.TemporaryDirectory() as tmp:
            base, spec, approval, env = self.full_base(Path(tmp))
            admitted = self.admit(base, spec, approval, env)
            self.assertEqual(admitted['_git_tree'], self.tree)
            self.assertEqual(set(admitted['_optics_links']), set(p.OPTICS))

    def test_actual_optics_materialize_and_exact_archived_overlay_survives(self):
        archive = os.environ.get('KSIGMA0_TEST_ARCHIVE')
        if not archive:
            self.skipTest('original external archive not configured')
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            base, spec, approval, env = self.full_base(directory)
            admitted = self.admit(base, spec, approval, env)
            before = p.inventory(base, admitted['_optics_links'])
            data = p.load_manifest()
            payloads = p.archive_payloads(Path(archive), data)
            setting = data['settings'][18]
            parent = directory / 'attempts'
            parent.mkdir()
            root = p.prepare(base, parent, setting, payloads, data, admitted, p.file_hash(approval))
            for alias, target in admitted['_optics_links'].items():
                self.assertFalse((root / alias).is_symlink())
                self.assertEqual((root / alias).read_bytes(), (base / target).read_bytes())
                self.assertEqual((root / alias).stat().st_nlink, 1)
                self.assertNotEqual((root / alias).stat().st_ino, (base / target).stat().st_ino)
            self.assertFalse((root / p.BASE_LAUNCHER).exists())
            original = payloads[setting['settings_member'].rsplit('/', 1)[0] + '/recon_hcana.C']
            self.assertEqual((root / p.ARCHIVED_MACRO).read_bytes(), original)
            self.assertEqual((root / p.RECON_MACRO).read_bytes(), p.reconstruction_macro(
                original, data['common_sources']['recon_hcana.C']['archive_sha256']))
            self.assertEqual((root / p.STAGE['recon_hcana.h']).read_bytes(),
                             payloads[setting['settings_member'].rsplit('/', 1)[0] + '/recon_hcana.h'])
            self.assertEqual(p.inventory(base, admitted['_optics_links']), before)
            receipt = json.loads((root / 'receipt.json').read_text())
            (root / 'shms/shms_forward.dat').write_text('private optics changed')
            with self.assertRaisesRegex(p.AdmissionError, 'immutable private dependency'):
                p.private_check(root, receipt)

    def test_reference_exceptions_do_not_allow_arbitrary_outputs_or_links(self):
        for name, mode in (('cosy/new.bin', '100644'), ('recon_hcana/new.C', '100644'),
                           ('worksim/stale.root', '100644'), ('shms/other.dat', '120000')):
            tree = copy.deepcopy(self.tree)
            tree[name] = [mode, '0' * 40]
            with self.subTest(name=name), self.assertRaises(p.AdmissionError):
                p.dependency_policy(tree)
        for name in [*p.REFERENCE_DEPENDENCIES, *p.OPTICS]:
            tree = copy.deepcopy(self.tree)
            tree[name][1] = '0' * 40
            with self.subTest(name=name), self.assertRaises(p.AdmissionError):
                p.dependency_policy(tree)
        for name in p.OPTICS:
            tree = copy.deepcopy(self.tree)
            tree[name][0] = '100644'
            with self.subTest(regular_alias=name), self.assertRaises(p.AdmissionError):
                p.dependency_policy(tree)

    def test_reference_untracked_file_even_if_manifest_lists_it_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base, spec, approval, env = self.full_base(Path(tmp))
            (base / 'extra.txt').write_text('untracked')
            spec['files']['extra.txt'] = p.file_hash(base / 'extra.txt')
            approval.write_text(json.dumps(spec))
            with self.assertRaisesRegex(p.AdmissionError, 'untracked'):
                self.admit(base, spec, approval, env)

    def test_reference_modified_dependency_or_optics_target_fails_git_identity(self):
        for name in ('cosy/COSY.bin', 'recon_hcana/recon_hcana.C',
                     'shms/' + p.OPTICS['shms/shms_forward.dat'][1]):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmp:
                base, spec, approval, env = self.full_base(Path(tmp))
                with (base / name).open('ab') as stream:
                    stream.write(b'changed')
                spec['files'] = p.inventory(base, p.dependency_policy(self.tree))
                approval.write_text(json.dumps(spec))
                with self.assertRaisesRegex(p.AdmissionError, 'modified tracked'):
                    self.admit(base, spec, approval, env)

    def test_reference_optics_retarget_external_chain_or_regular_alias_rejected(self):
        for failure in ('external', 'internal-other', 'chain', 'regular'):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                base, _, _, _ = self.full_base(Path(tmp))
                alias = base / 'shms/shms_forward.dat'
                alias.unlink()
                if failure == 'external':
                    outside = Path(tmp) / 'outside'
                    outside.write_text('outside')
                    alias.symlink_to(outside)
                elif failure == 'internal-other':
                    alias.symlink_to(p.OPTICS['shms/shms_recon.dat'][1])
                elif failure == 'chain':
                    target = base / 'shms' / p.OPTICS['shms/shms_forward.dat'][1]
                    target.unlink()
                    target.symlink_to(p.OPTICS['shms/shms_recon.dat'][1])
                    alias.symlink_to(target.name)
                else:
                    alias.write_text('materialization in external base is not a Git symlink')
                with self.assertRaises(p.AdmissionError):
                    p.inventory(base, p.dependency_policy(self.tree))

    def root_fixture(self, directory):
        base, spec, approval, env = self.full_base(directory)
        prefix = directory / 'root-runtime'
        (prefix / 'include').mkdir(parents=True)
        (prefix / 'lib').mkdir()
        module_path = prefix / 'lib/ROOT/__init__.py'
        module_path.parent.mkdir()
        module_path.write_text('PyROOT fixture')
        core = prefix / 'lib/libCore.so'
        core.write_text('ROOT core fixture')
        config = prefix / 'bin/root-config'
        executable(config, '#!/bin/bash\ncase "$*" in\n*--libs*) echo "-L' + str(prefix) + '/lib -lCore";;\n'
                   '*--incdir*) echo "' + str(prefix) + '/include";;\n--prefix) echo "' + str(prefix) + '";;\n'
                   '--version) echo "6.24/08";;\n*) exit 1;;\nesac\n')
        root = prefix / 'bin/root'
        executable(root, '#!/bin/bash\nprintf "KSIGMA0_VERSION=62408\\nKSIGMA0_MAP=000-111 r-xp 0 00:00 1 ' + str(core) + '\\n"\n')
        spec['tools']['make']['path'] = str(Path(shutil.which('make')).resolve())
        spec['tools']['root-config']['path'], spec['tools']['root']['path'] = str(config), str(root)
        spec['environment']['ROOTSYS'] = str(prefix)
        spec['runtime_files'].update({str(core): p.file_hash(core), str(module_path): p.file_hash(module_path)})
        report = {'version_int': 62408, 'module': str(module_path), 'core': str(core),
                  'module_sha256': p.file_hash(module_path), 'core_sha256': p.file_hash(core)}
        return base, spec, dict(spec['environment']), report

    def probe(self, base, spec, env, report, missing=False):
        real_command = p.command
        self.probe_calls = []
        def command(args, cwd, env=None, input_text=None):
            self.probe_calls.append([str(x) for x in args])
            if str(args[0]) == sys.executable:
                if missing:
                    raise subprocess.CalledProcessError(1, args, 'No module named ROOT')
                self.assertEqual(args, p.validator_argv('probe'))
                return 'KSIGMA0_VALIDATOR=' + json.dumps(report) + '\n'
            return real_command(args, cwd, env, input_text)
        with patch.object(p, 'command', side_effect=command):
            return p.root_preflight(base, spec, env)

    def test_actual_converter_makefile_effective_root_override_and_lightweight_probes(self):
        with tempfile.TemporaryDirectory() as tmp:
            base, spec, env, report = self.root_fixture(Path(tmp))
            self.assertEqual(p.git_blob((base / 'util/root_tree/Makefile').read_bytes()),
                             '028d69661ea692d6805b8b0f4e3a64151b8c5f8e')
            result = self.probe(base, spec, env, report)
            self.assertEqual(result['converter']['ROOTSYS'], env['ROOTSYS'])
            self.assertEqual(result['converter']['ROOTDIR'], env['ROOTSYS'])
            self.assertEqual(result['pyroot'], report)
            self.assertEqual(result['version_int'], 62408)
            self.assertFalse(any(Path(a[0]).name in ('simc', 'run_simc_tree', 'make_root_tree') for a in self.probe_calls))
            self.assertTrue(any('--prefix=' + env['ROOTSYS'] in a for a in self.probe_calls))

    def test_missing_or_mismatched_pyroot_blocks_before_generation(self):
        for failure in ('missing', 'version', 'core', 'unapproved-module'):
            with self.subTest(failure=failure), tempfile.TemporaryDirectory() as tmp:
                base, spec, env, report = self.root_fixture(Path(tmp))
                if failure == 'version': report['version_int'] = 63004
                if failure == 'core': report['core'] = '/different/libCore.so'
                if failure == 'unapproved-module': spec['runtime_files'].pop(report['module'])
                with self.assertRaises((p.AdmissionError, subprocess.CalledProcessError)):
                    self.probe(base, spec, env, report, missing=failure == 'missing')
                self.assertFalse(any(Path(a[0]).name == 'simc' for a in self.probe_calls))

    def test_root_tool_prefix_mismatch_and_makefile_path_drift_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            base, spec, env, _ = self.root_fixture(Path(tmp))
            wrong = copy.deepcopy(spec)
            wrong['tools']['root-config']['path'] = wrong['tools']['g++']['path']
            with self.assertRaisesRegex(p.AdmissionError, 'installation mismatch'):
                p.converter_identity(base, wrong, env)
            with self.assertRaisesRegex(p.AdmissionError, 'absolute canonical'):
                p.converter_identity(base, spec, {**env, 'ROOTSYS': env['ROOTSYS'] + '/'})
            makefile = base / 'util/root_tree/Makefile'
            makefile.write_bytes(makefile.read_bytes().replace(b'$(ROOTDIR)/bin/root-config', b'root-config'))
            with self.assertRaisesRegex(p.AdmissionError, 'unsupported converter'):
                p.converter_identity(base, spec, env)

    def test_pyroot_admission_failure_creates_no_attempt_or_build(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            base, spec, approval, env = self.full_base(directory)
            parent = directory / 'attempts'
            parent.mkdir()
            archive = directory / 'original.zip'
            archive.write_text('not consumed: admission-only ordering fixture')
            args = ['--sigma-only', '--setting', p.load_manifest()['settings'][18]['canonical_stem'],
                    '--archive', str(archive), '--base', str(base), '--base-manifest', str(approval),
                    '--base-manifest-sha256', p.file_hash(approval), '--attempt-parent', str(parent)]
            with patch.object(p, 'archive_payloads', return_value={}), \
                    patch.object(p, 'validate_base', side_effect=p.AdmissionError('validator PyROOT unavailable')), \
                    patch.object(p, 'prepare') as prepare, patch.object(p.subprocess, 'run') as launch, \
                    self.assertRaisesRegex(p.AdmissionError, 'PyROOT unavailable'):
                p.main(args)
            prepare.assert_not_called()
            launch.assert_not_called()
            self.assertFalse(list(parent.iterdir()))


@unittest.skipUnless(os.name == 'posix' and shutil.which('bash'), 'Linux shell fixtures')
class ShellTests(unittest.TestCase):
    def fixture(self, tmp):
        root = Path(tmp) / 'repo'
        for name in ('background_samples/set_BackgroundInput.sh', 'src/setup/run_simc_tree',
                     'background_samples/background_samples.conf'):
            dst = root / name
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.write_bytes((REPO / name).read_bytes())
            dst.chmod(0o700)
        executable(root / 'farm_env/print_ltsep_path_fields.py',
                   'import os\nprint(os.environ["STUB_PATH_FIELDS"])\n')
        executable(root / 'background_samples/prepare_ksigma0.py',
                   'import sys\nprint("SIGMA-DISPATCH",sys.argv[1:])\n')
        return root

    def test_sigma_dispatch_bypasses_ltsep_and_shared_cleanup(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.fixture(tmp)
            (root / 'farm_env/print_ltsep_path_fields.py').unlink()
            result = subprocess.run(['bash', root / 'background_samples/set_BackgroundInput.sh',
                                     '--sigma-only', '--setting', 'exact'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('SIGMA-DISPATCH', result.stdout)
            self.assertFalse((root / 'background_samples/simc_gfortran').exists())
            self.assertFalse((root / 'log').exists())

    def test_missing_opt_in_and_runner_mixed_reaction_fail_early(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.fixture(tmp)
            for args in (['--setting', 'x'], ['low', '4p4', '2p74', '--sigma-only']):
                result = subprocess.run(['bash', root / 'background_samples/set_BackgroundInput.sh', *args],
                                        capture_output=True)
                self.assertNotEqual(result.returncode, 0)
            for args in (['/no', 'x', 'production'], ['/no', 'x', 'kaon']):
                result = subprocess.run(['bash', root / 'src/setup/run_simc_tree',
                                         '--ksigma0-isolated', *args], capture_output=True)
                self.assertNotEqual(result.returncode, 0)

    def test_legacy_neutron_delta_sidis_routing_and_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = self.fixture(tmp)
            simc = Path(tmp) / 'simc'
            volatile = Path(tmp) / 'volatile'
            (volatile / 'worksim').mkdir(parents=True)
            simc.mkdir()
            (simc / 'recon_hcana').mkdir()
            executable(simc / 'simc', '#!/bin/bash\nread -r stem\nprintf bin > "worksim/$stem.bin"\n')
            executable(simc / 'util/root_tree/make_root_tree',
                       '#!/bin/bash\nread -r stem\nprintf root > "../../worksim/$stem.root"\n')
            tools = Path(tmp) / 'tools'
            executable(tools / 'make', '#!/bin/bash\nexit 0\n')
            executable(tools / 'root', '#!/bin/bash\nprintf "%s\\n" "$@" >> "$ROOT_CALLS"\n')
            fields = ['unused'] * 18
            for i, val in {0: volatile, 13: 'Kaon', 14: 'stubuser', 16: simc, 17: root}.items():
                fields[i] = str(val)
            env = {**os.environ, 'STUB_PATH_FIELDS': ','.join(fields),
                   'PATH': str(tools) + ':' + os.environ['PATH'], 'ROOT_CALLS': str(Path(tmp) / 'calls')}
            conf = root / 'background_samples/background_samples.conf'
            conf.write_text('BACKGROUND_SIMC_PATH="' + str(simc) + '"\nRUN_NEUTRON=true\n'
                            'RUN_DELTA=true\nRUN_SIDIS=true\nBACKGROUND_REACTION="production"\n')
            for channel in ('neutron', 'delta', 'sidis'):
                for arm in ('left', 'center', 'right'):
                    inp = root / 'input/backgrounds' / channel / ('Prod_Coin_Q4p4W2p74' + arm + '_highe.inp')
                    inp.parent.mkdir(parents=True, exist_ok=True)
                    inp.write_text(channel + arm)
            result = subprocess.run(['bash', root / 'background_samples/set_BackgroundInput.sh',
                                     'high', '4p4', '2p74'], env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            for channel in ('neutron', 'delta', 'sidis'):
                out = volatile / 'OUTPUT/Analysis/KaonLT/background_samples' / channel
                self.assertEqual(len(list(out.glob('*.root'))), 3)
            self.assertEqual((simc / 'input').readlink(), root / 'input/backgrounds/sidis')
            calls = (Path(tmp) / 'calls').read_text()
            self.assertEqual(calls.count('"production"'), 9)
            self.assertNotIn('0.493677', calls)


if __name__ == '__main__':
    unittest.main()
