#!/usr/bin/env python3
"""Isolated, default-off Sigma0 generation owner (Python 3.9+, Linux).

The public manifest is the sole setting lookup. No colleague script is executed.
An operator must supply an independently approved base manifest AND its SHA256.
Base schema kaonlt_ksigma0_base/v1:
  commit: full git commit; files: {relative file: sha256} (complete deployed tree,
  excluding .git); tracked: sorted git ls-files inventory; tools: {name: {path,
  sha256, version_sha256}} for make/gfortran/gcc/g++/root/root-config/ldd. Versions
  are stdout+stderr from --version. environment: exact PATH, LD_LIBRARY_PATH,
  ROOTSYS, LIBRARY_PATH, CPATH, PYTHONPATH values (empty permitted). runtime_files: {absolute
  library/include: sha256}; active_reads: [relative data paths].
  files contains only tracked dependencies; the two approved optics aliases
  hash their target bytes. Their exact Git link/target objects are checked
  independently. ROOTSYS must select one approved installation; runtime_files
  must include the imported ROOT module and loaded libCore library.
The manifest is external approval evidence, never generated/approved by this CLI.
No base is inferred from the interface-reference checkout. Build commands are
fixed, not supplied shell text. All mutable files live in a fresh owned attempt.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import uuid
import zipfile

REPO = Path(__file__).resolve().parents[1]
MANIFEST = Path(__file__).with_name('ksigma0_inputs_manifest.json')
MANIFEST_SHA = '4b795a56666b96986f135f54917e4f4c377f76740204d0afa5c016140c15e5ef'
CSV_SHA = '02773236a846d0d8512c8204c889f880f9f4c322f8117648add20266c0332b6d'
TOOLS = ('make', 'gfortran', 'gcc', 'g++', 'root', 'root-config', 'ldd')
ENV_KEYS = ('PATH', 'LD_LIBRARY_PATH', 'ROOTSYS', 'LIBRARY_PATH', 'CPATH', 'PYTHONPATH')
MAX_MEMBER = 16 * 1024 * 1024
MAX_ARCHIVE = 128 * 1024 * 1024
STAGE = {'physics_iterate.f': 'physics_iterate.f', 'par.pl': 'par.pl',
         'mc_shms_hut.f': 'shms/mc_shms_hut.f',
         'recon_hcana.C': 'recon_hcana/recon_hcana.C',
         'recon_hcana.h': 'recon_hcana/recon_hcana.h'}

ATTEMPT_DIRS = ('input', 'worksim', 'OUTPUTS', 'recon_hcana', 'logs', 'provenance', 'tmp')
RECON_MACRO = STAGE['recon_hcana.C']
ARCHIVED_MACRO = 'provenance/recon_hcana.C'
ARCHIVED_DESTRUCTOR = (
    b'recon_hcana::~recon_hcana() {\n'
    b'\tif (mTree) {\n\t\tdelete mTree;\n\t}\n'
    b'\tif (mNewTree) {\n\t\tdelete mNewTree;\n\t}\n}\n'
)
FILE_OWNED_DESTRUCTOR = b'recon_hcana::~recon_hcana() = default;\n'
REFERENCE_COMMIT = '641b8c07182fd045fd8a5c007cdc1b71a418f3ce'
REFERENCE_DEPENDENCIES = {
    'recon_hcana/recon_hcana.C': ('100755', '0f112e752286b71bf3ae20c1aa00521f7d7e6ab5'),
    'recon_hcana/recon_hcana.h': ('100755', 'bb22c91e701e458dd2504aa501a5c50961fbf53c'),
    'recon_hcana/run_simc_tree': ('100755', '440dc4b1245be8d21b636f6905c4338dd326e626'),
    'cosy/COSY.bin': ('100644', '0d07ac61ecba51370d32b9fb30b06607b59e7427'),
    'cosy/JLAB.bin': ('100644', '7587b274ffbe6aca5d7fbd1e58e6adf540e44c1e'),
}
OPTICS = {
    'shms/shms_forward.dat': ('f6cf9c5bc0e9f82920debdcdd08f6e29ff9d9e92',
        'shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_hb80_forward.dat',
        '12888798d61680a1be354d54c802e33164bdeab1'),
    'shms/shms_recon.dat': ('bf2c618dc36415bf8883443d0d326970dee4e37f',
        'shms-2017-26cm-monte_q1_1018_q2_1027_q3_1018_hb80_recon_fit.dat',
        '053d04a7223a5d3a0ece1dd19012802ca1246b01'),
}
BASE_LAUNCHER = 'recon_hcana/run_simc_tree'


class AdmissionError(ValueError):
    pass


def require(ok, message):
    if not ok:
        raise AdmissionError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def relative(name):
    require(isinstance(name, str) and name and '\\' not in name and
            not re.search(r'[:\x00-\x1f]', name), 'unsafe relative path')
    p = PurePosixPath(name)
    require(not p.is_absolute() and all(x not in ('', '.', '..') for x in name.split('/')),
            'unsafe relative path: ' + name)
    return p


def physical(path, exists=True):
    p = Path(path).absolute()
    # Reject aliases even when a later component is missing.
    for part in (p, *p.parents):
        require(not part.is_symlink(), 'symlink path: ' + str(part))
    require(p == p.resolve(strict=exists), 'ambiguous physical path: ' + str(p))
    return p


def contains(parent, child):
    return parent == child or parent in child.parents


def disjoint(a, b):
    require(not contains(a, b) and not contains(b, a), 'overlapping protected/private roots')


def owned(path):
    require(path.stat().st_uid == os.getuid(), 'private root owner mismatch')


def load_manifest():
    require(file_hash(MANIFEST) == MANIFEST_SHA, 'reviewed manifest hash mismatch')
    data = json.loads(MANIFEST.read_text())
    require(len(data['settings']) == 25, 'reviewed 25-setting inventory missing')
    return data


def select(data, stems):
    require(stems and len(set(stems)) == len(stems), 'empty/duplicate settings')
    lookup = {s['canonical_stem']: s for s in data['settings']}
    require(all(s in lookup for s in stems), 'unknown setting (including Right/low)')
    return [lookup[s] for s in stems]


def validate_input(raw, setting_raw, data):
    values = {}
    for line in raw.decode('ascii').splitlines():
        line = line.split(';', 1)[0].strip()
        match = re.fullmatch(r'([\w%()]+)\s*=\s*([-+\d.eEdD]+)', line)
        if match:
            key, value = match.groups()
            require(key.lower() not in values, 'duplicate input key')
            values[key.lower()] = float(value.replace('D', 'e').replace('d', 'e'))
    for key, expected in data['flags'].items():
        require(values.get(key.lower()) == expected, 'wrong input flag: ' + key)
    settings = json.loads(setting_raw)
    for key, expected in data['flags'].items():
        if key in settings['program_settings']:
            require(settings['program_settings'][key] == expected, 'wrong settings flag: ' + key)
    require(settings['recon_hcana_settings'] == {'mode': 'kaon', 'hadron_mass': 0.493677},
            'incompatible reconstruction settings')


def archive_payloads(archive, data):
    """Admit outer bytes before opening; bound and hash every consumed member."""
    require(file_hash(archive) == data['sha256'], 'archive SHA256 mismatch')
    payloads = {}
    with zipfile.ZipFile(archive) as z:
        infos = z.infolist()
        require(len(infos) <= 2000, 'archive inventory bound')
        names = [i.filename for i in infos]
        require(len(set(names)) == len(names), 'duplicate archive member')
        require(sum(i.file_size for i in infos) <= MAX_ARCHIVE, 'archive size bound')
        for i in infos:
            relative(i.filename.rstrip('/'))
            mode = i.external_attr >> 16
            require(not stat.S_ISLNK(mode) and (not stat.S_IFMT(mode) or
                    stat.S_ISREG(mode) or stat.S_ISDIR(mode)), 'nonregular archive member')
            require(not i.flag_bits & 1 and i.file_size <= MAX_MEMBER, 'archive member bound/encryption')

        def read(name, expected):
            require(name in names and not z.getinfo(name).is_dir(), 'missing member: ' + name)
            with z.open(name) as stream:
                b = stream.read(MAX_MEMBER + 1)
            require(len(b) <= MAX_MEMBER and digest(b) == expected, 'member hash mismatch: ' + name)
            payloads[name] = b
            return b

        # All 25 are admitted, even for one selected attempt. Historical reports
        # are verified provenance only; they are never staged as fresh output.
        for s in data['settings']:
            inp = read(s['input_member'], s['input_sha256'])
            cfg = read(s['settings_member'], s['settings_sha256'])
            read(s['historical_hist_member'], s['historical_hist_sha256'])
            validate_input(inp, cfg, data)
            prefix = s['settings_member'].rsplit('/', 1)[0] + '/'
            for name, spec in data['common_sources'].items():
                read(prefix + name, spec['archive_sha256'])
            for name, expected in data['common_tables'].items():
                read(prefix + name, expected)
    return payloads


def git_blob(data):
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def git_tree(base, commit):
    raw = subprocess.check_output(['git', '-C', str(base), 'ls-tree', '-rz', '--full-tree', commit])
    tree = {}
    for entry in raw.split(b'\0'):
        if not entry:
            continue
        meta, name = entry.split(b'\t')
        mode, kind, oid = meta.decode('ascii').split()
        name = name.decode('utf-8')
        relative(name)
        require(kind == 'blob' and mode in ('100644', '100755', '120000') and
                name not in tree, 'unsupported Git tree entry')
        tree[name] = [mode, oid]
    require(tree, 'empty Git base tree')
    return tree


def dependency_policy(tree):
    """Only the user-reviewed tracked exceptions; no suffix/directory exemption."""
    optics = {}
    for name, identity in tree.items():
        if name in REFERENCE_DEPENDENCIES:
            require(tuple(identity) == REFERENCE_DEPENDENCIES[name], 'reference dependency Git identity: ' + name)
        if identity[0] == '120000' or name in OPTICS:
            require(name in OPTICS, 'unapproved tracked symlink: ' + name)
            oid, target, target_oid = OPTICS[name]
            target_name = 'shms/' + target
            require(tuple(identity) == ('120000', oid) and
                    tuple(tree.get(target_name, ())) == ('100644', target_oid),
                    'optics Git link/target identity: ' + name)
            optics[name] = target_name
        require(Path(name).suffix.lower() not in ('.root', '.bin', '.hist', '.gen', '.geni') or
                name in REFERENCE_DEPENDENCIES, 'stale output in base: ' + name)
        require(name.split('/')[0] not in ATTEMPT_DIRS or name in REFERENCE_DEPENDENCIES,
                'base contains mutable job files: ' + name)
    return optics


def inventory(root, optics=None):
    files = {}
    optics = optics or {}
    for p in root.rglob('*'):
        if p.relative_to(root).parts[0] == '.git':
            continue
        name = p.relative_to(root).as_posix()
        if p.is_symlink():
            require(name in optics, 'link in deployed tree: ' + str(p))
            target = root / optics[name]
            require(os.readlink(p) == target.name and physical(target).is_file() and
                    contains(root, target) and p.resolve(strict=True) == target and
                    target.stat().st_nlink == 1, 'unsafe optics target: ' + name)
            files[name] = file_hash(target)
            continue
        require(name not in optics, 'approved optics path is not a symlink: ' + name)
        if p.is_file():
            require(p.stat().st_nlink == 1, 'hardlink in deployed tree')
            files[p.relative_to(root).as_posix()] = file_hash(p)
        else:
            require(p.is_dir(), 'nonregular deployed file')
    return files


def command(args, cwd, env=None, input_text=None):
    return subprocess.run([str(x) for x in args], cwd=cwd, env=env, input=input_text,
                          text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          check=True).stdout


def converter_make(spec):
    # Source ROOTSYS is hardcoded. A command-line assignment takes precedence;
    # ROOTDIR/root-config remain the source-defined expressions, checked below.
    return [spec['tools']['make']['path'], '-f', 'Makefile',
            'ROOTSYS=' + spec['environment']['ROOTSYS']]


def validator_argv(mode, *args):
    script = REPO / 'src/setup/validate_ksigma0_output.py'
    return [sys.executable, '-B', '-s', script, mode, *args]


def validator_report(output):
    reports = re.findall(r'^KSIGMA0_VALIDATOR=(.*)$', output, re.M)
    require(len(reports) == 1, 'validator identity/result report unavailable or ambiguous')
    report = json.loads(reports[0])
    require(isinstance(report, dict), 'invalid validator identity/result report')
    return report


def converter_core(output, expected):
    require(not re.search(r'=>\s*not found\b', output), 'converter runtime dependency unavailable')
    matches = re.findall(r'^[ \t]*libCore\.so(?:\.[^\s]+)?[ \t]+=>[ \t]+'
                         r'(/.*?)[ \t]+\(0x[0-9a-fA-F]+\)[ \t]*$',
                         output, re.M)
    require(len(matches) == 1, 'converter libCore resolution unavailable or ambiguous')
    path = Path(matches[0]).resolve(strict=True)
    identity = {'core': str(path), 'core_sha256': file_hash(path)}
    require(all(identity[key] == expected[key] for key in identity),
            'converter loaded libCore differs from approved ROOT/PyROOT')
    return identity


def converter_identity(base, spec, env):
    prefix = physical(env['ROOTSYS'])
    require(env['ROOTSYS'] == str(prefix), 'ROOTSYS must be an absolute canonical path')
    require(re.fullmatch(r'[A-Za-z0-9_/+.\-]+', str(prefix)) is not None,
            'ROOTSYS cannot be represented safely by the native Makefile')
    for tool in ('root-config', 'root'):
        require(physical(prefix / 'bin' / tool) == physical(spec['tools'][tool]['path']),
                'ROOTSYS/tool installation mismatch: ' + tool)
    makefile = base / 'util/root_tree/Makefile'
    source = ''.join(line.split('#', 1)[0] for line in makefile.read_text().splitlines())
    compact = re.sub(r'\s+', '', source)
    require('ROOTDIR=$(ROOTSYS)' in compact and
            '$(shell$(ROOTDIR)/bin/root-config--prefix=$(ROOTDIR)--libs)' in compact and
            '$(shell$(ROOTDIR)/bin/root-config--prefix=$(ROOTDIR)--incdir)' in compact,
            'unsupported converter ROOTSYS/root-config interface')
    query = '\n'.join('$(info KSIGMA0_' + key + '=$(' + key + '))'
                      for key in ('ROOTSYS', 'ROOTDIR', 'ROOTLIBS', 'ROOTINCLUDE'))
    query += '\n.PHONY: kaonlt-root-identity\nkaonlt-root-identity: ;\n'
    output = command(converter_make(spec) + ['--no-print-directory', '-f', '-', 'kaonlt-root-identity'],
                     makefile.parent, env, query)
    effective = {}
    for key in ('ROOTSYS', 'ROOTDIR', 'ROOTLIBS', 'ROOTINCLUDE'):
        values = re.findall(r'^KSIGMA0_' + key + r'=(.*)$', output, re.M)
        require(len(values) == 1, 'effective converter identity missing: ' + key)
        effective[key] = values[0].strip()
    require(effective['ROOTSYS'] == effective['ROOTDIR'] == str(prefix),
            'effective converter ROOTSYS mismatch')
    config = spec['tools']['root-config']['path']
    require(command([config, '--prefix'], base, env).strip() == str(prefix), 'root-config prefix mismatch')
    libs = command([config, '--prefix=' + str(prefix), '--libs'], base, env).strip()
    incdir = command([config, '--prefix=' + str(prefix), '--incdir'], base, env).strip()
    require(effective['ROOTLIBS'] == libs and effective['ROOTINCLUDE'] == '-I ' + incdir,
            'effective converter ROOT libraries/includes mismatch')
    require(contains(prefix, physical(incdir)), 'ROOT includes outside approved installation')
    version = command([config, '--version'], base, env).strip()
    match = re.fullmatch(r'(\d+)\.(\d+)[/.](\d+)', version)
    require(match, 'unrecognized ROOT version')
    major, minor, patch = map(int, match.groups())
    return prefix, major * 10000 + minor * 100 + patch, effective


def root_preflight(base, spec, env):
    """Read-only configuration/probes, before any build or generation.

    Do not run the native converter or reconstruction macro here. Probe homes
    and temporary files are disposable; the external base stays read-only.
    """
    with tempfile.TemporaryDirectory(prefix='ksigma0-root-probe-') as directory:
        probe_env = {**env, 'HOME': directory, 'TMPDIR': directory}
        prefix, version, effective = converter_identity(base, spec, probe_env)
        expression = ('{ printf("KSIGMA0_VERSION=%d\\n", gROOT->GetVersionInt()); '
                      'FILE *f=fopen("/proc/self/maps","r"); char b[4096]; '
                      'while(f && fgets(b,sizeof(b),f)) if(strstr(b,"/libCore.so")) '
                      'printf("KSIGMA0_MAP=%s",b); if(f) fclose(f); }')
        root_output = command([spec['tools']['root']['path'], '-n', '-l', '-b', '-q', '-e', expression],
                              directory, probe_env)
        versions = re.findall(r'^KSIGMA0_VERSION=(\d+)$', root_output, re.M)
        require(versions == [str(version)], 'reconstruction ROOT version mismatch')
        cores = {str(physical(line.split(None, 5)[5]))
                 for line in re.findall(r'^KSIGMA0_MAP=(.*)$', root_output, re.M)}
        require(len(cores) == 1, 'reconstruction ROOT core library unavailable/ambiguous')
        pyroot = validator_report(command(validator_argv('probe'), directory, probe_env))
        require(pyroot['version_int'] == version and pyroot['core'] in cores,
                'validator PyROOT/ROOT version or core mismatch')
        for name in ('module', 'core'):
            path = physical(pyroot[name])
            require(contains(prefix, path) and file_hash(path) == pyroot[name + '_sha256'] ==
                    spec['runtime_files'].get(str(path)), 'unapproved PyROOT/ROOT runtime: ' + name)
        return {'prefix': str(prefix), 'version_int': version, 'pyroot': pyroot,
                'converter': effective, 'makefile_sha256': file_hash(base / 'util/root_tree/Makefile'),
                'python_sha256': file_hash(Path(sys.executable).resolve())}


def validate_base(base, manifest_path, expected):
    require(re.fullmatch('[0-9a-f]{64}', expected or '') is not None, 'base approval digest required')
    require(file_hash(manifest_path) == expected, 'approved base manifest hash mismatch')
    spec = json.loads(manifest_path.read_text())
    require(spec['schema'] == 'kaonlt_ksigma0_base/v1', 'base manifest schema')
    require(re.fullmatch('[0-9a-f]{40}', spec['commit']) is not None, 'base commit required')
    require(command(['git', '-C', base, 'rev-parse', 'HEAD'], base).strip() == spec['commit'],
            'base commit mismatch')
    tree = git_tree(base, spec['commit'])
    optics = dependency_policy(tree)
    actual = inventory(base, optics)
    require(actual == spec['files'], 'deployed base inventory/hash mismatch')
    tracked = command(['git', '-C', base, 'ls-files', '-z'], base).rstrip('\0').split('\0')
    require(sorted(tracked) == spec['tracked'] == sorted(tree) and set(actual) == set(tree),
            'tracked base inventory mismatch/untracked file')
    for name, (mode, oid) in tree.items():
        raw = os.fsencode(os.readlink(base / name)) if mode == '120000' else (base / name).read_bytes()
        require(git_blob(raw) == oid, 'modified tracked base source: ' + name)
    for name in ('Makefile', 'modules.f', 'simulate.inc', 'dbase.f', 'simc.f',
                 'results_write.f', 'shms/mc_shms_hut.f', 'util/root_tree/Makefile',
                 'util/root_tree/make_root_tree.f'):
        require(name in actual, 'missing base dependency: ' + name)
    make = (base / 'Makefile').read_text()
    require(any('physics_iterate.o' in l for l in make.splitlines() if not l.lstrip().startswith('#')),
            'base does not build physics_iterate.o')
    db = (base / 'dbase.f').read_text().lower()
    require('which_kaon.eq.1' in db and 'sigma0' in db, 'Sigma0 recoil interface absent')
    reads = set(spec['active_reads'])
    require({'par.pl', 'weightc.dat', 'weightd.dat', 'saghai_proton.dat', 'saghai_sigma0.dat'} <= reads,
            'unapproved active data reads')
    for name in reads:
        relative(name)
    for name in ('weightc.dat', 'weightd.dat', 'saghai_proton.dat', 'saghai_sigma0.dat'):
        require("file='" + name + "'" in db.replace(' ', ''), 'startup data interface mismatch')
    require(set(spec['tools']) == set(TOOLS), 'toolchain identity incomplete')
    require(set(spec['environment']) == set(ENV_KEYS), 'runtime environment identity incomplete')
    env = {'HOME': '', 'LANG': 'C', 'LC_ALL': 'C'}
    for key in ENV_KEYS:
        require(os.environ.get(key, '') == spec['environment'][key], 'environment mismatch: ' + key)
        env[key] = spec['environment'][key]
    for name, record in spec['tools'].items():
        tool = physical(record['path'])
        require(file_hash(tool) == record['sha256'] and os.access(tool, os.X_OK), 'tool identity: ' + name)
        require(digest(command([tool, '--version'], base, env).encode()) == record['version_sha256'],
                'tool version: ' + name)
        require(Path(shutil.which(name, path=env['PATH']) or '').resolve() == tool,
                'PATH tool mismatch: ' + name)
    require(spec['runtime_files'], 'deployed runtime/library identity required')
    for name, expected_hash in spec['runtime_files'].items():
        require(file_hash(physical(name)) == expected_hash, 'runtime dependency mismatch')
    spec['_git_tree'] = tree
    spec['_optics_links'] = optics
    spec['_root_runtime'] = root_preflight(base, spec, env)
    return spec, env


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + '\n')


def reconstruction_macro(raw, expected):
    require(digest(raw) == expected, 'archived reconstruction macro hash mismatch')
    require(raw.count(ARCHIVED_DESTRUCTOR) == 1, 'archived tree destructor does not match')
    # Process() leaves both trees owned by its local TFiles.
    return raw.replace(ARCHIVED_DESTRUCTOR, FILE_OWNED_DESTRUCTOR, 1)


def prepare(base, parent, setting, payloads, data, spec, base_digest):
    """Called only after entire batch admission; mkdir is exclusive, never reused."""
    stem = setting['canonical_stem']
    root = parent / (stem + '-' + uuid.uuid4().hex)
    root.mkdir(mode=0o700)
    owned(root)
    try:
        for name, expected in spec['files'].items():
            if name == BASE_LAUNCHER:
                continue  # Admitted Git metadata, never a private executable.
            source = base / name
            if name in spec.get('_optics_links', {}):
                source = physical(base / spec['_optics_links'][name])
                require(os.readlink(base / name) == source.name and contains(base, source),
                        'optics alias changed before materialization')
            dst = root / str(relative(name))
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, dst)
            dst.chmod(source.stat().st_mode & 0o777)
            require(file_hash(dst) == expected, 'base changed during copy')
        for d in ATTEMPT_DIRS:
            (root / d).mkdir(exist_ok=True)
        prefix = setting['settings_member'].rsplit('/', 1)[0] + '/'
        loaded = {}
        mapping = dict(STAGE)
        mapping.update({name: name for name in data['common_tables']})
        for name, dest in mapping.items():
            b = payloads[prefix + name]
            if dest == RECON_MACRO:
                (root / ARCHIVED_MACRO).write_bytes(b)
                loaded[ARCHIVED_MACRO] = digest(b)
                b = reconstruction_macro(b, data['common_sources'][name]['archive_sha256'])
            (root / dest).write_bytes(b)
            loaded[dest] = digest(b)
        (root / 'input' / (stem + '.inp')).write_bytes(payloads[setting['input_member']])
        loaded['input/' + stem + '.inp'] = setting['input_sha256']
        for name in spec['active_reads']:
            require((root / name).is_file(), 'active data read unavailable: ' + name)
            require(name in loaded or name in spec['files'], 'unapproved active data path')
        for name in ('physics_iterate.f', 'par.pl', 'shms/mc_shms_hut.f'):
            shutil.copyfile(root / name, root / 'provenance' / Path(name).name)
        receipt = {'schema': 'kaonlt_ksigma0_attempt/v1', 'state': 'prepared',
                   'accepted_candidate': False, 'runtime_scientifically_accepted': False,
                   'root': str(root), 'owner_uid': os.getuid(), 'stem': stem,
                   'manifest_sha256': MANIFEST_SHA, 'archive_sha256': data['sha256'],
                   'base_manifest_sha256': base_digest, 'base_commit': spec['commit'],
                   'loaded': loaded, 'historical_reference': setting,
                   'base_files': spec['files'], 'tools': spec['tools'],
                   'environment': spec['environment'], 'runtime_files': spec['runtime_files'],
                   'base_git_tree': spec.get('_git_tree', {}),
                   'optics_links': spec.get('_optics_links', {}),
                   'root_runtime': spec.get('_root_runtime', {}),
                   'python': {'executable': sys.executable, 'version': sys.version,
                              'sha256': file_hash(Path(sys.executable).resolve())}}
        write_json(root / 'receipt.json', receipt)
        return root
    except BaseException as exc:
        write_json(root / 'failed-preparation.json', {'state': 'failed', 'error': str(exc)})
        raise


def private_check(root, receipt):
    require(str(physical(root)) == receipt['root'], 'attempt identity mismatch')
    owned(root)
    require(receipt['owner_uid'] == os.getuid(), 'attempt owner mismatch')
    for p in root.rglob('*'):
        require(not p.is_symlink() and contains(root, p.resolve()), 'escaping attempt path')
        if p.is_file():
            require(p.stat().st_nlink == 1, 'mutable hardlink')
    for name, expected in receipt['loaded'].items():
        require(file_hash(root / str(relative(name))) == expected, 'staged source changed: ' + name)
    dependencies = set(receipt.get('optics_links', {})) | set(receipt.get('optics_links', {}).values())
    dependencies.update(n for n in REFERENCE_DEPENDENCIES if n.startswith('cosy/') and n in receipt['base_files'])
    for name in dependencies:
        require(file_hash(root / str(relative(name))) == receipt['base_files'][name],
                'immutable private dependency changed: ' + name)


def execute(root, stem, base, approval, approved_digest, parent, archive):
    """Re-admit the external approval before launching the private build."""
    root, base, approval, parent, archive = map(physical, (root, base, approval, parent, archive))
    owned(root)
    owned(parent)
    require(root.parent == parent and stat.S_IMODE(root.stat().st_mode) == 0o700,
            'attempt must be an owned private child')
    for protected in (REPO, base, archive, approval):
        disjoint(parent, protected)
    disjoint(base, REPO)
    require(not any(k.startswith('LT_BG_SIGMA0_') and value for k, value in os.environ.items()),
            'consumer bindings must remain unset during generation')
    receipt_path = physical(root / 'receipt.json')
    require(receipt_path.stat().st_nlink == 1, 'attempt receipt must be a private file')
    receipt = json.loads(receipt_path.read_text())
    require(receipt['state'] == 'prepared' and receipt['stem'] == stem and
            re.fullmatch(r'Prod_Coin_Q[0-9]p[0-9]W[0-9]p[0-9]{2}(left|center|right)_(lowe|highe)', stem),
            'invalid/reused private attempt')
    require(re.fullmatch(re.escape(stem) + r'-[0-9a-f]{32}', root.name) is not None,
            'invalid attempt directory name')
    private_check(root, receipt)
    data = load_manifest()
    setting = select(data, [stem])[0]
    require(receipt['manifest_sha256'] == MANIFEST_SHA and
            receipt['archive_sha256'] == data['sha256'] == file_hash(archive) and
            receipt['historical_reference'] == setting, 'attempt provenance mismatch')
    spec, env = validate_base(base, approval, approved_digest)
    for name, value in {
        'base_manifest_sha256': approved_digest, 'base_commit': spec['commit'],
        'base_files': spec['files'], 'tools': spec['tools'],
        'environment': spec['environment'], 'runtime_files': spec['runtime_files'],
        'base_git_tree': spec.get('_git_tree', {}),
        'optics_links': spec.get('_optics_links', {}),
        'root_runtime': spec.get('_root_runtime', {}),
    }.items():
        require(receipt[name] == value, 'attempt approval mismatch: ' + name)
    expected_loaded = {dest: data['common_sources'][name]['archive_sha256']
                       for name, dest in STAGE.items()}
    archived_hash = expected_loaded[RECON_MACRO]
    expected_loaded[ARCHIVED_MACRO] = archived_hash
    expected_loaded[RECON_MACRO] = digest(reconstruction_macro(
        (root / ARCHIVED_MACRO).read_bytes(), archived_hash))
    expected_loaded.update(data['common_tables'])
    expected_loaded['input/' + stem + '.inp'] = setting['input_sha256']
    require(receipt['loaded'] == expected_loaded, 'unapproved staged file inventory')
    expected_files = {**spec['files'], **expected_loaded}
    expected_files.pop(BASE_LAUNCHER, None)
    for name in ('physics_iterate.f', 'par.pl', 'shms/mc_shms_hut.f'):
        expected_files['provenance/' + Path(name).name] = expected_loaded[name]
    expected_files['receipt.json'] = file_hash(receipt_path)
    require(not (root / '.git').exists() and inventory(root) == expected_files,
            'unapproved attempt files or existing outputs')
    for d in ATTEMPT_DIRS:
        require((root / d).is_dir(), 'missing attempt directory: ' + d)
    env.update(HOME=str(root), TMPDIR=str(root / 'tmp'))
    (root / 'logs' / 'claimed').touch(exist_ok=False)
    receipt.update(state='running', started_ns=time.time_ns(), commands=[])
    write_json(receipt_path, receipt)

    def run(argv, cwd, label, text=None):
        private_check(root, receipt)
        receipt['commands'].append({'argv': [str(x) for x in argv], 'cwd': str(cwd)})
        write_json(receipt_path, receipt)
        log = root / 'logs' / (label + '.log')
        with log.open('x') as stream:
            subprocess.run([str(x) for x in argv], cwd=cwd, env=env, input=text,
                           text=True, stdout=stream, stderr=subprocess.STDOUT, check=True)
        private_check(root, receipt)
        return log

    def artifact(path):
        require(path.is_file() and path.stat().st_size > 0 and
                path.stat().st_mtime_ns >= receipt['started_ns'], 'missing/stale artifact: ' + str(path))

    try:
        pyroot = receipt['root_runtime']['pyroot']
        probe = run(validator_argv('probe'), root, 'validator-preflight')
        require(validator_report(probe.read_text()) == pyroot,
                'private validator PyROOT differs from approved runtime before generation')
        make = receipt['tools']['make']['path']
        run([make, 'clean'], root, 'generator-clean')
        run([make], root, 'generator-build')
        artifact(root / 'physics_iterate.o')
        artifact(root / 'simc')
        require(os.access(root / 'simc', os.X_OK), 'private generator not executable')
        receipt['build_artifacts'] = {n: file_hash(root / n) for n in ('simc', 'physics_iterate.o')}
        run([root / 'simc'], root, 'generator', stem + '\n')
        artifact(root / 'worksim' / (stem + '.bin'))
        for ext in ('.hist', '.gen', '.geni'):
            artifact(root / 'OUTPUTS' / (stem + ext))
        converter = root / 'util/root_tree'
        run(converter_make(spec) + ['clean'], converter, 'converter-clean')
        run(converter_make(spec), converter, 'converter-build')
        executable = converter / 'make_root_tree'
        artifact(executable)
        require(os.access(executable, os.X_OK), 'private converter not executable')
        receipt['build_artifacts']['util/root_tree/make_root_tree'] = file_hash(executable)
        loader_log = run([receipt['tools']['ldd']['path'], '--', executable], converter, 'converter-loader')
        receipt['converter_runtime'] = converter_core(loader_log.read_text(), pyroot)
        require(file_hash(executable) == receipt['build_artifacts']['util/root_tree/make_root_tree'],
                'private converter changed after loader check')
        run([executable], converter, 'converter', stem + '\n')
        raw = root / 'worksim' / (stem + '.root')
        artifact(raw)
        recon = root / 'recon_hcana'
        shutil.copyfile(raw, recon / raw.name)
        shutil.copyfile(root / 'OUTPUTS' / (stem + '.hist'), recon / (stem + '.hist'))
        invocation = 'recon_hcana.C("' + stem + '", "kaon", 0.493677)'
        run([receipt['tools']['root']['path'], '-n', '-q', '-l', '-b', invocation], recon, 'reconstruction')
        artifact(recon / (stem + '_recon.root'))
        validated = run(validator_argv('validate', root), root, 'validator')
        result = validator_report(validated.read_text())
        receipt.update(state='checked', accepted_candidate=True, validation=result)
    except BaseException as exc:
        receipt.update(state='failed', accepted_candidate=False, error=str(exc))
        raise
    finally:
        receipt['digests'] = {p.relative_to(root).as_posix(): file_hash(p)
                              for d in ('logs', 'OUTPUTS', 'worksim', 'recon_hcana', 'provenance')
                              for p in (root / d).rglob('*') if p.is_file() and not p.is_symlink()}
        write_json(receipt_path, receipt)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sigma-only', action='store_true', required=True)
    p.add_argument('--setting', action='append', required=True, help='exact canonical stem; repeat for explicit batch')
    p.add_argument('--archive', required=True)
    p.add_argument('--base', required=True)
    p.add_argument('--base-manifest', required=True)
    p.add_argument('--base-manifest-sha256', required=True)
    p.add_argument('--attempt-parent', required=True)
    args = p.parse_args(argv)
    require(os.name == 'posix', 'Linux/WSL required')
    require(not any(k.startswith('LT_BG_SIGMA0_') and v for k, v in os.environ.items()),
            'consumer bindings must remain unset during generation')
    data = load_manifest()
    settings = select(data, args.setting)
    archive, base, base_manifest, parent = [physical(x) for x in
                                          (args.archive, args.base, args.base_manifest, args.attempt_parent)]
    owned(parent)
    require(parent.is_dir() and base.is_dir(), 'base/attempt parent must exist')
    for protected in (REPO, base, archive, base_manifest):
        disjoint(parent, protected)
    disjoint(base, REPO)
    payloads = archive_payloads(archive, data)
    spec, _ = validate_base(base, base_manifest, args.base_manifest_sha256)
    for s in settings:
        root = prepare(base, parent, s, payloads, data, spec, args.base_manifest_sha256)
        try:
            subprocess.run(['bash', str(REPO / 'src/setup/run_simc_tree'), '--ksigma0-isolated', str(root),
                            s['canonical_stem'], 'kaon', str(base), str(base_manifest),
                            args.base_manifest_sha256, str(parent), str(archive)], check=True)
        finally:
            require(inventory(base, spec.get('_optics_links')) == spec['files'],
                    'protected base changed during attempt')
        print(root / 'receipt.json')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (AdmissionError, OSError, KeyError, ValueError, subprocess.CalledProcessError,
            zipfile.BadZipFile) as exc:
        print('KSigma0 BLOCKED: ' + str(exc), file=sys.stderr)
        sys.exit(1)
