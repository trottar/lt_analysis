#!/usr/bin/env python3
"""One matched-cohort acquisition owner. No automatic retry or science acceptance.

The private-output candidate must pass the complete writer audit before a child
is enabled. This module deliberately exposes that gate separately from capture.
"""
from __future__ import annotations

import ast
import argparse
from contextlib import redirect_stdout
from copy import deepcopy
import hashlib
import glob
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time
import uuid
import zipfile

try:
    from . import run_pion_t_template_comparison_gate as ab
    from . import attest_pion_chain_ltsep_default_root as resolver
    from . import pion_local_t_event_validation_capture as capture
    from . import run_pion_local_t_event_validation_gate as event_owner
except ImportError:
    import run_pion_t_template_comparison_gate as ab
    import attest_pion_chain_ltsep_default_root as resolver
    import pion_local_t_event_validation_capture as capture
    import run_pion_local_t_event_validation_gate as event_owner

isolation = ab.isolation
BASE = "b23d231ce41f96759a2d9cddd39bfa20f8953ad2"
AUTH_ENV = "LT_ANALYSIS_PION_MATCHED_COHORT_MANIFEST"
PROFILE_ID = "pion_local_t_matched_cohort/v1"
OWNER = "testing/run_pion_local_t_matched_cohort_gate.py"
SETTINGS = ("Center", "Left")

# Observe Python opens of only the six fixed diagnostic targets. A/B JSON has
# two deliberate writes (pre-render and final PDF identity); every other target
# has one. Aggregate across all child Python processes to detect repeated low
# runs or observers overwriting an earlier result. Native ROOT writes are not
# covered by this hook; those use byte-preserving private input/model copies.
OBSERVER_AUDIT = r'''
import atexit, hashlib, json, os, pathlib, sys, uuid
manifest_path = os.environ.get("LT_ANALYSIS_PION_MATCHED_COHORT_MANIFEST")
if manifest_path:
    raw = pathlib.Path(manifest_path).read_bytes()
    if hashlib.sha256(raw).hexdigest() != os.environ.get("LT_ANALYSIS_PION_MATCHED_COHORT_SHA256"):
        raise RuntimeError("observer_audit_manifest_hash")
    manifest = json.loads(raw)
    output = pathlib.Path(manifest["private_outpath"])
    limits = {}
    for setting in ("Center", "Left"):
        base = setting + "_kaon_Q4p4W2p74_lowe_"
        limits[base + "pion-t-template-comparison.json"] = 2
        limits[base + "pion-t-template-comparison.pdf"] = 1
        limits[base + "pion-local-t-events.json"] = 1
    counts = {}
    def audit(event, args):
        if event != "open" or not isinstance(args[0], (str, bytes, os.PathLike)):
            return
        mode, flags = args[1:3]
        writing = (isinstance(mode, str) and any(c in mode for c in "wax+")) or (isinstance(flags, int) and flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND))
        if not writing:
            return
        p = pathlib.Path(os.fsdecode(args[0]))
        if p.name in limits:
            if p.resolve().parent != output:
                raise RuntimeError("observer_output_escape")
            counts[p.name] = counts.get(p.name, 0) + 1
            if counts[p.name] > limits[p.name]:
                raise RuntimeError("observer_output_overwrite:" + p.name)
    def save():
        path = output / ("observer-open-audit-" + uuid.uuid4().hex + ".json")
        with path.open("x", encoding="utf-8") as stream:
            json.dump({"manifest_sha256": hashlib.sha256(raw).hexdigest(), "counts": counts}, stream)
    sys.addaudithook(audit)
    atexit.register(save)
'''


def packet_name(setting):
    return setting + "_kaon_Q4p4W2p74_lowe_pion-local-t-events.json"


def comparison_name(setting):
    return setting + "_kaon_Q4p4W2p74_lowe_pion-t-template-comparison.json"


def artifacts():
    return [name for s in SETTINGS for name in
            (comparison_name(s), comparison_name(s).replace('.json', '.pdf'), packet_name(s))]


def identity(path):
    path = Path(path)
    info = path.stat()
    require(stat.S_ISREG(info.st_mode), "nonregular_input:" + str(path))
    result = capture.file_identity(path)
    result.update(device=info.st_dev, inode=info.st_ino, mode=stat.S_IMODE(info.st_mode),
                  mtime_ns=info.st_mtime_ns, resolved_path=str(path.resolve(strict=True)))
    return result


def inventory(root):
    """Do not follow directory symlinks; record every referent actually hashed."""
    root = Path(root)
    require(root.is_dir(), "inventory_root_missing:" + str(root))
    result = {}
    for parent, dirs, files in os.walk(root, followlinks=False):
        for name in sorted(dirs + files):
            p = Path(parent) / name
            key = p.relative_to(root).as_posix()
            if p.is_symlink():
                result[key] = {"link": os.readlink(p), "resolved_path": str(p.resolve(strict=True))}
                if p.is_file():
                    result[key]["referent"] = identity(p)
            elif p.is_file():
                result[key] = identity(p)
            else:
                result[key] = {"directory": True, "mode": stat.S_IMODE(p.stat().st_mode)}
    return result


def prepare_private_output(worktree, paths):
    """Only owner-created siblings; no replacement of existing links or trees."""
    worktree = Path(worktree)
    require(worktree.is_absolute() and worktree.resolve() == worktree and
            worktree.name == 'analysis' and not worktree.is_symlink(), "private_worktree_invalid")
    require(Path(paths['LTANAPATH']).resolve() == worktree and paths['ANATYPE'] == 'Kaon',
            "private_ltanapath_or_anatype")
    lexical = worktree / 'OUTPUT/Analysis/KaonLT'
    require(Path(paths['OUTPATH']) == lexical, "private_lexical_outpath")
    private = worktree.parent / 'matched-output'
    link = worktree / 'OUTPUT'
    require(not os.path.lexists(private) and not os.path.lexists(link), "private_output_collision")
    if hasattr(os, 'geteuid'):
        require(worktree.parent.stat().st_uid == os.geteuid() and
                stat.S_IMODE(worktree.parent.stat().st_mode) & 0o077 == 0, 'private_parent_permissions')
    private.mkdir(mode=0o700)
    output = private / 'Analysis/KaonLT'
    output.mkdir(parents=True, mode=0o700)
    link.symlink_to(private, target_is_directory=True)
    require(link.resolve() == private and lexical.resolve() == output, "private_output_resolution")
    shared = (Path(paths['VOLATILEPATH']) / 'OUTPUT').resolve()
    require(not private.is_relative_to(shared) and not shared.is_relative_to(private), 'private_output_overlaps_shared_output')
    return output


def probe_private_callers(repo, worktree, env, baseline, outdir):
    """No scientific imports: interrogate actual no-mode and Plot Root calls."""
    probes = isolation.probe_paths(worktree, env, baseline)
    records = []
    callers = [(role, caller, modes, cwd) for role, caller, modes, cwd, _ in resolver.CALLERS]
    known = {(caller, tuple(modes)) for _, caller, modes, _ in callers}
    callers += [(caller, caller, modes, cwd) for caller, modes, cwd in reached_resolver_callers(repo)
                if (caller, tuple(modes)) not in known]
    for role, caller, modes, cwd in callers:
        record = resolver.probe(worktree, caller, modes, env, cwd)
        base = resolver.probe(repo, caller, modes, os.environ, cwd)
        fields = record['paths']
        for key in resolver.FIELDS:
            if key not in {'LTANAPATH', 'OUTPATH'}:
                require(fields.get(key) == base['paths'].get(key), 'stable_input_role_changed:' + role + ':' + key)
        raw = str(fields['OUTPATH']).replace('None', fields['ANATYPE'] + 'LT') if not modes else fields['OUTPATH']
        require(Path(raw) == Path(worktree) / 'OUTPUT/Analysis/KaonLT' and
                Path(raw).resolve() == outdir and Path(raw).is_absolute(),
                'default_root_private_outpath_mismatch:' + role)
        require(resolver.module_identity(record) == worktree.parent / 'runtime-python/ltsep',
                'copied_ltsep_import_mismatch:' + role)
        records.append(record)
    return {'plot': probes, 'all_callers': records}


def reached_resolver_callers(repo):
    """Source-only import closure from the unchanged child entry points.

    Resolve unique repository modules, including function-local imports. Never
    import production Python or execute its ROOT/Fortran side effects. A local
    import ambiguity or a new nonliteral Root mode stops this narrow audit.
    """
    repo = Path(repo)
    modules = {}
    for path in (repo / 'src').rglob('*.py'):
        modules.setdefault(path.stem, []).append(path)
    pending = [(repo / name, cwd) for name, cwd in (
        ('src/main.py', 'src'), ('src/setup/findEffectiveCharge.py', 'src/setup'),
        ('src/setup/calcTotalEffectiveCharge.py', 'src/setup'), ('src/setup/set_sig_fortran.py', 'src/setup'))]
    seen, callers = set(), []
    while pending:
        path, cwd = pending.pop()
        if (path, cwd) in seen:
            continue
        seen.add((path, cwd))
        tree = ast.parse(path.read_text(encoding='utf-8'))
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id in {'Root', 'create_ltsep_root'}:
                require(node.args and not node.keywords, 'resolver_call_layout:' + str(path))
                modes = []
                for mode in node.args[1:]:
                    require(isinstance(mode, ast.Constant) and isinstance(mode.value, str), 'nonliteral_resolver_mode:' + str(path))
                    modes.append(mode.value)
                require(modes in ([], ['Plot_LTSep']), 'unapproved_resolver_mode:' + str(path))
                callers.append((path.relative_to(repo).as_posix(), tuple(modes), cwd))
            names = ([node.module] if isinstance(node, ast.ImportFrom) and node.module else
                     [alias.name for alias in node.names] if isinstance(node, ast.Import) else [])
            for name in names:
                candidates = modules.get(name.split('.')[-1], [])
                require(len(candidates) <= 1, 'ambiguous_source_import:' + name)
                if candidates:
                    pending.append((candidates[0], cwd))
    return sorted(callers)


def stage_input_files(shared, outdir):
    """Source-derived physical inputs used by low/high shared preflight.

    Private copies are intentional: ordinary analysed ROOTs can be opened by
    legacy helpers and must never provide a write-through path to originals.
    No historical diagnostics, weights, candidate models or fits are staged.
    """
    shared, outdir = Path(shared), Path(outdir)
    records = {}
    for eps, settings in (('lowe', ('Center', 'Left', 'Right')), ('highe', ('Center', 'Left', 'Right'))):
        for setting in settings:
            names = [f'{setting}_kaon_Analysed_{role}_Q4p4W2p74_{eps}.root' for role in ('Data', 'Dummy')]
            names += [f'Prod_Coin_Q4p4W2p74{setting.lower()}_{eps}{suffix}'
                      for suffix in (('.root', '.hist') if eps == 'lowe' else ('.root',))]
            # Preserve source availability: a missing setting is not invented.
            present = [(shared / name).is_file() for name in names[:2]]
            if not any(present):
                require(setting not in SETTINGS or eps != 'lowe', 'required_lowe_setting_missing:' + setting)
                continue
            require(all(present), 'incomplete_data_dummy_pair:' + setting + ':' + eps)
            for name in names:
                source, target = shared / name, outdir / name
                if name.startswith('Prod_Coin_') and not source.is_file() and not (eps == 'lowe' and setting in SETTINGS):
                    continue  # Preserve the ordinary shift_prep availability decision.
                require(not os.path.lexists(target), 'staged_input_collision:' + name)
                before = identity(source)
                with source.open('rb') as src, target.open('xb') as dst:
                    shutil.copyfileobj(src, dst, 1024 * 1024)
                require(identity(source) == before, 'input_changed_during_copy:' + name)
                require(capture.file_identity(target) == {k: before[k] for k in ('sha256', 'bytes')},
                        'private_copy_identity:' + name)
                require(target.stat().st_ino != source.stat().st_ino or target.stat().st_dev != source.stat().st_dev,
                        'private_copy_hardlink:' + name)
                records[name] = {'source': str(source), 'original': before, 'private': str(target),
                                 'private_mutation': 'frozen_shiftprep_derived_branches' if '_Analysed_' in name else 'none'}
    # DiamondPlot reads shortest matching filenames, including mid-epsilon and
    # center fallback candidates. Preserve its actual candidate namespace; do
    # not silently change the common cut by staging just four lowe ROOT files.
    for source in map(Path, glob.glob(str(shared / '*kaon*Q4p4W2p74*.root'))):
        if 'simc' in str(source).lower() or 'shiftprep' in source.name.lower() or source.name in records:
            continue
        target = outdir / source.name
        require(not os.path.lexists(target), 'diamond_input_collision:' + source.name)
        before = identity(source)
        with source.open('rb') as src, target.open('xb') as dst:
            shutil.copyfileobj(src, dst, 1024 * 1024)
        require(identity(source) == before and capture.file_identity(target) ==
                {k: before[k] for k in ('sha256', 'bytes')}, 'diamond_copy_identity:' + source.name)
        records[source.name] = {'source': str(source), 'original': before, 'private': str(target),
                                'private_mutation': 'ordinary_private_output_replacement'}
    require(diamond_selections(shared) == diamond_selections(outdir), 'diamond_input_selection_changed')
    return records


def diamond_selections(root):
    """Source-derived shortest-file choices; preserve case and first-tie rules."""
    root = Path(root)
    def choose(tokens, eps, fallback=False, main=False):
        best = None
        for tok in tokens:
            for path in glob.glob(str(root / ('*' + tok + '*kaon*Q4p4W2p74*.root'))):
                low = path.lower()
                if 'simc' in low or 'shiftprep' in Path(path).name.lower() or eps not in (path if main else low):
                    continue
                if fallback and ('left' in low or 'right' in low):
                    continue
                if best is None or len(path) < len(best):
                    best = path
        return None if best is None else Path(best).name
    return {f'{kind}:{s}:{eps}': choose(tokens, eps, fallback=kind == 'fallback', main=kind == 'main')
            for kind, settings in (('main', [(s, [s]) for s in ('Center', 'Left', 'Right')]),
                                   ('common_overlay', [(s, [s, s.lower()]) for s in ('Center', 'Left', 'Right')]),
                                   ('fallback', [('Center', [''])]))
            for s, tokens in settings for eps in ('low', 'mid', 'high')}


def runtime_input_identities(outdir):
    return {s: {role: capture.file_identity(Path(outdir) / f'{s}_kaon_Analysed_{title}_Q4p4W2p74_lowe.root')
                for role, title in (('data', 'Data'), ('dummy', 'Dummy'))} for s in SETTINGS}


def skim_paths(repo, paths):
    """Evaluate only the two existing pure path normalizers, not setup scripts."""
    source = (Path(repo) / 'src/setup/findEffectiveCharge.py').read_text(encoding='utf-8')
    tree = ast.parse(source)
    function = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == 'resolve_skim_outpath']
    require(len(function) == 1, 'effective_skim_normalizer_missing')
    namespace = {'os': os}
    exec(compile(ast.Module(body=function, type_ignores=[]), '<frozen-skim-normalizer>', 'exec'), namespace)
    effective = namespace['resolve_skim_outpath'](paths['SKIMPATH'], paths['ANATYPE'])
    launcher = (Path(repo) / 'run_Prod_Analysis.sh').read_text(encoding='utf-8')
    require(launcher.count('normalize_ltsep_dir() {') == 1, 'launcher_skim_normalizer_missing')
    body = 'normalize_ltsep_dir() {' + launcher.split('normalize_ltsep_dir() {', 1)[1].split('\n}', 1)[0] + '\n}\n'
    result = subprocess.run(['bash', '-c', body + '\nANATYPE="$1"\nnormalize_ltsep_dir "$2"',
                             'matched-skim-path', paths['ANATYPE'], paths['SKIMPATH']], capture_output=True, text=True)
    require(result.returncode == 0 and result.stdout.strip() and not result.stderr, 'launcher_skim_normalizer_failed')
    return {'launcher': result.stdout.strip(), 'effective_charge': effective}


def external_input_inventory(worktree, paths, private, env):
    """Declared stable read roles for the exact unchanged Q4p4W2p74 -d child.

    Background ROOT consumers use TFile.Open(..., READ); skim paths are tested
    for existence by check_runs_in_effcharge, never merged under plain -d.
    PID logs, efficiency tables and timing parameters are read-only consumers.
    SIMC external symlink repair is independently prohibited before this call.
    """
    records = {}
    def add(path, role, required=True):
        path = Path(path)
        if not path.is_file():
            require(not required, 'required_external_input_missing:' + role + ':' + str(path))
            records[str(path)] = {'role': role, 'missing': True}
        else:
            records[str(path)] = {'role': role, 'identity': identity(path)}
    launcher = (worktree / 'run_Prod_Analysis.sh').read_text(encoding='utf-8')
    tables = re.findall(r'(?m)^EffData="([^"\r\n]+)"', launcher)
    require(len(tables) == 1, 'efficiency_table_source_not_unique')
    add(Path(paths['UTILPATH']) / 'scripts/efficiency/OUTPUTS' / tables[0], 'efficiency_table')
    add(Path(paths['PARAMPATH']) / 'Timing_Parameters.csv', 'timing_parameters')
    shared = Path(paths['VOLATILEPATH']) / 'OUTPUT/Analysis/KaonLT'
    skim = skim_paths(worktree, paths)['effective_charge']
    for eps, settings in (('low', ('Center', 'Left', 'Right')), ('high', ('Center', 'Left', 'Right'))):
        for s in settings:
            if not (private / f'{s}_kaon_Analysed_Data_Q4p4W2p74_{eps}e.root').is_file():
                continue
            for dummy in ('', '_dummy'):
                runlist = Path(paths['UTILPATH']) / 'run_list/KaonLT' / f'Q4p4W2p74{s.lower()}_{eps}e{dummy}'
                add(runlist, 'run_list')
                runs = re.split(r'\s+', runlist.read_text(encoding='utf-8').strip())
                require(all(re.fullmatch(r'\d+', run) for run in runs), 'run_list_not_integer')
                # Skim existence controls effective-charge selection. Preserve
                # every observed input and every absent-input decision.
                for run in runs:
                    add(Path(skim) / f'kaon_{run}_-1_Raw_Data.root', 'skim_existence', required=False)
                    add(Path(paths['VOLATILEPATH']) / 'log' / f'{s}_kaon_{run}_Q4p4W2p74_{eps}e.log',
                        'PID_log', required=False)
            # Full low analysis uses these read-only external ROOT templates;
            # high canonical raw prepass exits before model preparation.
            if eps == 'low' and (private / f'Prod_Coin_Q4p4W2p74{s.lower()}_lowe.root').is_file():
                for component in ('neutron', 'delta', 'sidis'):
                    add(shared / 'background_samples' / component /
                        f'Prod_Coin_Q4p4W2p74{s.lower()}_{eps}e.root', 'pion_SIMC_' + component)
                # The unchanged launcher permits unconfigured Sigma0; its
                # runtime manifest strips blanks and retains the existing
                # unavailable-template fallback. Verify every configured file.
                sigma = str(env.get(f'LT_BG_SIGMA0_{s.upper()}_{eps.upper()}_ROOT', '')).strip()
                if sigma:
                    add(sigma, 'sigma0_SIMC')
    # External cut/config inputs are stable resolver roles, not new guessed
    # paths. Record their complete regular-file inventory without link repair.
    for role in ('CUTPATH', 'PARAMPATH'):
        root = Path(paths[role])
        if root.is_dir():
            for name in inventory(root):
                p = root / name
                if p.is_file():
                    add(p, role)
    return records


def verify_external_inputs(records):
    for path, record in records.items():
        if record.get('missing'):
            require(not Path(path).exists(), 'external_missing_input_changed:' + path)
        else:
            require(identity(path) == record['identity'], 'external_input_changed:' + path)


def verify_external_links(record):
    for item in record['links']:
        link = Path(item['path'])
        require(link.is_symlink() and link.exists() and os.readlink(link) == item['literal_target']
                and str(link.resolve()) == item['resolved_target'], 'external_link_changed:' + str(link))


def preservation_checks(checks):
    """A failed check must not suppress independent preservation checks."""
    failures = []
    for name, check in checks:
        try:
            check()
        except Exception as exc:
            failures.append(name + ':' + str(exc))
    require(not failures, 'preservation_failed:' + ';'.join(failures))


def toolchain():
    python = shutil.which('python3')
    pdfinfo = shutil.which('pdfinfo')
    require(python and Path(python).resolve() == Path(sys.executable).resolve(), 'child_python_interpreter_mismatch')
    require(pdfinfo, 'pdf_structure_backend_unavailable')
    return {'python3': str(Path(python).resolve()), 'pdfinfo': str(Path(pdfinfo).resolve())}


def run_owned_child(worktree, log, env):
    copied = Path(worktree).parent / 'runtime-python/ltsep'
    require(copied.is_dir() and not copied.is_symlink(), 'copied_ltsep_directory_missing')
    before = resolver.tree_bytes(copied)
    try:
        return ab.debug.run_analysis(worktree, log, env)
    finally:
        require(resolver.tree_bytes(copied) == before, 'copied_ltsep_changed_during_child')


def child_environment(env, manifest_path=None):
    result = dict(env)
    for key in ('LT_ANALYSIS_DEBUG_LEFT_LOW', 'LT_ANALYSIS_CANONICAL_PREPASS_CAPTURE',
                'LT_ANALYSIS_ALLOW_UNPAIRED_CANONICAL_BINNING'):
        result.pop(key, None)
    require(not any(k.startswith('LT_ANALYSIS_') and v for k, v in result.items()),
            'untrusted_analysis_environment')
    if manifest_path is None:
        return result
    result[ab.ENABLE_ENV] = result[capture.ENABLE_ENV] = '1'
    result[AUTH_ENV] = str(manifest_path)
    result[capture.MATCHED_HASH_ENV] = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    result['PYTHONDONTWRITEBYTECODE'] = '1'
    return result


def verify_joint(outdir, manifest, source_identities):
    """Accept neither a successful child nor file existence without joins."""
    outdir = Path(outdir)
    if manifest.get('private_input_policy') == 'frozen_shiftprep_private_copy':
        require(str(outdir.resolve()) == manifest['private_outpath'], 'joint_private_output_path')
    counts = {name: 0 for name in artifacts()}
    for audit in outdir.glob('observer-open-audit-*.json'):
        require(not audit.is_symlink(), 'observer_audit_symlink')
        record = ab.strict_json(audit.read_bytes())
        require(record['manifest_sha256'] == capture.digest(manifest), 'observer_audit_wrong_attempt')
        for name, count in record['counts'].items():
            require(name in counts and type(count) is int and count > 0, 'observer_audit_invalid')
            counts[name] += count
    require(counts == {name: 2 if name.endswith('comparison.json') else 1 for name in artifacts()},
            'observer_output_missing_duplicate_or_overwritten')
    records, docs, packets = {}, {}, {}
    for name in artifacts():
        p = outdir / name
        require(p.is_file() and not p.is_symlink() and p.stat().st_mtime_ns >= manifest['started_ns'],
                'joint_artifact_missing_or_stale:' + name)
        require(p.stat().st_size > 0, 'joint_artifact_empty:' + name)
        if p.suffix == '.pdf':
            with p.open('rb') as stream:
                require(stream.read(5) == b'%PDF-', 'joint_PDF_signature')
        records[name] = capture.file_identity(p)
    observed_runtime_inputs = (runtime_input_identities(outdir) if manifest.get('private_input_policy') ==
                               'frozen_shiftprep_private_copy' else manifest['inputs'])
    for s in SETTINGS:
        doc = ab.strict_json((outdir / comparison_name(s)).read_bytes())
        pdfname = comparison_name(s).replace('.json', '.pdf')
        require(ab.pdf_pages(outdir / pdfname) == 6, 'joint_AB_pdf_pages')
        ab.verify_payload(doc, s, manifest['source_HEAD'], source_identities,
                          dict(records[pdfname], basename=pdfname))
        require(doc['status'] == 'available', 'joint_AB_unavailable:' + s)
        if manifest.get('private_input_policy') == 'frozen_shiftprep_private_copy':
            require(Path(doc['pdf']['path']).resolve() == (outdir / pdfname).resolve(), 'joint_PDF_path')
        packet = ab.strict_json((outdir / packet_name(s)).read_bytes())
        capture.validate_packet(packet)
        require(packet['source_HEAD'] == manifest['source_HEAD'] and packet['setting'] == s,
                'joint_packet_source_or_setting')
        expected_inputs = observed_runtime_inputs[s]
        require(packet['input_identity'] == expected_inputs, 'joint_packet_input_identity')
        require(packet['source_files'] == manifest['source_files'], 'joint_packet_source_files')
        authority = packet.get('cohort_authority', {})
        require(authority.get('mode') == 'matched_cohort' and authority.get('attempt_id') == manifest['attempt_id']
                and authority.get('manifest_sha256') == capture.digest(manifest)
                and authority.get('AB') == records[comparison_name(s)], 'joint_packet_authority')
        if manifest.get('private_input_policy') == 'frozen_shiftprep_private_copy':
            require(authority.get('original_input_identity') == manifest['inputs'][s] and
                    authority.get('runtime_input_identity') == expected_inputs, 'joint_original_runtime_input_join')
            for role, title in (('data', 'Data'), ('dummy', 'Dummy')):
                expected_path = outdir / f'{s}_kaon_Analysed_{title}_Q4p4W2p74_lowe.root'
                require(Path(packet['source_input_paths'][role]).resolve() == expected_path.resolve(),
                        'joint_backing_file_path')
                require(Path(doc['original_analysis_input_identity'][role]['path']).resolve() == expected_path.resolve(),
                        'joint_AB_backing_file_path')
        require(packet['t_edges'] == doc['canonical_t_edges'] and packet['phi_edges'] == doc['canonical_phi_edges'],
                'joint_coordinate_edges')
        for role in ('data', 'dummy'):
            require(all(doc['original_analysis_input_identity'][role][k] == expected_inputs[role][k]
                        for k in ('sha256', 'bytes')), 'joint_AB_input_identity')
        for p, parent in zip(packet['parents'], doc['parents']):
            for key in ('t_bin_index', 't_edges', 'pion_parent_id', 'canonical_interval_pair_id',
                        'canonical_interval_pair_hash', 'coordinate_fingerprint', 'baseline_parent_fingerprint'):
                require(p[key] == parent[key], 'joint_parent_identity:' + key)
            require(packet['coordinate_fingerprint'] == parent['coordinate_fingerprint'], 'joint_coordinate')
            for side, key in (('pion', 'pion_input'), ('kaon', 'kaon_input')):
                require(p[key] == parent['frozen_inputs'][side], 'joint_histogram_identity')
            # Reconstruct the unweighted input spectra independently of C fits.
            try:
                from . import pion_local_t_event_validation_replay as replay
            except ImportError:
                import pion_local_t_event_validation_replay as replay
            for pop, key in (('pion', 'pion_input'), ('kaon', 'kaon_input')):
                rows = [r for r in packet['events'][pop] if r['t_index'] == p['t_bin_index'] and r['nommcuts']]
                spectrum = replay.split_spectra(rows, p[key]['edges'])['parent']
                for observed, saved in (('signed', 'values'), ('sumw2', 'stored_sumw2')):
                    import math
                    require(all(math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)
                                for a, b in zip(spectrum[observed], p[key][saved])),
                            'joint_source_histogram_closure:' + s + ':' + pop)
        docs[s], packets[s] = doc, packet
    require(docs['Center']['canonical_t_edges'] == docs['Left']['canonical_t_edges'], 'joint_mixed_canonical_t')
    return records, docs, packets


def profile(source_commit):
    return {'schema_version': 'pion_hgcer_validation_bundle_profile/v4', 'validation_profile': PROFILE_ID,
            'collection_mode': 'generic_artifacts', 'settings': deepcopy(ab.SETTINGS),
            'artifacts': {'global': [{'key': 'artifact_' + str(i), 'basename_template': name,
                                    'kind': 'json' if name.endswith('.json') else 'file', 'required': True}
                                   for i, name in enumerate(artifacts() + ['matched-cohort-authority.json'])],
                          'settings': []},
            'source_identity': {'required_analysis_commit': source_commit, 'allowed_committed_files': [],
                                'allowed_non_analysis_path_prefixes': ['docs/memory/']}}


def verify_archive(path, expected_head=None):
    """Read ZIP members directly, never extract or trust member paths."""
    with zipfile.ZipFile(path) as z:
        require(z.testzip() is None and len(z.namelist()) == len(set(z.namelist())), 'matched_zip_crc_or_duplicates')
        m = ab.strict_json(z.read('manifest.json'))
        head = m.get('required_analysis_commit')
        require(re.fullmatch('[0-9a-f]{40}', head or '') is not None and
                m.get('git_head') == head and (expected_head is None or head == expected_head)
                and m.get('complete') is True and not m.get('errors')
                and m.get('validation_profile') == PROFILE_ID, 'matched_zip_manifest')
        require(m.get('requested_kinematic') == 'Q4p4W2p74' and m.get('requested_settings') == ab.SETTINGS,
                'matched_zip_setting_scope')
        require(m.get('required_analysis_commit_is_ancestor') is True and m.get('allowed_committed_files') == []
                and m.get('committed_files_after_required_analysis_commit') == []
                and m.get('unexpected_committed_files_after_required_analysis_commit') == [], 'matched_zip_source_identity')
        checks = m.get('source_checks', [])
        require(isinstance(checks, list) and checks and all(c.get('returncode') == 0 for c in checks),
                'matched_zip_source_checks')
        require(m.get('settings') == [dict(s, kinematic='Q4p4W2p74', artifacts={}) for s in ab.SETTINGS],
                'matched_zip_metadata_settings')
        declarations = profile(head)['artifacts']['global']
        expected = {'global/' + row['basename_template'] for row in declarations}
        expected_members = expected | {'manifest.json', 'source_state.txt', 'source_checks.txt', 'global/'}
        expected_members.update(s['phi'] + '_' + s['epsilon'] + '/' for s in ab.SETTINGS)
        require(set(z.namelist()) == expected_members, 'matched_zip_member_inventory')
        require({x for x in z.namelist() if x.startswith('global/') and not x.endswith('/')} == expected,
                'matched_zip_artifact_inventory')
        saved = m.get('global_artifacts', {})
        require(set(saved) == {row['key'] for row in declarations}, 'matched_zip_manifest_inventory')
        raw = {}
        for row in declarations:
            name = row['basename_template']
            data = z.read('global/' + name)
            entry = saved[row['key']]
            require(entry.get('status') == 'exists' and entry.get('archive_path') == 'global/' + name and
                    entry.get('sha256') == hashlib.sha256(data).hexdigest() and entry.get('byte_size') == len(data),
                    'matched_zip_hash:' + name)
            if row['kind'] == 'json':
                require(entry.get('json_status') == 'valid', 'matched_zip_JSON_status:' + name)
            raw[name] = data
        authority = ab.strict_json(raw['matched-cohort-authority.json'])
        require(authority['source_HEAD'] == head and authority.get('schema_version') == 'pion_matched_cohort_authority/v1'
                and authority.get('mode') == 'matched_cohort' and set(authority.get('inputs', {})) == set(SETTINGS),
                'matched_zip_authority_source')
        return raw, authority, {'archive': capture.file_identity(path), 'manifest_sha256': hashlib.sha256(z.read('manifest.json')).hexdigest()}


def execute_gate(repo, outdir, output, source_commit):
    repo, outdir, output = map(lambda p: Path(p).absolute(), (repo, outdir, output))
    require(repo.resolve() == isolation.accepted.REPO.resolve() and
            outdir.resolve() == isolation.accepted.ARTIFACTS.resolve() and
            output.parent.resolve() == isolation.accepted.GLOBUS.resolve() and output.suffix == '.zip',
            'owner_path_boundary_invalid')
    before = isolation.snapshot(repo)
    require(before['head'] == before['origin_test'] == source_commit and before['branch'] == 'test',
            'owner_source_not_published_test')
    require(source_commit != BASE, 'owner_not_committed')
    ab.git(repo, 'merge-base', '--is-ancestor', BASE, source_commit)
    allowed = {OWNER, 'testing/test_run_pion_local_t_matched_cohort_gate.py',
               'testing/pion_local_t_event_validation_capture.py', 'testing/test_pion_local_t_event_validation.py',
               'testing/pion_local_t_transfer_variance_comparison.py', 'testing/pion_local_t_event_validation_replay.py',
               'testing/test_pion_local_t_transfer_variance_comparison.py', 'testing/test_pion_local_t_matched_cohort_science.py',
               'src/binning/calculate_yield.py', 'src/cuts/rand_sub.py'}
    require(all(p in allowed or p.startswith('docs/memory/') for p in
                ab.git(repo, 'diff', '--name-only', BASE, source_commit).splitlines()), 'unexpected_committed_source')
    for line in before['porcelain'].splitlines():
        require(line[3:] in isolation.FARM_OUTPUTS or line[3:].startswith('OUTPUT/'), 'dirty_gate_source:' + line)
    destinations = ab.destinations(outdir, output)
    require(all(not os.path.lexists(p) for p in destinations.values()), 'owner_destination_collision')
    child_environment(os.environ)  # Reject inherited opt-ins before any mutation.
    status = ab.GateStatus(outdir, output, source_commit)
    status.update(schema_version='pion_local_t_matched_cohort_gate_status/v1',
                  boundaries={'diagnostic_only': True, 'scientifically_accepted': False,
                              'production_applied': False, 'private_child_outputs_required': True,
                              'shared_inputs_outputs_preserved': 'NOT VERIFIED',
                              'pdf_legibility_verified': False, 'method_a_detached': True,
                              'method_b_numerically_excluded': True})
    overlay = installed = None
    shared_before = None
    try:
        tools = toolchain()
        audit = audit_default_output_writer(repo)
        status.update('source_preflight', writer_audit=audit, toolchain=tools)
        with isolation.owned_worktree(repo, source_commit, 'source') as source:
            with isolation.collection_module(source) as module, redirect_stdout(sys.stderr):
                _, checks = module.collect_source_checks(source, required_analysis_commit=source_commit,
                                                         allowed_committed_files=[])
                require(all(c['returncode'] == 0 for c in checks), 'collector_source_check_failed')
        original_import = isolation.import_identity(repo)
        package = Path(original_import['package_file']).resolve().parent
        installed_package_before = resolver.tree_bytes(package)
        installed = {'original_source_sha256': {str(p): ab.sha256(p) for p in
                     (Path(original_import['package_file']), Path(original_import['pathing_file']),
                      *sorted((package / 'PATH_TO_DIR').glob('*.path')))}}
        with isolation.owned_worktree(repo, source_commit, 'analysis') as worktree:
            env, overlay = isolation.prepare_runtime_overlay(repo, worktree)
            paths = overlay['baseline_paths']
            shared = Path(paths['VOLATILEPATH']) / 'OUTPUT/Analysis' / (paths['ANATYPE'] + 'LT')
            require(shared.resolve() == outdir.resolve() == Path(paths['OUTPATH']).resolve(),
                    'ordinary_output_authority_mismatch')
            # All shared outputs, not just historical A/B filenames. Exclude only
            # this owner's explicitly named log/status/summary destinations.
            excluded = {p.name for p in destinations.values() if p.parent == outdir}
            def shared_state():
                return {k: v for k, v in inventory(shared).items() if k not in excluded}
            shared_before = shared_state()
            external = isolation.external_symlink_preflight(worktree, paths)
            if external['background_simc_path']:
                require(Path(external['background_simc_path']).is_dir(), 'background_simc_directory_missing')
            require((shared / 'background_samples').is_dir(), 'external_background_mkdir_required')
            # Launcher mkdir of the stable skim role must already be a no-op.
            skim = skim_paths(worktree, paths)['launcher']
            require(Path(skim).is_dir(), 'external_skim_mkdir_required')
            output_paths = isolation.probe_paths(worktree, env, paths)[0]['paths']
            private = prepare_private_output(worktree, output_paths)
            probes = probe_private_callers(repo, worktree, env, paths, private)
            bindings = stage_input_files(shared, private)
            external_inputs = external_input_inventory(worktree, paths, private, env)
            audit_path = worktree.parent / 'runtime-python/sitecustomize.py'
            require(not os.path.lexists(audit_path), 'observer_audit_code_collision')
            with audit_path.open('x', encoding='utf-8') as stream:
                stream.write(OBSERVER_AUDIT)
            inputs = {s: {role: capture.file_identity(private / f'{s}_kaon_Analysed_{title}_Q4p4W2p74_lowe.root')
                          for role, title in (('data', 'Data'), ('dummy', 'Dummy'))} for s in SETTINGS}
            started = time.time_ns()
            authority = {'schema_version': 'pion_matched_cohort_authority/v1', 'mode': 'matched_cohort',
                         'attempt_id': str(uuid.uuid4()), 'source_HEAD': source_commit, 'started_ns': started,
                         'private_outpath': str(private), 'inputs': inputs, 'owner_source': capture.file_identity(worktree / OWNER),
                         'input_bindings': bindings, 'source_files': event_owner.committed_identities(repo, source_commit),
                         'external_inputs': external_inputs, 'observer_audit_sha256': ab.sha256(audit_path)}
            authority['private_input_policy'] = 'frozen_shiftprep_private_copy'
            authority['input_identity_roles'] = {
                'inputs': 'original shared ROOT container SHA256/bytes before private copying; never rewritten',
                'AB.original_analysis_input_identity': 'private runtime ROOT container after frozen derived-branch writes',
                'packet.input_identity': 'must equal same-run A/B private runtime container; no historical pin substitution'}
            authority['private_writer_source_files'] = {p: capture.file_identity(worktree / p) for p in
                ('src/cuts/shift_prep.py', 'src/setup/shift_MM.py')}
            auth_path = private / 'matched-cohort-authority.json'
            with auth_path.open('xb') as stream:
                stream.write(capture.encoded(authority))
            env = child_environment(env, auth_path)
            # The unchanged launcher uses TMPDIR for its canonical NPZ prepass;
            # plotting libraries also need private writable cache directories.
            for key, leaf in (('TMPDIR', 'runtime-tmp'), ('MPLCONFIGDIR', 'runtime-mpl'),
                              ('XDG_CACHE_HOME', 'runtime-cache')):
                target = worktree.parent / leaf
                require(not os.path.lexists(target), 'runtime_cache_collision:' + key)
                target.mkdir(mode=0o700)
                env[key] = str(target)
            # All gates above precede the one unchanged destructive/setup child.
            status.update('analysis', analysis_started=True, probes=probes, input_bindings=bindings,
                          external_symlink_preflight=external, authority_sha256=capture.digest(authority))
            code = run_owned_child(worktree, destinations['log'], env)
            status.update('completion_markers', analysis_returncode=code,
                          log_identity=capture.file_identity(destinations['log']))
            ab.verify_completion(destinations['log'], code)
            require('High Epsilon Completed!' not in destinations['log'].read_text(encoding='utf-8', errors='replace'),
                    'unexpected_full_high_analysis')
            status.update('verify_joint', analysis_completed=True)
            records, _, packets = verify_joint(private, authority, ab.committed_identities(repo, source_commit))
            runtime_inputs = {s: packets[s]['input_identity'] for s in SETTINGS}
            for rec in bindings.values():
                require(identity(rec['source']) == rec['original'], 'original_input_changed:' + rec['source'])
                if rec['private_mutation'] == 'none':
                    require(capture.file_identity(rec['private']) == {k: rec['original'][k] for k in ('sha256', 'bytes')},
                            'private_archived_SIMC_changed:' + rec['private'])
            verify_external_inputs(external_inputs)
            require(shared_state() == shared_before, 'shared_outputs_changed')
            verify_external_links(external)
            isolation.verify_ltsep_preservation(installed)
            isolation.verify_preservation(repo, before)
            status.update('collection', artifact_verification_completed=True)
            with tempfile.TemporaryDirectory(prefix='matched-profile-') as directory:
                pf = Path(directory) / 'profile.json'
                pf.write_bytes(capture.encoded(profile(source_commit)))
                with isolation.owned_worktree(repo, source_commit, 'source') as source:
                    with isolation.collection_module(source) as module, redirect_stdout(sys.stderr):
                        result = module.collect_validation_bundle(outdir=private, kinematic='Q4p4W2p74',
                            output=output, profile_path=pf, repo_root=source)
                    require(result['returncode'] == 0, 'collection_failed')
            verify_archive(output, source_commit)
            status.update(zip_verification_completed=True, collection_completed=True)
        isolation.verify_preservation(repo, before)

        require(shared_state() == shared_before, 'shared_outputs_changed_after_cleanup')
        verify_external_inputs(external_inputs)
        verify_external_links(external)
        require(resolver.tree_bytes(package) == installed_package_before, 'installed_ltsep_package_changed')
        isolation.verify_ltsep_preservation(installed)
        summary = {'schema_version': PROFILE_ID, 'source_HEAD': source_commit, 'artifacts': records,
                   'authority_sha256': capture.digest(authority), 'diagnostic_only': True,
                   'scientifically_accepted': False, 'zip_identity': capture.file_identity(output),
                   'log_identity': capture.file_identity(destinations['log']),
                   'original_input_identity': inputs, 'runtime_input_identity': runtime_inputs,
                   'private_input_policy': authority['private_input_policy']}
        with destinations['summary'].open('xb') as stream:
            stream.write(capture.encoded(summary))
    except Exception as exc:
        status.update(status='failed', failure_reason=str(exc), evidence_admissible=False)
        raise
    finally:
        checks = [('ordinary_checkout', lambda: isolation.verify_preservation(repo, before))]
        if shared_before is not None:
            checks.append(('shared_outputs', lambda: require(shared_state() == shared_before, 'shared_outputs_changed')))
        if installed is not None:
            checks.append(('installed_ltsep', lambda: isolation.verify_ltsep_preservation(installed)))
        if 'installed_package_before' in locals():
            checks.append(('installed_package', lambda: require(resolver.tree_bytes(package) == installed_package_before,
                                                                 'installed_ltsep_package_changed')))
        if 'external_inputs' in locals():
            checks.append(('external_inputs', lambda: verify_external_inputs(external_inputs)))
        if 'external' in locals():
            checks.append(('external_links', lambda: verify_external_links(external)))
        try:
            preservation_checks(checks)
        except Exception as exc:
            status.update(status='failed', evidence_admissible=False, preservation_failure_reason=str(exc))
            raise
    # Delivery is reachable only after the finally preservation barrier passes.
    try:
        status.update(status='success', artifacts=records, scientifically_accepted=False,
                      shared_inputs_outputs_preserved=True, worktree_cleanup_completed=True)
        isolation.deliver_evidence(status, destinations['log'], destinations['summary'], output)
    except Exception as exc:
        status.update(status='failed', failure_reason=str(exc), evidence_admissible=False)
        raise
    return output


def offline_science(bundle, output_dir, expected_head):
    """One explicit new-cohort route; original archive/JSON CLIs stay strict."""
    try:
        from . import pion_local_t_transfer_variance_comparison as c
        from . import pion_local_t_event_validation_replay as replay
    except ImportError:
        import pion_local_t_transfer_variance_comparison as c
        import pion_local_t_event_validation_replay as replay
    output_dir = Path(output_dir).absolute()
    require(output_dir.is_dir() and not output_dir.is_symlink(), 'offline_output_directory')
    repo = Path(__file__).resolve().parents[1]
    require(not output_dir.resolve().is_relative_to(repo), 'offline_output_inside_checkout')
    targets = {key: output_dir / name for key, name in
               (('candidate', 'matched-cohort-C.json'), ('replay', 'matched-cohort-science.json'),
                ('pdf', 'matched-cohort-science.pdf'))}
    require(all(not os.path.lexists(p) for p in targets.values()), 'offline_output_collision')
    raw, authority, provenance = verify_archive(bundle, expected_head)
    docs = {s: ab.strict_json(raw[comparison_name(s)]) for s in SETTINGS}
    packets = {s: ab.strict_json(raw[packet_name(s)]) for s in SETTINGS}
    for s in SETTINGS:
        provenance[s] = {'member': 'global/' + comparison_name(s),
                         'sha256': hashlib.sha256(raw[comparison_name(s)]).hexdigest(),
                         'canonical_json_sha256': c.digest(docs[s]),
                         'original_input_identity': authority['inputs'][s],
                         'runtime_input_identity': packets[s]['input_identity']}
        require(packets[s].get('cohort_authority', {}).get('AB') ==
                {'sha256': provenance[s]['sha256'], 'bytes': len(raw[comparison_name(s)])}, 'offline_AB_packet_join')
    candidate = c.compare(docs, provenance, matched_authority=authority)
    # Retain failed attempts for diagnosis, never render/replay them as science.
    if candidate['status'] != 'available':
        candidate['evidence_admissible'] = False
        with targets['candidate'].open('xb') as stream:
            stream.write(capture.encoded(candidate))
        raise ValueError('offline_C_first_failed_gate:' + str(candidate.get('first_failed_invariant')))
    result = replay.replay(packets, candidate, expected_head=expected_head, matched_authority=authority)
    result['input_bundle'] = provenance
    result['fresh_candidate_sha256'] = c.digest(candidate)
    render_science_pdf(result, targets['pdf'])
    result['pdf'] = capture.file_identity(targets['pdf'])
    for key, document in (('candidate', candidate), ('replay', result)):
        with targets[key].open('xb') as stream:
            stream.write(capture.encoded(document))
    return targets


def render_science_pdf(result, path, *, preview_dir=None):
    """Six pages; signed fit, residual, application, variance and domain views."""
    import math
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    require(len(result['parents']) == 6, 'science_six_parent_inventory')
    with Path(path).open('xb') as stream:
        with PdfPages(stream, metadata={'Title': 'Matched-cohort A/B/C signed-event comparison', 'CreationDate': None}) as pdf:
            for p in result['parents']:
                fig, axes = plt.subplots(4, 3, figsize=(16, 15), layout='constrained')
                metrics = '; '.join(f'{name}={value:.3g}' for name, value in p.get('descriptive_data_chi2', {}).items())
                phi = list(p['populations']['pion']['nommcuts']['children'])
                fig.suptitle(f"{p['setting']}/lowe t{p['t_bin']} | detached full-sample template A/B/C, seven unselected strengths; occupied phi bins {phi}\n"
                             'Fixed-weight errors only; held-out fitted estimator, bias truth and MSE selection unavailable; '
                             'setting-local phi support, no cross-setting normalization claim\n'
                             'Common-fit-bin descriptive chi2 (not calibrated goodness-of-fit): ' + metrics, fontsize=10)
                x = [(a + b) / 2 for a, b in zip(p['edges'], p['edges'][1:])]
                observed = p['populations']['pion']['nommcuts']['parent']
                for name, y in p['pion_fit_models'].items():
                    axes[0, 0].plot(x, y[1:-1], label=name, lw=1)
                    axes[0, 1].plot(x, [(a - b) / e if e > 0 else math.nan
                                        for a, b, e in zip(observed['signed'][1:-1], y[1:-1],
                                                           observed['fixed_weight_error'][1:-1])], lw=1)
                    axes[0, 2].plot(x, y[1:-1], lw=1)
                for ax in (axes[0, 0], axes[0, 2]):
                    ax.errorbar(x, observed['signed'][1:-1], observed['fixed_weight_error'][1:-1], fmt='k.', ms=2)
                axes[0, 0].set(title='Pion control: fitted MM', xlim=(p['edges'][0], p['edges'][-1]))
                axes[0, 1].set(title='(control - model) / fixed input error', xlim=(.7, 1.3))
                axes[0, 1].axhline(0, color='black', lw=.5)
                axes[0, 2].set(title='Pion control: Lambda detail', xlim=(1.10, 1.16))
                visible_control = [y for name, values in p['pion_fit_models'].items()
                                   for xx, y in zip(x, values[1:-1]) if 1.10 <= xx <= 1.16]
                visible_control += [y + sign*e for xx, y, e in zip(x, observed['signed'][1:-1],
                                    observed['fixed_weight_error'][1:-1]) if 1.10 <= xx <= 1.16 for sign in (-1, 1)]
                if visible_control:
                    low, high = min(0., min(visible_control)), max(0., max(visible_control))
                    pad = max((high-low)*.1, 1e-10)
                    axes[0, 2].set_ylim(low-pad, high+pad)
                views = ((axes[1, 0], 'nommcuts', (p['edges'][0], p['edges'][-1])),
                         (axes[1, 1], 'allcuts', (1.10, 1.16)), (axes[1, 2], 'allcuts', (1.105, 1.135)),
                         (axes[2, 0], 'nommcuts', (.94, 1.03)), (axes[2, 1], 'nommcuts', (.7, 1.3)))
                for ax, cut, window in views:
                    for name, model in p['models'].items():
                        applied = model['production_convention_application'][cut]['parent']
                        ax.errorbar(x, applied['signed'][1:-1], applied['fixed_weight_error'][1:-1],
                                    lw=.7, elinewidth=.4, label=name)
                        if model['convention_differs']:
                            ax.plot(x, model['application'][cut]['parent']['signed'][1:-1], ':', lw=.8,
                                    label=name + ' raw K/P')
                    for label, y in p['kaon_fit_models'].items():
                        ax.plot(x, y[1:-1], '--', lw=1.3, label='fitted K_' + label)
                    ax.set(title=cut + ': event replay vs fitted K; dotted=raw K/P', xlim=window)
                    visible = [float(y) for line in ax.lines for xx, y in zip(line.get_xdata(), line.get_ydata())
                               if window[0] <= xx <= window[1] and math.isfinite(y)]
                    visible += [y + sign*e for model in p['models'].values()
                                for branch in ('application', 'production_convention_application')
                                for xx, y, e in zip(x, model[branch][cut]['parent']['signed'][1:-1],
                                                   model[branch][cut]['parent']['fixed_weight_error'][1:-1])
                                if window[0] <= xx <= window[1] for sign in (-1, 1)]
                    if visible:
                        low, high = min(0, min(visible)), max(0, max(visible))
                        pad = max((high - low) * .1, 1e-10)
                        ax.set_ylim(low - pad, high + pad)
                names = list(p['models'])
                for cut, marker in (('nommcuts', 'o-'), ('allcuts', 'x--')):
                    axes[2, 2].plot(range(len(names)), [p['models'][n]['application'][cut]['parent']['domains']
                                                      ['signal_1.10_1.16']['sumw2'] for n in names], marker, label=cut)
                    if any(p['models'][n]['convention_differs'] for n in names):
                        axes[2, 2].plot(range(len(names)), [p['models'][n]['production_convention_application'][cut]['parent']
                            ['domains']['signal_1.10_1.16']['sumw2'] for n in names], ':', label=cut + ' floor/clip')
                axes[2, 2].set_xticks(range(len(names)), names, rotation=60, fontsize=6)
                axes[2, 2].set_title('Lambda fixed-weight Sumw2 (no optimal lambda)')
                axes[2, 2].set_yscale('symlog', linthresh=1e-12)
                axes[2, 2].legend(fontsize=7)
                for cut, marker in (('nommcuts', 'o-'), ('allcuts', 'x--')):
                    differences = [p['models'][n]['paired_difference_from_B']['production_convention'][cut]['parent']
                                   ['domains']['signal_1.10_1.16'] for n in names]
                    axes[3, 0].errorbar(range(len(names)), [v['signed'] for v in differences],
                                       [v['fixed_weight_error'] for v in differences], fmt=marker, label=cut)
                    axes[3, 1].plot(range(len(names)), [p['models'][n]['production_convention_application'][cut]['parent']
                        ['domains']['normal_full_MM']['largest_one_percent_Sumw2_fraction'] for n in names], marker, label=cut)
                for ax in (axes[3, 0], axes[3, 1]):
                    ax.set_xticks(range(len(names)), names, rotation=60, fontsize=6)
                    ax.legend(fontsize=7)
                axes[3, 0].axhline(0., color='black', lw=.5)
                axes[3, 0].set_title('Lambda-window difference vs B; paired fixed-weight error')
                axes[3, 1].set_title('Full normal MM: largest 1% contribution to Sumw2')
                axes[3, 2].axis('off')
                notes = ['Source and parent/phi-child closures: PASS', 'No positive unsupported denominator admitted',
                         'Allcuts and nommcuts remain separate', 'Raw K/P and production floor/clip both retained',
                         'Window totals use inclusive normal-bin centers', 'Paired difference is NOT an extracted yield',
                         'Fits, priors and K_B use the full cohort', 'No held-out estimator / bias truth / MSE optimum',
                         'No covariance from refitting / repeated training']
                axes[3, 2].text(0., 1., '\n'.join(notes), va='top', fontsize=8)
                axes[0, 0].legend(fontsize=6, ncol=3)
                axes[1, 0].legend(fontsize=6, ncol=3)
                for ax in axes.flat:
                    ax.tick_params(labelsize=7)
                    ax.set_xlabel('model' if ax in (axes[2, 2], axes[3, 0], axes[3, 1]) else 'shifted MM (GeV)')
                    ax.set_ylabel('signed normalized bin content')
                axes[0, 1].set_ylabel('standardized residual')
                axes[2, 2].set_ylabel('fixed-weight Sumw2')
                axes[3, 0].set_ylabel('signed normalized window total difference')
                axes[3, 1].set_ylabel('fraction of quadratic support')
                pdf.savefig(fig)
                if preview_dir is not None:
                    fig.savefig(Path(preview_dir) / f"{p['setting']}-t{p['t_bin']}.png", dpi=90)
                plt.close(fig)
            require(pdf.get_pagecount() == 6, 'science_PDF_page_count')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-commit', required=True)
    parser.add_argument('--offline-bundle', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--repo', type=Path)
    parser.add_argument('--outdir', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    require(re.fullmatch('[0-9a-f]{40}', args.source_commit) is not None, 'source_commit_requires_full_sha')
    try:
        if args.offline_bundle:
            require(args.output_dir and not any((args.repo, args.outdir, args.output)), 'offline_argument_scope')
            outputs = offline_science(args.offline_bundle, args.output_dir, args.source_commit)
            print(json.dumps({k: str(v) for k, v in outputs.items()}, sort_keys=True))
        else:
            require(not args.output_dir and all((args.repo, args.outdir, args.output)), 'capture_argument_scope')
            print(execute_gate(args.repo, args.outdir, args.output, args.source_commit).as_posix())
    except Exception as exc:
        print('Matched-cohort gate failed: ' + str(exc), file=sys.stderr)
        if args.outdir and args.output:
            print('Owner gate status: ' + str(isolation.accepted.gate_status_path(args.outdir, args.output)), file=sys.stderr)
        return 1
    return 0



def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def audit_default_output_writer(repo):
    """Locate the first non-Plot output writer on the unchanged -d path.

    This is source evidence, not an installed resolver attestation. In
    particular a Plot_LTSep probe cannot certify a no-mode Root invocation.
    """
    repo = Path(repo)
    launcher = (repo / "run_Prod_Analysis.sh").read_text(encoding="utf-8")
    caller = repo / "src/setup/findEffectiveCharge.py"
    text = caller.read_text(encoding="utf-8")
    tree = ast.parse(text)
    calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call)
             and isinstance(node.func, ast.Name) and node.func.id == "create_ltsep_root"]
    require(len(calls) == 1 and len(calls[0].args) == 1 and not calls[0].keywords,
            "effective_charge_default_root_call_changed")
    require(re.search(r"python3 findEffectiveCharge\.py", launcher) is not None,
            "effective_charge_launcher_call_missing")
    require('foutcsv = OUTPATH + "/" + OutFilename + ".csv"' in text,
            "effective_charge_output_expression_changed")
    helper = (repo / "src/setup/getDataTable.py").read_text(encoding="utf-8")
    require('data_to_csv(foutcsv, "Total Efficiency", tot_efficiency, runNum)' in helper,
            "effective_charge_csv_writer_changed")
    utility = (repo / "src/utility/utility.py").read_text(encoding="utf-8")
    require('def data_to_csv(' in utility, "effective_charge_csv_sink_missing")
    return {"caller": caller.relative_to(repo).as_posix(), "mode": [],
            "writer": "src/setup/getDataTable.py -> src/utility/utility.py::data_to_csv",
            "destination": "OUTPATH/table_kaon_<KIN>.csv",
            "required_invariant": "effective_charge_default_root_private_outpath_verified",
            "private_input_writer": "src/cuts/shift_prep.py::_apply_shift_branches -> src/setup/shift_MM.py::add_derived_branches_to_file",
            "private_input_scope": "MM_shift and t_shift in private prompt trees; original ROOT container identities remain distinct",
            "root_semantics": "source reviewed only; real ROOT clone/update behavior requires the user farm gate"}


if __name__ == '__main__':
    raise SystemExit(main())
