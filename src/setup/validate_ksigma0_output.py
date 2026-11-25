#!/usr/bin/env python3
"""Lossless generation adapter only; no production weights, bindings or cuts.

ROOT is imported only by runtime probe or validation. Pure helpers support local
deterministic tests; they do not establish ROOT, detector or physical acceptance.
Coordinates follow the unchanged archived macro: Q2/t GeV^2 (t signed), W/MM
GeV, phipq radians, target/focal positions cm, slopes radians, delta percent.
"""
from array import array
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import struct


def require(ok, message):
    if not ok:
        raise ValueError(message)


def number(text, pattern):
    matches = re.findall(pattern + r'\s*=\s*([-+0-9.eEdD]+)', text, re.I)
    require(len(matches) == 1, 'missing/ambiguous sidecar field: ' + pattern)
    value = float(matches[0].replace('D', 'e').replace('d', 'e'))
    require(math.isfinite(value), 'nonfinite sidecar field')
    return value


def parse_hist(text, requested=200000):
    require('producing a SIGMA0' in text and "H(e,e'K)" in text, 'wrong reaction sidecar')
    for flag, expected in {'doing_kaon': 'T', 'doing_pion': 'F', 'doing_delta': 'F',
                           'doing_semi': 'F', 'doing_phsp': 'F', 'doing_hydkaon': 'T',
                           'doing_decay': 'T', 'using_rad': 'T', 'mc_smear': 'T',
                           'using_Eloss': 'T', 'correct_Eloss': 'T'}.items():
        values = re.findall(r'\b' + flag + r'\s*=\s*([TF])\b', text)
        require(values == [expected], 'wrong/missing reaction flag: ' + flag)
    require(number(text, r'\bwhich_kaon') == 1, 'wrong recoil selector')
    ngen = number(text, r'\bNgen\s*\(request\)')
    ncontribute = number(text, r'\bNcontribute')
    normfac = number(text, r'\bnormfac')
    require(ngen == requested and ncontribute == requested and normfac > 0,
            'incomplete generation/normalization fields')
    mass = re.findall(r'Mrec_struc\s*=\s*([-+0-9.eEdD]+)\s*MeV', text)
    require(len(mass) == 1 and abs(float(mass[0].replace('D', 'e')) - 1192.64) < 0.02,
            'Sigma0 recoil mass not verified')
    return {'Ngen_requested': int(ngen), 'Ncontribute': int(ncontribute), 'normfac': normfac,
            'recoil_MeV': float(mass[0]), 'Ebeam_MeV': number(text, r'\bEbeam')}


FP = ('hsdelta', 'hsxptar', 'hsyptar', 'hsytar', 'hsxfp', 'hsxpfp', 'hsyfp', 'hsypfp',
      'ssdelta', 'ssxptar', 'ssyptar', 'ssytar', 'ssxfp', 'ssxpfp', 'ssyfp', 'ssypfp')
IMMUTABLE = FP + ('Weight', 'fry', 'radphot', 'Q2i', 'Wi', 'ti', 'Mhadron')
REQUIRED = FP + ('Weight', 'Q2', 'W', 't', 'phipq', 'missmass',
                 'phgcer_x_det', 'phgcer_y_det', 'paero_x_det', 'paero_y_det')


def f32(x):
    return struct.unpack('f', struct.pack('f', x))[0]


def detector_axis(position, slope, z):
    # All three operands in the archived common buffer are Float_t. Match
    # multiplication and addition rounding, not a double-precision surrogate.
    return f32(position + f32(f32(z) * slope))


def tray(row, stem):
    # Archived buffer stores these derived axes as Float_t before comparison.
    x = detector_axis(row['ssxfp'], row['ssxpfp'], 231.0)
    y = detector_axis(row['ssyfp'], row['ssypfp'], 231.0)
    tight = any(k in stem for k in ('Q4p4W2p74', 'Q3p0W3p14', 'Q5p5W3p02'))
    xlim, ylim = (45., 30.) if tight else (55., 50.)
    return -xlim < x < xlim and -ylim < y < ylim


def check_rows(raw, reconstructed, stem):
    """Independent surviving membership/order and unchanged detector/Weight check."""
    reconstructed = iter(reconstructed)
    weights = []
    membership = hashlib.sha256()
    raw_count = 0
    for i, before in enumerate(raw):
        raw_count += 1
        if not tray(before, stem):
            continue
        after = next(reconstructed, None)
        require(after is not None, 'surviving count mismatch/empty')
        for name in IMMUTABLE:
            require(math.isfinite(before[name]) and before[name] == after[name],
                    'surviving event identity/Weight mismatch: ' + name)
        membership.update(struct.pack('<Q', i))
        weights.append(after['Weight'])
    require(weights and next(reconstructed, None) is None, 'surviving count mismatch/empty')
    sumw2 = math.fsum(w * w for w in weights)
    require(math.isfinite(sumw2) and sumw2 > 0, 'invalid/zero Weight Sumw2')
    return {'raw_entries': raw_count, 'surviving_entries': len(weights),
            'membership_sha256': membership.hexdigest(), 'sum_Weight': math.fsum(weights),
            'sumw2_Weight': sumw2, 'native_Weight_preserved': True}


def branch_names(tree):
    names = [b.GetName() for b in tree.GetListOfBranches()]
    require(len(set(names)) == len(names), 'duplicate ROOT branches')
    return names


def scalar_types(tree, names):
    for name in names:
        leaf = tree.GetLeaf(name)
        require(leaf and leaf.GetLenStatic() == 1 and not leaf.GetLeafCount() and
                leaf.GetTypeName() in ('Float_t', 'Double_t'), 'missing/incompatible scalar branch: ' + name)


def rows(tree, names):
    """Stream entries; do not materialize hundreds of MB of per-row dictionaries."""
    scalar_types(tree, names)
    for i in range(tree.GetEntries()):
        require(tree.GetEntry(i) > 0, 'unreadable ROOT entry')
        row = {n: float(getattr(tree, n)) for n in names}
        require(all(math.isfinite(v) for v in row.values()), 'nonfinite ROOT branch')
        yield row


def open_tree(ROOT, path):
    handle = ROOT.TFile.Open(str(path), 'READ')
    require(handle and not handle.IsZombie(), 'unreadable ROOT artifact')
    tree = handle.Get('h10')
    require(tree and tree.InheritsFrom('TTree') and tree.GetEntries() > 0, 'missing/empty ROOT h10')
    return handle, tree


def clone_identity(ROOT, source, target):
    """Separate complete clone, one exact Weight identity branch; never reweight."""
    require(not target.exists(), 'consumer candidate already exists')
    require('iter_weight' not in branch_names(source), 'unexpected existing iter_weight')
    output = ROOT.TFile.Open(str(target), 'CREATE')
    require(output and not output.IsZombie(), 'cannot create separate consumer candidate')
    try:
        output.cd()
        cloned = source.CloneTree(0)
        value = array('d', [0.])
        branch = cloned.Branch('iter_weight', value, 'iter_weight/D')
        require(branch, 'iter_weight creation failed')
        for i in range(source.GetEntries()):
            require(source.GetEntry(i) > 0, 'source entry unreadable')
            value[0] = float(source.Weight)
            require(cloned.Fill() > 0, 'clone fill failed')
        require(cloned.Write() > 0, 'clone write failed')
    finally:
        output.Close()


def root_identity(ROOT):
    """Identity of the imported validator runtime, without reading event data."""
    paths = set()
    with open('/proc/self/maps', encoding='utf-8') as maps:
        for line in maps:
            if '/libCore.so' in line:
                paths.add(str(Path(line.split(None, 5)[5].strip()).resolve(strict=True)))
    require(len(paths) == 1, 'PyROOT loaded core library unavailable/ambiguous')
    core = paths.pop()
    module = str(Path(ROOT.__file__).resolve(strict=True))
    def sha(name):
        h = hashlib.sha256()
        with open(name, 'rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                h.update(block)
        return h.hexdigest()
    return {'version_int': int(ROOT.gROOT.GetVersionInt()), 'module': module,
            'module_sha256': sha(module), 'core': core, 'core_sha256': sha(core)}


def validate_attempt(root, receipt):
    import ROOT  # Actual facility runtime only. Local tests use deterministic stubs.
    require(root_identity(ROOT) == receipt['root_runtime']['pyroot'],
            'validator ROOT identity changed since pre-generation admission')
    stem = receipt['stem']
    recon_dir = root / 'recon_hcana'
    hist = parse_hist((root / 'OUTPUTS' / (stem + '.hist')).read_text())
    raw_file, raw = open_tree(ROOT, recon_dir / (stem + '.root'))
    recon_file, reconstructed = open_tree(ROOT, recon_dir / (stem + '_recon.root'))
    candidate = root / 'OUTPUTS' / (stem + '_consumer.root')
    candidate_file = None
    try:
        result = check_rows(rows(raw, IMMUTABLE), rows(reconstructed, IMMUTABLE), stem)
        require(result['raw_entries'] == hist['Ncontribute'], 'converter/generator count mismatch')
        support = {n: [math.inf, -math.inf] for n in ('Q2', 'W', 't', 'phipq', 'missmass')}
        for row in rows(reconstructed, REQUIRED):
            require(0 < row['Q2'] < 4 * (hist['Ebeam_MeV'] / 1000.) ** 2 and
                    0 < row['W'] < hist['Ebeam_MeV'] / 1000. + 0.938272 and
                    row['t'] < 0 and -math.pi <= row['phipq'] <= math.pi and
                    row['missmass'] > 0, 'incompatible kinematic coordinates/support')
            for axis in ('x', 'y'):
                pos, slope = ('ssxfp', 'ssxpfp') if axis == 'x' else ('ssyfp', 'ssypfp')
                require(row['phgcer_' + axis + '_det'] == detector_axis(row[pos], row[slope], 156.27) and
                        row['paero_' + axis + '_det'] == detector_axis(row[pos], row[slope], 231.0),
                        'detector coordinate definition mismatch')
            for name, bounds in support.items():
                bounds[0], bounds[1] = min(bounds[0], row[name]), max(bounds[1], row[name])
        clone_identity(ROOT, reconstructed, candidate)
        candidate_file, copy = open_tree(ROOT, candidate)
        names = branch_names(reconstructed)
        require(set(branch_names(copy)) == set(names) | {'iter_weight'}, 'clone branch inventory mismatch')
        original_rows = rows(reconstructed, names)
        copied_rows = rows(copy, names + ['iter_weight'])
        require(reconstructed.GetEntries() == copy.GetEntries(), 'identity clone count mismatch')
        scalar_types(copy, ['iter_weight'])
        for name in names:
            require(reconstructed.GetLeaf(name).GetTypeName() == copy.GetLeaf(name).GetTypeName(),
                    'identity clone branch type mismatch')
        copied_weights = []
        for orig, copied in zip(original_rows, copied_rows):
            require(all(copied[n] == v for n, v in orig.items()) and
                    copied['iter_weight'] == orig['Weight'], 'lossy clone/weight mapping')
            copied_weights.append(copied['iter_weight'])
        require(math.fsum(w ** 2 for w in copied_weights) == result['sumw2_Weight'],
                'identity Sumw2 mismatch')
        result.update(histogram=hist, iter_weight_identity=True,
                      units={'Q2': 'GeV^2', 't': 'signed GeV^2', 'W': 'GeV',
                             'missmass': 'GeV; kaon hypothesis 0.493677', 'phipq': 'radians'},
                      support=support,
                      scientific_acceptance=False, consumer_binding_exported=False)
        return result
    finally:
        if candidate_file:
            candidate_file.Close()
        recon_file.Close()
        raw_file.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest='mode', required=True)
    modes.add_parser('probe')
    validate = modes.add_parser('validate')
    validate.add_argument('root', type=Path)
    args = parser.parse_args()
    if args.mode == 'probe':
        import ROOT
        result = root_identity(ROOT)
    else:
        receipt = json.loads((args.root / 'receipt.json').read_text())
        result = validate_attempt(args.root, receipt)
    print('KSIGMA0_VALIDATOR=' + json.dumps(result, sort_keys=True, allow_nan=False))


if __name__ == '__main__':
    main()
