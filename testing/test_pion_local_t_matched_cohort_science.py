"""Synthetic signed cohort -> seven fresh fits -> six-parent replay.

No physical truth, ROOT, farm launch or accepted historical artifact is used.
"""
from copy import deepcopy
import json
import hashlib
from pathlib import Path
import tempfile
import subprocess
import unittest
from unittest.mock import patch
import zipfile

import numpy as np
from scipy.optimize import least_squares

from testing import pion_local_t_event_validation_capture as cap
from testing import pion_local_t_event_validation_replay as replay
from testing import pion_local_t_transfer_variance_comparison as c
from testing import run_pion_local_t_matched_cohort_gate as owner
from testing.test_pion_local_t_event_validation import histogram, fixture, comparison_fixture
from testing.test_pion_local_t_transfer_variance_comparison import synthetic

HEAD = 'b' * 40


def scientific_fixture():
    docs, packets = {}, {}
    q = synthetic()
    # Exact signed source counts, including cancellation and quadratic Sumw2.
    coeff = cap.coefficients(2., .5, 2.)
    counts = dict(zip(cap.SOURCES, (3, 2, 1, 1)))
    signed = sum(coeff[s] * counts[s] for s in cap.SOURCES)
    variance = sum(coeff[s] ** 2 * counts[s] for s in cap.SOURCES)
    q['y'] = np.r_[0., np.full(12, signed), 0.]
    q['s2'] = np.r_[0., np.full(12, variance), 0.]
    r = least_squares(lambda b: c.base_residual(q, b), q['initial'], bounds=(np.zeros(3), np.full(3, np.inf)), method='trf')
    base = c.base_residual(q, r.x)
    def h(y, v=None):
        # Existing fixture fingerprint builder uses its declared EDGES; patch
        # only that helper's input so the production validator remains real.
        with patch('testing.test_pion_local_t_event_validation.EDGES', q['edges'].tolist()):
            return histogram(list(map(float, y)), None if v is None else list(map(float, v)))
    for s in c.SETTINGS:
        packet, _, _ = fixture(s)
        packet['source_HEAD'] = HEAD
        packet['events'] = {'pion': [], 'kaon': []}
        parents = []
        for ti in range(3):
            parent = deepcopy(packet['parents'][ti])
            p = {k: deepcopy(parent[k]) for k in ('t_bin_index', 't_edges', 'pion_parent_id',
                 'canonical_interval_pair_id', 'canonical_interval_pair_hash', 'coordinate_fingerprint', 'baseline_parent_fingerprint')}
            data = h(q['y'], q['s2'])
            fit = {'pion_component_alignment': {'fixed': True}, 'resolved_subtraction_config': {},
                   'H_pion_control_input': data, 'H_kaon_nosub_input': deepcopy(data),
                   'H_pion_fit_total_refined': h(q['t'] @ r.x),
                   'H_kaon_pion_bg_fit_total_refined': h(q['k']),
                   'diagnostics': {'pion': {'fit_min': .7, 'fit_max': 1.3, 'exclude_windows': [],
                       'staged_component_amplitudes': dict(zip(q['names'], q['initial'].tolist())),
                       'joint_refinement': {'success': True, 'status': 'success', 'regularization_enabled': True,
                           'coefficients': dict(zip(q['names'], r.x.tolist())), 'fit_bin_indices': q['ix'].tolist(),
                           'regularization_widths': dict(zip(q['names'], q['tau'].tolist())),
                           'data_chi2_contribution': float(base[:-3] @ base[:-3]),
                           'regularization_contribution': float(base[-3:] @ base[-3:])}}}}
            templates = {}
            for j, n in enumerate(q['names']):
                column = q['t'][:, j]
                normalized = h(column / sum(column))
                templates[n] = {'status': 'genuinely_local', 'normalization_count': 1,
                    'canonical_selection_closure': True, 'weighted_child_to_parent_closure': True,
                    'production_global_selection_closure': True, 't_edges': p['t_edges'],
                    'provenance': {'status': 'available', 'sha256': 'f'*64, 'tree_name': 'synthetic', 'missing_branches': []},
                    'histogram': normalized, 'normalized_values': normalized['values'],
                    'normalized_variances': normalized['stored_sumw2'], 'normalized_integral': 1.}
                fit['H_simc_shape_' + n] = h(column)
                fit['H_kaon_fit_' + n + '_scaled_refined'] = h(column * .1 if j == 0 else column * (.2 if j == 1 else .1))
            p.update(status='available', authoritative_parent_unchanged=True, production_applied=False,
                     rejection_reasons=[], alignment_mode='fixed_parent_alignment_template_shape_test',
                     baseline_A=deepcopy(fit), candidate_B=fit, templates=templates,
                     frozen_inputs={'pion': data, 'kaon': deepcopy(data)})
            parents.append(p)
            parent.update(pion_input=deepcopy(data), kaon_input=deepcopy(data))
            packet['parents'][ti] = parent
            for population in ('pion', 'kaon'):
                for source in cap.SOURCES:
                    for bi, mm in enumerate((q['edges'][:-1] + q['edges'][1:]) / 2):
                        for j in range(counts[source]):
                            role = 'dummy' if source.startswith('dummy') else 'data'
                            row = {'population': population, 'source_label': source,
                                   'entry_index': ti*1000 + bi*10 + j, 'source_sha256': packet['input_identity'][role]['sha256'],
                                   'coefficient': coeff[source], 'post_proton_factor': 1., 'proton_rf_accept': True,
                                   'MM': float(mm), 't': ti+.5, 'phi': 30., 't_index': ti, 'phi_index': 0,
                                   'allcuts': bool(1.10 <= mm <= 1.16), 'nommcuts': True,
                                   'coordinate_fingerprint': packet['coordinate_fingerprint']}
                            row['fold'] = cap.fold(row)
                            packet['events'][population].append(row)
        for population, rows in packet['events'].items():
            packet['source_counts'][population] = {s: sum(r['source_label'] == s for r in rows) for s in cap.SOURCES}
            packet['event_hashes'][population] = cap.digest(rows)
        cap.validate_packet(packet)
        packets[s] = packet
        docs[s] = {'schema_version': 'pion_t_template_comparison/v1', 'source_HEAD': HEAD,
                   'status': 'available', 'diagnostic_only': True, 'production_applied': False,
                   'source_identity_rejection_reasons': [], 'phi_setting': s,
                   'setting': {'EPSSET': 'low', 'ParticleType': 'kaon', 'Q2': '4p4', 'W': '2p74'},
                   'canonical_t_edges': packet['t_edges'], 'canonical_phi_edges': packet['phi_edges'],
                   'source_files': {}, 'original_input_population_accounting': {},
                   'original_analysis_input_identity': packet['input_identity'], 'parents': parents}
    authority = {'schema_version': 'pion_matched_cohort_authority/v1', 'mode': 'matched_cohort',
                 'source_HEAD': HEAD, 'attempt_id': 'synthetic-only', 'inputs': {s: packets[s]['input_identity'] for s in c.SETTINGS}}
    provenance = {s: {'sha256': c.digest(docs[s]), 'canonical_json_sha256': c.digest(docs[s])} for s in c.SETTINGS}
    for s in c.SETTINGS:
        packets[s]['cohort_authority'] = {'mode': 'matched_cohort', 'attempt_id': authority['attempt_id'],
            'manifest_sha256': cap.digest(authority), 'AB': {'sha256': provenance[s]['sha256'], 'bytes': len(c.encoded(docs[s]))}}
    return docs, packets, authority, provenance


class MatchedScience(unittest.TestCase):
    def test_unchanged_collector_builds_the_actual_matched_profile(self):
        docs, packets, authority, _ = scientific_fixture()
        raw = {owner.comparison_name(s): c.encoded(docs[s]) for s in c.SETTINGS}
        raw.update({owner.packet_name(s): cap.encoded(packets[s]) for s in c.SETTINGS})
        raw.update({owner.comparison_name(s).replace('.json', '.pdf'): b'%PDF-synthetic fixture' for s in c.SETTINGS})
        raw['matched-cohort-authority.json'] = cap.encoded(authority)
        collector = owner.ab.collector
        with tempfile.TemporaryDirectory() as d:
            directory = Path(d)
            evidence = directory / 'evidence'
            evidence.mkdir()
            for name, content in raw.items():
                (evidence / name).write_bytes(content)
            profile = directory / 'profile.json'
            profile.write_bytes(cap.encoded(owner.profile(HEAD)))
            self.assertEqual(collector.load_validation_profile(profile), owner.profile(HEAD))
            archive = directory / 'collected.zip'
            with patch.object(collector, 'collect_source_state', return_value=('synthetic source only', HEAD)), patch.object(
                collector, 'collect_source_checks', return_value=('synthetic checks only', [{'name': 'synthetic_only', 'returncode': 0}])), patch.object(
                collector, '_committed_identity', return_value=(True, [], [])):
                result = collector.collect_validation_bundle(outdir=evidence, kinematic='Q4p4W2p74', output=archive,
                    profile_path=profile, repo_root=directory)
            self.assertEqual(result['returncode'], 0)
            checked, saved, _ = owner.verify_archive(archive, HEAD)
            self.assertEqual(checked, raw)
            self.assertEqual(saved, authority)

    def test_historical_replay_and_C_math_match_unchanged_HEAD_source(self):
        root = Path(__file__).resolve().parents[1]
        def baseline(name, active):
            raw = subprocess.check_output(['git', 'show', 'HEAD:testing/' + name], cwd=root)
            namespace = {'__name__': 'testing._historical_control', '__package__': 'testing', '__file__': active.__file__}
            exec(compile(raw, 'HEAD:testing/' + name, 'exec'), namespace)
            return namespace
        packets, candidate = comparison_fixture()
        old = baseline('pion_local_t_event_validation_replay.py', replay)
        self.assertEqual(cap.encoded(old['replay'](packets, candidate, expected_head='a'*40)),
                         cap.encoded(replay.replay(packets, candidate, expected_head='a'*40)))
        docs, _, _, provenance = scientific_fixture()
        for setting, document in docs.items():
            document['source_HEAD'] = c.SOURCE
            provenance[setting] = {'sha256': c.digest(document), 'canonical_json_sha256': c.digest(document)}
        old = baseline('pion_local_t_transfer_variance_comparison.py', c)
        # Both versions report the same active implementation-file hash here;
        # every numerical/schema field is otherwise compared without projection.
        self.assertEqual(c.encoded(old['compare'](docs, provenance)), c.encoded(c.compare(docs, provenance)))

    def test_joint_capture_closure_and_negative_identity_gates(self):
        docs, packets, authority, _ = scientific_fixture()
        authority.update(started_ns=0, source_files={})
        with tempfile.TemporaryDirectory() as d:
            output = Path(d)
            for s in c.SETTINGS:
                abpath = output / owner.comparison_name(s)
                abpath.write_bytes(c.encoded(docs[s]))
                (output / abpath.with_suffix('.pdf').name).write_bytes(b'%PDF-synthetic')
                packets[s]['cohort_authority'].update(manifest_sha256=cap.digest(authority), AB=cap.file_identity(abpath))
                (output / owner.packet_name(s)).write_bytes(cap.encoded(packets[s]))
            (output / 'observer-open-audit-fixture.json').write_bytes(cap.encoded({
                'manifest_sha256': cap.digest(authority),
                'counts': {name: 2 if name.endswith('comparison.json') else 1 for name in owner.artifacts()}}))
            with patch.object(owner.ab, 'pdf_pages', return_value=6), patch.object(owner.ab, 'verify_payload'):
                records, _, _ = owner.verify_joint(output, authority, {})
                self.assertEqual(len(records), 6)
                bad = deepcopy(packets['Left'])
                bad['phi_edges'] = [0., 50., 120.]
                (output / owner.packet_name('Left')).write_bytes(cap.encoded(bad))
                with self.assertRaisesRegex(ValueError, 'joint_coordinate_edges'):
                    owner.verify_joint(output, authority, {})
                (output / owner.packet_name('Left')).write_bytes(cap.encoded(packets['Left']))
                (output / owner.packet_name('Center')).unlink()
                with self.assertRaisesRegex(ValueError, 'joint_artifact_missing'):
                    owner.verify_joint(output, authority, {})

    def test_six_parent_fresh_sweep_and_signed_replay(self):
        docs, packets, authority, provenance = scientific_fixture()
        before = cap.digest([docs, packets])
        candidate = c.compare(docs, provenance, matched_authority=authority)
        self.assertEqual(candidate['status'], 'available')
        self.assertEqual(len(candidate['parents']), 6)
        for parent in candidate['parents']:
            self.assertEqual([v['lambda'] for v in parent['C']], list(c.LAMBDAS))
            self.assertTrue(all(v['selected'] is False for v in parent['C']))
            self.assertEqual(len({v['fit_identity'] for v in parent['C']}), 7)
            self.assertEqual(len({v['model_identity'] for v in parent['C']}), 7)
        result = replay.replay(packets, candidate, expected_head=HEAD, matched_authority=authority)
        self.assertEqual(len(result['parents']), 6)
        self.assertIsNone(result['selected_lambda'])
        self.assertEqual(cap.digest([docs, packets]), before)
        for parent in result['parents']:
            for cut in ('allcuts', 'nommcuts'):
                spectra = parent['models']['C_1000']['application'][cut]
                self.assertEqual(spectra['closure'], 'PASS')
                self.assertGreater(spectra['parent']['domains']['signal_1.10_1.16']['sumw2'], 0)
            paired = parent['models']['C_0']['paired_difference_from_B']['ratio']['allcuts']['parent']
            self.assertAlmostEqual(paired['domains']['signal_1.10_1.16']['signed'], 0., places=7)
            self.assertLess(paired['domains']['signal_1.10_1.16']['sumw2'], 1e-12)
        with self.assertRaisesRegex(ValueError, 'accepted_C_policy'):
            replay.validate_candidate(candidate)
        with self.assertRaisesRegex(ValueError, 'source/schema'):
            c.validate_document(docs['Center'], 'Center')

    def test_mixed_cohort_and_histogram_rejection(self):
        docs, packets, authority, provenance = scientific_fixture()
        candidate = c.compare(docs, provenance, matched_authority=authority)
        for key, value in (('attempt_id', 'old'), ('manifest_sha256', '0'*64)):
            bad = deepcopy(packets)
            bad['Left']['cohort_authority'][key] = value
            with self.assertRaisesRegex(ValueError, 'matched_replay_authority'):
                replay.replay(bad, candidate, expected_head=HEAD, matched_authority=authority)
        bad = deepcopy(packets)
        bad['Left']['parents'][0]['pion_input']['values'][1] += 1
        with self.assertRaisesRegex(ValueError, 'fingerprint'):
            replay.replay(bad, candidate, expected_head=HEAD, matched_authority=authority)

    def test_private_runtime_source_is_not_relabelled_as_original_container(self):
        docs, packets, authority, provenance = scientific_fixture()
        runtime = deepcopy(authority['inputs'])
        authority['private_input_policy'] = 'frozen_shiftprep_private_copy'
        authority['inputs'] = {s: {role: {'sha256': ('e' if role == 'data' else 'f')*64, 'bytes': 90}
                                  for role in ('data', 'dummy')} for s in c.SETTINGS}
        for s in c.SETTINGS:
            provenance[s].update(original_input_identity=authority['inputs'][s], runtime_input_identity=runtime[s])
            packets[s]['cohort_authority'].update(manifest_sha256=cap.digest(authority),
                original_input_identity=authority['inputs'][s], runtime_input_identity=runtime[s])
        candidate = c.compare(docs, provenance, matched_authority=authority)
        result = replay.replay(packets, candidate, expected_head=HEAD, matched_authority=authority)
        self.assertEqual(len(result['parents']), 6)
        self.assertNotEqual(candidate['provenance']['Center']['original_input_identity'],
                            candidate['input_metadata']['Center']['original_analysis_input_identity'])
        bad = deepcopy(packets)
        bad['Center']['cohort_authority']['original_input_identity'] = runtime['Center']
        with self.assertRaisesRegex(ValueError, 'matched_replay_original_inputs'):
            replay.replay(bad, candidate, expected_head=HEAD, matched_authority=authority)

    def test_bundle_validation_offline_fresh_fits_report_and_collision(self):
        docs, packets, authority, provenance = scientific_fixture()
        raw = {owner.comparison_name(s): c.encoded(docs[s]) for s in c.SETTINGS}
        raw.update({owner.packet_name(s): cap.encoded(packets[s]) for s in c.SETTINGS})
        raw.update({owner.comparison_name(s).replace('.json', '.pdf'): b'%PDF-synthetic fixture only' for s in c.SETTINGS})
        raw['matched-cohort-authority.json'] = cap.encoded(authority)
        manifest = {'complete': True, 'errors': [], 'git_head': HEAD, 'required_analysis_commit': HEAD,
                    'validation_profile': owner.PROFILE_ID, 'requested_kinematic': 'Q4p4W2p74',
                    'requested_settings': owner.ab.SETTINGS, 'global_artifacts': {},
                    'required_analysis_commit_is_ancestor': True, 'allowed_committed_files': [],
                    'committed_files_after_required_analysis_commit': [],
                    'unexpected_committed_files_after_required_analysis_commit': [],
                    'source_checks': [{'name': 'synthetic_only', 'returncode': 0}],
                    'settings': [dict(s, kinematic='Q4p4W2p74', artifacts={}) for s in owner.ab.SETTINGS]}
        for row in owner.profile(HEAD)['artifacts']['global']:
            name = row['basename_template']
            manifest['global_artifacts'][row['key']] = {'archive_path': 'global/' + name,
                'sha256': hashlib.sha256(raw[name]).hexdigest(), 'byte_size': len(raw[name]),
                'status': 'exists', 'json_status': 'valid' if name.endswith('.json') else None}
        with tempfile.TemporaryDirectory() as d:
            directory = Path(d)
            archive = directory / 'synthetic.zip'
            with zipfile.ZipFile(archive, 'w') as z:
                z.writestr('manifest.json', cap.encoded(manifest))
                z.writestr('source_state.txt', 'synthetic source fixture only')
                z.writestr('source_checks.txt', 'synthetic source checks only')
                z.writestr('global/', b'')
                for s in owner.ab.SETTINGS:
                    z.writestr(s['phi'] + '_' + s['epsilon'] + '/', b'')
                for name, data in raw.items():
                    z.writestr('global/' + name, data)
            checked, saved, _ = owner.verify_archive(archive, HEAD)
            self.assertEqual(saved, authority)
            self.assertEqual(checked, raw)
            outputs = owner.offline_science(archive, directory, HEAD)
            result = json.loads(outputs['replay'].read_bytes())
            self.assertEqual(len(result['parents']), 6)
            self.assertTrue(outputs['pdf'].read_bytes().startswith(b'%PDF-'))
            self.assertIsNone(result['selected_lambda'])
            with self.assertRaisesRegex(ValueError, 'offline_output_collision'):
                owner.offline_science(archive, directory, HEAD)
            with self.assertRaisesRegex(ValueError, 'matched_zip_manifest'):
                owner.verify_archive(archive, 'c'*40)
            with zipfile.ZipFile(archive, 'a') as z:
                # A duplicated member cannot hide altered archived inputs.
                import warnings
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore', UserWarning)
                    z.writestr('global/' + owner.packet_name('Center'), b'{}')
            with self.assertRaisesRegex(ValueError, 'duplicates'):
                owner.verify_archive(archive, HEAD)


if __name__ == '__main__':
    unittest.main()
