"""Real ledger writer and extracted main consumers; no full main/ROOT execution."""
import ast
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src/utility'))
import frozen_manifest as manifests
from epsilon_correction_compare import load_ledgers_and_compare


def main_tree():return ast.parse((ROOT/'src/main.py').read_text(encoding='utf-8'))


def high_block(tree):
    return next(n for n in tree.body if isinstance(n,ast.If)and any(isinstance(c,ast.Call)
        and isinstance(c.func,ast.Name)and c.func.id=='build_frozen_manifest_payload'for c in ast.walk(n)))


def input_nodes(tree):
    body=high_block(tree).body
    start=next(i for i,n in enumerate(body)if isinstance(n,ast.Assign)and any(isinstance(t,ast.Name)and t.id=='zeroth_bundle'for t in n.targets))
    end=next(i for i,n in enumerate(body)if isinstance(n,ast.With)and any(isinstance(a,ast.Assign)
        and any(isinstance(t,ast.Name)and t.id=='high_bg_summary'for t in a.targets)for a in n.body))
    return body[start:end+1]


def step6_nodes(tree):
    start=next(i for i,n in enumerate(tree.body)if isinstance(n,ast.Assign)and any(isinstance(t,ast.Name)and t.id=='correction_ledger_payload'for t in n.targets))
    end=next(i for i,n in enumerate(tree.body[start:],start)if isinstance(n,ast.Expr)and isinstance(n.value,ast.Call)
        and isinstance(n.value.func,ast.Name)and n.value.func.id=='debug_file_status'
        and n.value.args and isinstance(n.value.args[0],ast.Constant)and n.value.args[0].value=='Correction ledger JSON')
    return tree.body[start:end+1]


def execute(nodes,context):exec(compile(ast.Module(body=nodes,type_ignores=[]),'actual-main-ledger-statements','exec'),context)


def ledger_module():
    def unused_hist_integrator(*args,**kwargs):raise AssertionError('ROOT histogram integration is outside this fixture')
    with patch.dict(sys.modules,{'utility':types.SimpleNamespace(integrate_hist_range=unused_hist_integrator)}):
        spec=importlib.util.spec_from_file_location('_actual_correction_ledger',ROOT/'src/utility/correction_ledger.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


class MainCorrectionLedgerPaths(unittest.TestCase):
    @contextmanager
    def written_pair(self,profile):
        ledger=ledger_module()
        with tempfile.TemporaryDirectory()as directory:
            out=Path(directory);payloads={};written={};recorded={}
            for eps in ('low','high'):
                filename='FullAnalysis_Q4p4W2p74_'+eps+'e'
                inp=dict(ParticleType='kaon',Q2='4p4',W='2p74',EPSSET=eps,OutFilename=filename,
                         mm_min=1.09,mm_max=1.15,bg_active_profile=profile)
                manifests.write_zeroth_iteration_input_bundle(str(out),'kaon','4p4','2p74',eps,[],inp,active_profile=profile)
                context=dict(histlist=[],inpDict=inp,OUTPATH=str(out),ParticleType='kaon',OutFilename=filename,
                    output_file_lst=[],build_correction_ledger=ledger.build_correction_ledger,
                    write_correction_ledger=ledger.write_correction_ledger,get_correction_ledger_paths=manifests.get_correction_ledger_paths,
                    debug_file_status=lambda *args:None)
                execute(step6_nodes(main_tree()),context)
                payloads[eps]=json.loads(json.dumps(context['correction_ledger_payload']));written[eps]=context['correction_ledger_paths'];recorded[eps]=inp['correction_ledger_path']
                (out/('kaon_'+filename+'_bg_opt_'+eps+'.json')).write_text(json.dumps({'mode':eps,'active_profile':profile}))
            context=dict(OUTPATH=str(out),ParticleType='kaon',Q2='4p4',W='2p74',OutFilename='FullAnalysis_Q4p4W2p74_highe',
                inpDict={'bg_active_profile':profile},get_correction_ledger_paths=manifests.get_correction_ledger_paths,
                load_zeroth_iteration_input_bundle=manifests.load_zeroth_iteration_input_bundle,json=json)
            yield out,context,payloads,written,recorded

    def test_actual_writer_step6_and_both_step7_consumers_agree_for_profiles(self):
        for profile in ('no_empirical_residual','nominal_weighted',None,''):
            with self.subTest(profile=profile),self.written_pair(profile)as(out,context,payloads,written,recorded):
                for eps in ('low','high'):
                    self.assertEqual(recorded[eps],written[eps][0])
                    self.assertTrue(Path(recorded[eps]).is_file())
                    generic=manifests.get_correction_ledger_paths(str(out),'kaon','FullAnalysis_Q4p4W2p74_'+eps+'e',profile)['json']
                    if profile=='no_empirical_residual':self.assertFalse(Path(generic).exists())
                execute(input_nodes(main_tree()),context)
                for eps in ('low','high'):
                    self.assertEqual(context[eps+'_ledger_path'],written[eps][0])
                    self.assertEqual(context[eps+'_ledger_payload'],payloads[eps])
                    self.assertEqual(context[eps+'_bg_summary'],dict(mode=eps,active_profile=profile))
                compare=load_ledgers_and_compare(context['low_ledger_path'],context['high_ledger_path'])
                self.assertEqual(compare['active_profile'],profile)
                self.assertEqual(compare['failure_bins'],[])

    def test_stale_generic_ledger_is_ignored(self):
        with self.written_pair('no_empirical_residual')as(out,context,payloads,written,recorded):
            for eps in ('low','high'):
                generic=manifests.get_correction_ledger_paths(str(out),'kaon','FullAnalysis_Q4p4W2p74_'+eps+'e','no_empirical_residual')['json']
                Path(generic).write_text(json.dumps({'stale_wrong_profile':True}))
            execute(input_nodes(main_tree()),context)
            for eps in ('low','high'):self.assertEqual(context[eps+'_ledger_payload'],payloads[eps])

    def test_missing_or_malformed_profile_ledgers_fail_without_generic_fallback(self):
        for eps in ('low','high'):
            for failure in ('missing','malformed'):
                with self.subTest(eps=eps,failure=failure),self.written_pair('no_empirical_residual')as(out,context,payloads,written,recorded):
                    primary=Path(written[eps][0]);generic=manifests.get_correction_ledger_paths(str(out),'kaon','FullAnalysis_Q4p4W2p74_'+eps+'e','no_empirical_residual')['json']
                    Path(generic).write_text(json.dumps(payloads[eps]))
                    if failure=='missing':primary.unlink()
                    else:primary.write_text('{')
                    with self.assertRaises(FileNotFoundError if failure=='missing'else json.JSONDecodeError):
                        execute(input_nodes(main_tree()),context)

    def test_required_finalization_inputs_precede_expensive_separation(self):
        tree=main_tree();inputs=input_nodes(tree)
        separation=next(n for n in ast.walk(high_block(tree))if isinstance(n,ast.Call)and isinstance(n.func,ast.Name)
            and n.func.id=='run_bash_script'and 'run_xsect.sh'in ast.unparse(n))
        self.assertLess(inputs[-1].end_lineno,separation.lineno)
        finalizer_start=next(n for n in high_block(tree).body if isinstance(n,ast.Assign)
            and any(isinstance(t,ast.Name)and t.id=='manifest_payload'for t in n.targets))
        body=high_block(tree).body;timer=body[body.index(finalizer_start)-1]
        self.assertEqual(ast.unparse(timer),'stage_start = perf_counter()')


if __name__=='__main__':unittest.main()
