"""Narrow AST audit of the Step-6 presentation consumer ordering."""
import ast
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


def canonical_main(source):
    """Normalize only owned path keys/input relocation and the prior SIMC move."""
    block = 'stage_start = perf_counter()\nyieldDict.update(find_yield_simc(histlist, inpDict))\nrecord_stage_time("Step 6 simc yields", stage_start)\n'
    assert source.count(block) == 1
    tree = ast.parse(source.replace(block, ""))
    ledgers = [n for n in ast.walk(tree) if isinstance(n, ast.Subscript)
        and isinstance(n.value, ast.Call) and isinstance(n.value.func, ast.Name)
        and n.value.func.id == "get_correction_ledger_paths"]
    assert len(ledgers) == 3
    for node in ledgers:
        if isinstance(node.slice, ast.Constant) and node.slice.value in ("json", "json_profile"):
            node.slice = ast.Constant(value="json")
    high = next(n for n in tree.body if isinstance(n, ast.If) and any(isinstance(c, ast.Call)
        and isinstance(c.func, ast.Name) and c.func.id == "build_frozen_manifest_payload" for c in ast.walk(n)))
    start = next(i for i,n in enumerate(high.body) if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "zeroth_bundle" for t in n.targets))
    end = next(i for i,n in enumerate(high.body) if isinstance(n, ast.With) and any(isinstance(a, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "high_bg_summary" for t in a.targets) for a in n.body))
    inputs = high.body[start:end+1]
    del high.body[start:end+1]
    position = next(i for i,n in enumerate(high.body) if isinstance(n, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "manifest_payload" for t in n.targets))
    high.body[position:position] = inputs
    return ast.dump(tree)


class MainOrderTests(unittest.TestCase):
    def test_single_existing_producers_precede_single_per_kaon_finalizer(self):
        tree = ast.parse((ROOT / "src/main.py").read_text(encoding="utf-8"))
        calls = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                calls.setdefault(node.func.id, []).append(node)
        names = ("find_yield_data", "find_yield_simc", "finalize_full_background_subtraction_e8_2", "build_correction_ledger")
        for name in names:
            self.assertEqual(len(calls[name]), 1, name)
        self.assertEqual(sorted(calls[name][0].lineno for name in names), [calls[name][0].lineno for name in names])
        for name in names[:2]:
            self.assertEqual([ast.unparse(arg) for arg in calls[name][0].args], ["histlist", "inpDict"])
        finalizer = calls[names[2]][0]
        loops = [node for node in tree.body if isinstance(node, ast.For) and any(child is finalizer for child in ast.walk(node))]
        self.assertEqual(len(loops), 1)
        self.assertEqual(ast.unparse(loops[0].iter), "histlist")
        self.assertEqual(ast.unparse(loops[0].target), "hist")
        self.assertIsInstance(loops[0].body[0], ast.If)
        self.assertIn('!= \'kaon\'', ast.unparse(loops[0].body[0].test))
        self.assertIsInstance(loops[0].body[0].body[0], ast.Continue)

    def test_main_ast_diff_is_only_scoped_paths_input_order_and_prior_simc_move(self):
        before = subprocess.check_output(["git", "show", "HEAD:src/main.py"], cwd=ROOT).decode()
        after = (ROOT / "src/main.py").read_text(encoding="utf-8")
        block = 'stage_start = perf_counter()\nyieldDict.update(find_yield_simc(histlist, inpDict))\nrecord_stage_time("Step 6 simc yields", stage_start)\n'
        self.assertEqual(before.count(block), 1)
        self.assertEqual(after.count(block), 1)
        self.assertEqual(canonical_main(before), canonical_main(after))


    def test_main_ast_guard_still_rejects_unrelated_numerical_change(self):
        source = (ROOT / "src/main.py").read_text(encoding="utf-8")
        self.assertEqual(source.count("find_yield_data(histlist, inpDict)"), 1)
        changed = source.replace("find_yield_data(histlist, inpDict)", "find_yield_data([], inpDict)")
        self.assertNotEqual(canonical_main(source), canonical_main(changed))


if __name__ == "__main__":
    unittest.main()
