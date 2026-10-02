"""Narrow AST audit of the Step-6 presentation consumer ordering."""
import ast
from pathlib import Path
import subprocess
import unittest

ROOT = Path(__file__).resolve().parents[1]


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

    def test_main_ast_diff_is_only_moving_existing_simc_block(self):
        before = subprocess.check_output(["git", "show", "HEAD:src/main.py"], cwd=ROOT).decode()
        after = (ROOT / "src/main.py").read_text(encoding="utf-8")
        block = 'stage_start = perf_counter()\nyieldDict.update(find_yield_simc(histlist, inpDict))\nrecord_stage_time("Step 6 simc yields", stage_start)\n'
        self.assertEqual(before.count(block), 1)
        self.assertEqual(after.count(block), 1)
        self.assertEqual(ast.dump(ast.parse(before.replace(block, ""))), ast.dump(ast.parse(after.replace(block, ""))))


if __name__ == "__main__":
    unittest.main()
