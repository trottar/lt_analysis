"""Non-ROOT checks for the shared Lambda-gate regression-table writer."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src" / "utility"))

from lambda_gate_regression import upsert_sorted_csv  # noqa: E402


class LambdaGateRegressionCsvTests(unittest.TestCase):
    def test_upsert_is_keyed_sorted_and_does_not_expose_temporary_csv(self):
        fields = ("kinematic", "epsilon", "phi_setting", "status")
        key_fields = ("kinematic", "epsilon", "phi_setting")
        with tempfile.TemporaryDirectory() as directory:
            path = str(Path(directory) / "proton_cleaning_lambda_gate_regression_summary.csv")
            upsert_sorted_csv(path, {
                "kinematic": "Q4p4W2p74", "epsilon": "high", "phi_setting": "Right", "status": "fail",
            }, fields, key_fields)
            upsert_sorted_csv(path, {
                "kinematic": "Q4p4W2p74", "epsilon": "low", "phi_setting": "Center", "status": "pass",
            }, fields, key_fields)
            upsert_sorted_csv(path, {
                "kinematic": "Q4p4W2p74", "epsilon": "high", "phi_setting": "Right", "status": "pass",
            }, fields, key_fields)

            with open(path, newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(
                [(row["epsilon"], row["phi_setting"], row["status"]) for row in rows],
                [("high", "Right", "pass"), ("low", "Center", "pass")],
            )
            self.assertFalse(list(Path(directory).glob("*.tmp")))


if __name__ == "__main__":
    unittest.main()
