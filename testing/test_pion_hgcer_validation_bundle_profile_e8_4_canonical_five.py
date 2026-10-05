"""Compatibility of the isolated owner profile with the unchanged collector."""
from copy import deepcopy
from pathlib import Path
import unittest

from testing import run_e8_4_fix5_canonical_five_plot_gate as owner

ROOT = Path(__file__).resolve().parents[1]


class ProfileTests(unittest.TestCase):
    def test_exact_inventory_and_generic_artifact_keys(self):
        profile = owner.resolved_profile(ROOT, "a" * 40)
        owner.validate_settings(profile["settings"])
        self.assertEqual({(s["phi"], s["epsilon"]) for s in profile["settings"]},
                         {("Left", "lowe"), ("Center", "lowe"), ("Left", "highe"),
                          ("Center", "highe"), ("Right", "highe")})
        self.assertEqual(owner.collector.resolve_settings(profile=profile),
                         tuple((s["phi"], s["epsilon"]) for s in profile["settings"]))
        with self.assertRaises(ValueError):
            owner.collector.resolve_settings("Right", "lowe", profile)
        self.assertEqual(profile["collection_mode"], "generic_artifacts")
        self.assertEqual({e["key"] for e in profile["artifacts"]["global"]},
                         {"candidate_f3", "candidate_f4", "run_summary"})
        self.assertEqual({e["key"] for e in profile["artifacts"]["settings"]},
                         {"procedure_pdf", "page_manifest", "full_analysis",
                          "correction_ledger_json", "correction_ledger_csv"})
        names = owner.artifact_names(profile)
        self.assertIn(owner.SUMMARY, names)
        self.assertNotIn(owner.accepted.SUMMARY, names)
        self.assertEqual(len(names), 19)  # 3 globals + 10 per-phi + 6 per-epsilon
        self.assertTrue(all(e["required"] for scope in profile["artifacts"].values() for e in scope))

    def test_effective_source_pin_does_not_rewrite_template(self):
        raw = (ROOT / owner.PROFILE).read_bytes()
        template = owner.collector.load_validation_profile(ROOT / owner.PROFILE)
        effective = owner.resolved_profile(ROOT, "a" * 40)
        self.assertEqual(effective["source_identity"]["required_analysis_commit"], "a" * 40)
        effective["source_identity"]["required_analysis_commit"] = owner.BASE_HEAD
        self.assertEqual(effective, template)
        self.assertEqual((ROOT / owner.PROFILE).read_bytes(), raw)

    def test_invalid_setting_inventories_fail_closed(self):
        valid = list(owner.CANONICAL_SETTINGS)
        for invalid in (valid[:-1], valid + [valid[0]], valid + [{"phi": "Right", "epsilon": "lowe"}],
                        [{"phi": "unknown", "epsilon": "lowe"}] + valid[1:]):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, "inventory_invalid"):
                owner.validate_settings(deepcopy(invalid))


if __name__ == "__main__":
    unittest.main()
