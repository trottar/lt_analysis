"""Declarative inventory and pushed-source binding for the narrow plot gate."""
import unittest
from pathlib import Path
from testing import run_e8_4_fix5_left_lowe_plot_gate as owner


class ProfileTests(unittest.TestCase):
    def test_canonical_five_inventory_and_empty_source_range_allowlist(self):
        profile = owner.collector.load_validation_profile(Path(__file__).resolve().parents[1] / owner.PROFILE)
        self.assertEqual(profile["settings"], [
            {"phi": "Left", "epsilon": "lowe"}, {"phi": "Left", "epsilon": "highe"},
            {"phi": "Center", "epsilon": "lowe"}, {"phi": "Center", "epsilon": "highe"},
            {"phi": "Right", "epsilon": "highe"}])
        self.assertEqual(owner.collector.resolve_settings("Left", "lowe", profile), (("Left", "lowe"),))
        self.assertEqual(profile["collection_mode"], "generic_artifacts")
        names = {declaration["basename_template"].format(kinematic=owner.KINEMATIC, phi="Left", epsilon="lowe")
                 for scope in ("global", "settings") for declaration in profile["artifacts"][scope]}
        self.assertEqual(names, set(owner.CANDIDATES) | {owner.PDF, owner.PAGE_MANIFEST, owner.SUMMARY,
            "kaon_FullAnalysis_Q4p4W2p74_lowe.json",
            "kaon_FullAnalysis_Q4p4W2p74_lowe_correction_ledger_no_empirical_residual.json",
            "kaon_FullAnalysis_Q4p4W2p74_lowe_correction_ledger_no_empirical_residual.csv"})
        self.assertTrue(all(entry["required"] for scope in ("global", "settings") for entry in profile["artifacts"][scope]))
        self.assertEqual(profile["source_identity"], {"required_analysis_commit": owner.BASE_HEAD,
            "allowed_committed_files": [], "allowed_non_analysis_path_prefixes": ["docs/memory/"]})

    def test_effective_profile_binds_only_exact_source_commit(self):
        repo = Path(__file__).resolve().parents[1]
        template_bytes = (repo / owner.PROFILE).read_bytes()
        before = owner.collector.load_validation_profile(repo / owner.PROFILE)
        resolved = owner.resolved_profile(repo, "a" * 40)
        self.assertEqual(resolved["source_identity"]["required_analysis_commit"], "a" * 40)
        resolved["source_identity"]["required_analysis_commit"] = owner.BASE_HEAD
        self.assertEqual(resolved, before)
        self.assertEqual((repo / owner.PROFILE).read_bytes(), template_bytes)


if __name__ == "__main__":
    unittest.main()
