"""Deterministic schema-2/schema-3 memory-health compatibility coverage."""

from __future__ import annotations

import importlib
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = REPO_ROOT / "tools"
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

health = importlib.import_module("check_memory_health")
manifest_tool = importlib.import_module("update_memory_manifest")


class MemoryHealthCompatibilityTests(unittest.TestCase):
    def assert_error(self, errors: list[str], fragment: str) -> None:
        self.assertTrue(any(fragment in error for error in errors), errors)

    def run_checks(self, root: Path, *, resolver=None) -> list[str]:
        errors, warnings = health.run_checks(
            root,
            health.SOFT_LIMIT_BYTES,
            health.HARD_LIMIT_BYTES,
            check_manifest=False,
            commit_resolver=resolver or (lambda _root, _commit: None),
        )
        self.assertEqual(warnings, [])
        return errors

    def write_schema2(self, root: Path) -> None:
        health.write_fixture(root)

    def write_schema3(self, root: Path) -> None:
        health.write_schema3_fixture(root)

    def test_schema2_healthy_fixture_passes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema2(root)
            self.assertEqual(self.run_checks(root), [])

    def test_schema2_mirrored_active_state_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema2(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8").replace("Fixture work item", "other"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "active-state metadata differs")

    def test_schema2_malformed_frontmatter_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema2(root)
            (root / "docs/memory/CURRENT.md").write_text("# no frontmatter\n", encoding="utf-8")
            self.assert_error(self.run_checks(root), "invalid active-state frontmatter")

    def test_schema2_missing_commit_field_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema2(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(
                health.fixture_current().replace("source_commit: bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb\n", ""),
                encoding="utf-8",
            )
            self.assert_error(self.run_checks(root), "invalid active-state frontmatter")

    def test_schema2_invalid_commit_syntax_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema2(root)
            current = root / "docs/memory/CURRENT.md"
            replacement = health.fixture_current().replace("baseline_commit: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "baseline_commit: invalid")
            current.write_text(replacement, encoding="utf-8")
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(health.fixture_frontmatter().replace("baseline_commit: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "baseline_commit: invalid") + "# Handoff\n", encoding="utf-8")
            roadmap = root / "docs/memory/roadmap/CURRENT.md"
            roadmap.write_text(health.fixture_frontmatter().replace("baseline_commit: aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", "baseline_commit: invalid") + "# Roadmap\n", encoding="utf-8")
            self.assert_error(self.run_checks(root), "not a 40-character lowercase SHA")

    def test_schema2_unresolved_commit_fails_when_resolver_reports_false(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema2(root)
            self.assert_error(self.run_checks(root, resolver=lambda _root, _commit: False), "does not resolve")

    def test_schema2_duplicate_or_out_of_order_current_heading_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema2(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8") + "\n## Next Action\n\nduplicate\n", encoding="utf-8")
            self.assert_error(self.run_checks(root), "exactly the ordered")

    def test_schema2_runtime_claim_without_evidence_link_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema2(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(health.fixture_current().replace("[evidence](evidence/accepted.md)", "no evidence link", 1), encoding="utf-8")
            self.assert_error(self.run_checks(root), "runtime-closed claim lacks")

    def test_schema3_valid_target_fixture_passes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            self.assertEqual(self.run_checks(root), [])

    def test_schema3_extra_frontmatter_field_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("memory_schema: 3", "memory_schema: 3\nactive_status: forbidden"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "schema-3 frontmatter fields differ")

    def test_schema3_missing_schema_marker_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("memory_schema: 3\n", ""), encoding="utf-8")
            self.assert_error(self.run_checks(root), "invalid active-state frontmatter")

    def test_schema3_unsupported_schema_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("memory_schema: 3", "memory_schema: 4"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "memory_schema must be 2 or 3")

    def test_schema3_duplicate_heading_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8") + "\n## Next Action\n\nsecond\n", encoding="utf-8")
            self.assert_error(self.run_checks(root), "exactly the ordered")

    def test_schema3_out_of_order_heading_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            text = current.read_text(encoding="utf-8")
            text = text.replace("## Active Objective", "## Temporary", 1).replace("## Current Work Item", "## Active Objective", 1).replace("## Temporary", "## Current Work Item", 1)
            current.write_text(text, encoding="utf-8")
            self.assert_error(self.run_checks(root), "exactly the ordered")

    def test_schema3_empty_required_section_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("## Blockers\n\nfixture", "## Blockers\n"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "section is missing or empty: Blockers")

    def test_schema3_zero_next_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("NEXT — Complete the one fixture action.", "fixture"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "exactly one active NEXT")

    def test_schema3_two_next_statements_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("fixture", "NEXT — another action", 1), encoding="utf-8")
            self.assert_error(self.run_checks(root), "exactly one active NEXT")

    def test_schema3_next_outside_next_action_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            text = current.read_text(encoding="utf-8").replace("NEXT — Complete the one fixture action.", "fixture")
            current.write_text(text.replace("## Active Objective\n\nfixture", "## Active Objective\n\nNEXT — misplaced"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "must occur in ## Next Action")

    def test_schema3_runtime_claim_without_evidence_link_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("[evidence](evidence/accepted.md)", "no evidence link", 1), encoding="utf-8")
            self.assert_error(self.run_checks(root), "runtime-closed claim lacks")

    def test_schema3_empty_transfer_state_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8").replace("No exceptional transfer state is recorded.\n\n", ""), encoding="utf-8")
            self.assert_error(self.run_checks(root), "Transfer State must not be empty")

    def test_schema3_handoff_requires_current_sole_authority(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8").replace("CURRENT.md is the sole authoritative resumable state. ", ""), encoding="utf-8")
            self.assert_error(self.run_checks(root), "must state CURRENT.md sole authority")

    def test_schema3_handoff_frontmatter_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text("---\nmemory_schema: 3\n---\n" + handoff.read_text(encoding="utf-8"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "must not contain active-state frontmatter")

    def test_schema3_handoff_requires_non_override_statement(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8").replace("The handoff cannot override CURRENT.md. ", ""), encoding="utf-8")
            self.assert_error(self.run_checks(root), "cannot override CURRENT.md")

    def test_schema3_handoff_resume_must_point_to_current(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            text = handoff.read_text(encoding="utf-8")
            resume = text.index("## Resume")
            handoff.write_text(text[:resume] + text[resume:].replace("CURRENT.md", "the state record"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "Resume must point to CURRENT.md")

    def test_schema3_handoff_routine_current_heading_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8").replace("## Transfer State", "## Active Objective"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "must not contain routine CURRENT headings")

    def test_schema3_bootstrap_order_drift_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            agents = root / "docs/memory/AGENTS.md"
            agents.write_text(agents.read_text(encoding="utf-8").replace("1. `AGENTS.md`\n2. `CURRENT.md`", "1. `CURRENT.md`\n2. `AGENTS.md`"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "bootstrap order differs")

    def test_schema3_bootstrap_requires_full_core_read(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            agents = root / "docs/memory/AGENTS.md"
            agents.write_text(agents.read_text(encoding="utf-8").replace("Read these five files in full, in this exact order:", "Read these files:"), encoding="utf-8")
            self.assert_error(self.run_checks(root), "full five-file core read")

    def test_schema3_bootstrap_prohibits_eager_expansion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            agents = root / "docs/memory/AGENTS.md"
            agents.write_text(agents.read_text(encoding="utf-8").replace("Do not eagerly load the whole memory hierarchy.", "Expand without restriction."), encoding="utf-8")
            self.assert_error(self.run_checks(root), "prohibit eager/unbounded expansion")

    def test_manifest_parser_preserves_schema2_result(self):
        with tempfile.TemporaryDirectory() as directory:
            current = Path(directory) / "CURRENT.md"
            current.write_text(manifest_tool.fixture_current(), encoding="utf-8")
            parsed = manifest_tool.parse_active_state(current)
            self.assertEqual(parsed["memory_schema"], 2)
            self.assertEqual(set(parsed), set(manifest_tool.ACTIVE_STATE_KEYS))

    def test_manifest_parser_accepts_minimal_schema3(self):
        with tempfile.TemporaryDirectory() as directory:
            current = Path(directory) / "CURRENT.md"
            current.write_text("---\nmemory_schema: 3\n---\n# Current\n", encoding="utf-8")
            self.assertEqual(manifest_tool.parse_active_state(current), {"memory_schema": 3})

    def test_manifest_parser_rejects_extra_schema3_field(self):
        with tempfile.TemporaryDirectory() as directory:
            current = Path(directory) / "CURRENT.md"
            current.write_text("---\nmemory_schema: 3\nactive_status: forbidden\n---\n", encoding="utf-8")
            with self.assertRaisesRegex(manifest_tool.ActiveStateError, "schema-3 frontmatter fields differ"):
                manifest_tool.parse_active_state(current)

    def test_manifest_parser_rejects_unsupported_and_malformed_frontmatter(self):
        with tempfile.TemporaryDirectory() as directory:
            current = Path(directory) / "CURRENT.md"
            current.write_text("---\nmemory_schema: 4\n---\n", encoding="utf-8")
            with self.assertRaisesRegex(manifest_tool.ActiveStateError, "memory_schema must be 2 or 3"):
                manifest_tool.parse_active_state(current)
            current.write_text("---\nmemory_schema: 3\n", encoding="utf-8")
            with self.assertRaisesRegex(manifest_tool.ActiveStateError, "frontmatter is not closed"):
                manifest_tool.parse_active_state(current)

    def test_transitional_manifest_represents_schema3_with_schema_version_two(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_schema3(root)
            manifest_tool.write_manifest(root)
            payload = json.loads((root / "docs/memory/manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(payload["schema_version"], 2)
            self.assertEqual(payload["active_state"], {"memory_schema": 3})
            self.assertEqual(manifest_tool.check_manifest(root), [])


if __name__ == "__main__":
    unittest.main()
