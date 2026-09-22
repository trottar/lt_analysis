"""Deterministic strict schema-3 memory-health coverage."""

from __future__ import annotations

import importlib
import json
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIR = REPO_ROOT / "tools"
if str(TOOLS_DIR) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIR))

health = importlib.import_module("check_memory_health")
manifest_tool = importlib.import_module("update_memory_manifest")
bootstrap = importlib.import_module("memory_bootstrap")


class StrictMemoryHealthTests(unittest.TestCase):
    def assert_error(self, errors: list[str], fragment: str) -> None:
        self.assertTrue(any(fragment in error for error in errors), errors)

    def write_fixture(self, root: Path) -> None:
        health.write_fixture(root)

    def run_health(self, root: Path) -> tuple[list[str], list[str]]:
        return health.run_checks(root, check_manifest=False)

    def errors(self, root: Path) -> list[str]:
        errors, warnings = self.run_health(root)
        self.assertEqual(warnings, [])
        return errors

    def test_valid_strict_schema3_fixture_passes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            self.assertEqual(self.errors(root), [])

    def test_schema2_current_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("memory_schema: 3", "memory_schema: 2"), encoding="utf-8")
            self.assert_error(self.errors(root), "memory_schema must be 3")

    def test_missing_schema_marker_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("memory_schema: 3\n", ""), encoding="utf-8")
            self.assert_error(self.errors(root), "schema-3 frontmatter fields differ")

    def test_unsupported_schema_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("memory_schema: 3", "memory_schema: 4"), encoding="utf-8")
            self.assert_error(self.errors(root), "memory_schema must be 3")

    def test_extra_frontmatter_field_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("memory_schema: 3", "memory_schema: 3\nactive_status: forbidden"), encoding="utf-8")
            self.assert_error(self.errors(root), "schema-3 frontmatter fields differ")

    def test_health_owned_frontmatter_parser_rejects_malformed_forms(self):
        with tempfile.TemporaryDirectory() as directory:
            current = Path(directory) / "CURRENT.md"
            for payload in (
                "# no frontmatter\n",
                "---\nmemory_schema: 3\n",
                "---\nmemory_schema:\n---\n",
                "---\nmemory_schema: 3\nmemory_schema: 3\n---\n",
            ):
                with self.subTest(payload=payload):
                    current.write_text(payload, encoding="utf-8")
                    with self.assertRaises(health.CurrentFrontmatterError):
                        health.parse_current_schema3(current)

    def test_current_heading_order_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            text = current.read_text(encoding="utf-8")
            text = text.replace("## Active Objective", "## Temporary", 1).replace("## Current Work Item", "## Active Objective", 1).replace("## Temporary", "## Current Work Item", 1)
            current.write_text(text, encoding="utf-8")
            self.assert_error(self.errors(root), "exactly the ordered active-state sections")

    def test_current_empty_section_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("## Blockers\n\nfixture", "## Blockers\n"), encoding="utf-8")
            self.assert_error(self.errors(root), "section is missing or empty: Blockers")

    def test_current_second_next_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("fixture", "NEXT — duplicate", 1), encoding="utf-8")
            self.assert_error(self.errors(root), "exactly one active NEXT")

    def test_current_next_outside_next_action_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            text = current.read_text(encoding="utf-8").replace("NEXT — Complete the one fixture action.", "fixture")
            current.write_text(text.replace("## Active Objective\n\nfixture", "## Active Objective\n\nNEXT — misplaced"), encoding="utf-8")
            self.assert_error(self.errors(root), "must occur only in ## Next Action")

    def test_current_runtime_claim_without_evidence_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("[accepted evidence](evidence/accepted.md)", "accepted evidence"), encoding="utf-8")
            self.assert_error(self.errors(root), "runtime-closed claim lacks")

    def test_current_broken_local_link_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("evidence/accepted.md", "evidence/missing.md"), encoding="utf-8")
            self.assert_error(self.errors(root), "CURRENT.md link target does not exist")

    def test_current_unapproved_status_label_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8").replace("`CLOSED / RUNTIME VALIDATED`", "`BLOCKED / DEFERRED`"), encoding="utf-8")
            self.assert_error(self.errors(root), "unapproved uppercase status label")

    def test_missing_core_file_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            (root / "docs/memory/USER.md").unlink()
            self.assert_error(self.errors(root), "missing required file: docs/memory/USER.md")

    def test_agents_startup_order_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            agents = root / "docs/memory/AGENTS.md"
            agents.write_text(agents.read_text(encoding="utf-8").replace("1. `AGENTS.md`\n2. `CURRENT.md`", "1. `CURRENT.md`\n2. `AGENTS.md`"), encoding="utf-8")
            self.assert_error(self.errors(root), "AGENTS.md five-file startup order")

    def test_agents_incomplete_selective_expansion_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            agents = root / "docs/memory/AGENTS.md"
            agents.write_text(agents.read_text(encoding="utf-8").replace("required canonical evidence/decision/phase records", "other records"), encoding="utf-8")
            self.assert_error(self.errors(root), "AGENTS.md selective-expansion sources")

    def test_readme_startup_order_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            readme = root / "docs/memory/README.md"
            readme.write_text(readme.read_text(encoding="utf-8").replace("1. [AGENTS.md](AGENTS.md)\n2. [CURRENT.md](CURRENT.md)", "1. [CURRENT.md](CURRENT.md)\n2. [AGENTS.md](AGENTS.md)"), encoding="utf-8")
            self.assert_error(self.errors(root), "README.md five-file startup order")

    def test_maintenance_startup_order_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            maintenance = root / "docs/memory/MAINTENANCE.md"
            maintenance.write_text(maintenance.read_text(encoding="utf-8").replace("1. `AGENTS.md`\n2. `CURRENT.md`", "1. `CURRENT.md`\n2. `AGENTS.md`"), encoding="utf-8")
            self.assert_error(self.errors(root), "MAINTENANCE.md five-file startup order")

    def test_handoff_frontmatter_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text("---\nmemory_schema: 3\n---\n" + handoff.read_text(encoding="utf-8"), encoding="utf-8")
            self.assert_error(self.errors(root), "handoff must not contain active-state frontmatter")

    def test_handoff_structure_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8").replace("## Resume", "## Wrong"), encoding="utf-8")
            self.assert_error(self.errors(root), "handoff must contain exactly")

    def test_handoff_stable_sentence_drift_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8").replace("No exceptional transfer state is recorded.", "No exceptional transfer state is recorded today."), encoding="utf-8")
            self.assert_error(self.errors(root), "stable Transfer State")

    def test_handoff_authority_and_non_override_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8").replace("sole authoritative resumable state", "ordinary state").replace("handoff cannot override CURRENT.md", "handoff is supplementary"), encoding="utf-8")
            errors = self.errors(root)
            self.assert_error(errors, "sole authority")
            self.assert_error(errors, "cannot override CURRENT.md")

    def test_handoff_next_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8") + "NEXT — wrong\n", encoding="utf-8")
            self.assert_error(self.errors(root), "handoff must not contain active NEXT")

    def test_memory_active_heading_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            memory = root / "docs/memory/MEMORY.md"
            memory.write_text(memory.read_text(encoding="utf-8").replace("## Authority and evidence", "## Active Objective"), encoding="utf-8")
            self.assert_error(self.errors(root), "durable-knowledge sections")

    def test_memory_next_and_hashes_and_status_fail(self):
        cases = (
            ("\nNEXT — bad\n", "MEMORY.md must not contain active NEXT"),
            ("\naaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n", "40-character"),
            ("\naaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\n", "64-character"),
            ("\nSOURCE REVIEWED\n", "phase-ledger status"),
        )
        for extra, expected in cases:
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.write_fixture(root)
                memory = root / "docs/memory/MEMORY.md"
                memory.write_text(memory.read_text(encoding="utf-8") + extra, encoding="utf-8")
                self.assert_error(self.errors(root), expected)

    def test_roadmap_missing_and_legacy_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            (root / "docs/memory/roadmap/STATUS.md").unlink()
            self.assert_error(self.errors(root), "missing required file: docs/memory/roadmap/STATUS.md")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            (root / "docs/memory/roadmap/CURRENT.md").write_text("# Legacy\n", encoding="utf-8")
            self.assert_error(self.errors(root), "legacy roadmap/CURRENT.md")

    def test_roadmap_frontmatter_and_next_fail(self):
        for text, expected in (
            ("---\nmemory_schema: 3\n---\n# Approved KaonLT roadmap status\n", "roadmap STATUS must not contain active-state frontmatter"),
            ("# Approved KaonLT roadmap status\n\n## NEXT\n\nfixture\n", "roadmap STATUS must not contain ## NEXT"),
            ("# Approved KaonLT roadmap status\n\nNEXT — fixture\n", "roadmap STATUS must not contain active NEXT"),
        ):
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.write_fixture(root)
                (root / "docs/memory/roadmap/STATUS.md").write_text(text, encoding="utf-8")
                self.assert_error(self.errors(root), expected)

    def test_roadmap_runtime_evidence_and_status_presentation_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            roadmap = root / "docs/memory/roadmap/STATUS.md"
            roadmap.write_text("# Approved KaonLT roadmap status\n\n### F\n\n`CLOSED / RUNTIME VALIDATED` — no evidence.\n\nCURRENT.md owns the active objective, blockers, and next action.\n", encoding="utf-8")
            self.assert_error(self.errors(root), "runtime-closed claim lacks")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            roadmap = root / "docs/memory/roadmap/STATUS.md"
            roadmap.write_text("# Approved KaonLT roadmap status\n\n### F\n\n`CLOSED / RUNTIME VALIDATED` — [missing](../evidence/missing.md).\n\nCURRENT.md owns the active objective, blockers, and next action.\n", encoding="utf-8")
            errors = self.errors(root)
            self.assert_error(errors, "link target does not exist")
            self.assert_error(errors, "runtime-closed claim lacks")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            roadmap = root / "docs/memory/roadmap/STATUS.md"
            roadmap.write_text(roadmap.read_text(encoding="utf-8").replace("`ACTIVE`", "ACTIVE"), encoding="utf-8")
            self.assert_error(self.errors(root), "phase-state claim")

    def test_duplicate_roles_fail(self):
        cases = (
            ("docs/memory/CURRENT_copy.md", "duplicate or missing top-level CURRENT role"),
            ("docs/memory/handoffs/CURRENT_HANDOFF_old.md", "duplicate or missing handoff CURRENT_HANDOFF role"),
            ("docs/memory/roadmap/STATUS_v2.md", "duplicate or missing roadmap STATUS role"),
        )
        for relative, expected in cases:
            with self.subTest(relative=relative), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.write_fixture(root)
                path = root / relative
                path.write_text("# Variant\n", encoding="utf-8")
                self.assert_error(self.errors(root), expected)

    def test_readme_navigation_and_other_control_link_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            readme = root / "docs/memory/README.md"
            readme.write_text(readme.read_text(encoding="utf-8").replace("[x](manifest.json)\n", ""), encoding="utf-8")
            self.assert_error(self.errors(root), "README missing required control link: manifest.json")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            tools = root / "docs/memory/TOOLS.md"
            tools.write_text(tools.read_text(encoding="utf-8") + "[missing](missing.md)\n", encoding="utf-8")
            self.assert_error(self.errors(root), "TOOLS.md link target does not exist")

    def test_representation_failures(self):
        cases = (
            ("# Wrong\n", "must contain exactly one H1"),
            ("# KaonLT tools and operational commands", "must end with a newline"),
            ("# KaonLT tools and operational commands\n\n```\n", "unbalanced fenced"),
        )
        for text, expected in cases:
            with self.subTest(expected=expected), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                self.write_fixture(root)
                (root / "docs/memory/TOOLS.md").write_text(text, encoding="utf-8")
                self.assert_error(self.errors(root), expected)

    def test_each_role_size_threshold(self):
        for relative, (soft, hard) in health.SIZE_LIMITS.items():
            with self.subTest(relative=relative, level="at_soft"):
                errors, warnings = health.check_size_limit(relative, soft, soft, hard)
                self.assertEqual((errors, warnings), ([], []))
            with self.subTest(relative=relative, level="warn"):
                errors, warnings = health.check_size_limit(relative, soft + 1, soft, hard)
                self.assertEqual(errors, [])
                self.assertEqual(len(warnings), 1)
            with self.subTest(relative=relative, level="hard"):
                errors, warnings = health.check_size_limit(relative, hard + 1, soft, hard)
                self.assertEqual(warnings, [])
                self.assertEqual(len(errors), 1)

    def test_manifest_integrity_shape_and_drift(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            manifest_tool.write_manifest(root)
            payload = json.loads((root / "docs/memory/manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(set(payload), {"files", "schema_version"})
            self.assertEqual(payload["schema_version"], 3)
            self.assertNotIn("docs/memory/manifest.json", [entry["path"] for entry in payload["files"]])
            self.assertEqual(manifest_tool.check_manifest(root), [])
            current = root / "docs/memory/CURRENT.md"
            current.write_text(current.read_text(encoding="utf-8") + "changed\n", encoding="utf-8")
            self.assertTrue(manifest_tool.check_manifest(root))
            manifest_tool.write_manifest(root)
            extra = root / "docs/memory/EXTRA.md"
            extra.write_text("extra\n", encoding="utf-8")
            self.assertTrue(manifest_tool.check_manifest(root))
            extra.unlink()
            self.assertEqual(manifest_tool.check_manifest(root), [])
            for key, value in (
                ("active_state", {"memory_schema": 3}),
                ("generated_date_utc", "2026-09-22"),
                ("observed_git_head", "a" * 40),
                ("unexpected", True),
            ):
                legacy = dict(manifest_tool.build_manifest(root))
                legacy[key] = value
                (root / "docs/memory/manifest.json").write_bytes(manifest_tool.encoded_manifest(legacy))
                self.assertTrue(any("top-level keys" in problem for problem in manifest_tool.check_manifest(root)))
            unsupported = manifest_tool.build_manifest(root)
            unsupported["schema_version"] = 2
            (root / "docs/memory/manifest.json").write_bytes(manifest_tool.encoded_manifest(unsupported))
            self.assertTrue(any("schema_version" in problem for problem in manifest_tool.check_manifest(root)))

    def test_manifest_does_not_depend_on_current(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "docs/memory/SOURCE.md"
            source.parent.mkdir(parents=True)
            source.write_text("source\n", encoding="utf-8")
            manifest_tool.write_manifest(root)
            self.assertEqual(manifest_tool.check_manifest(root), [])

    def test_bootstrap_core_git_handoff_and_reference_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.write_fixture(root)
            clean_git = {
                "branch": "fixture",
                "head": "a" * 40,
                "worktree": {"available": True, "clean": True, "status_short": []},
            }
            summary = bootstrap.collect_summary(
                root,
                lambda _: (0, "MEMORY HEALTH: PASS"),
                lambda _: clean_git,
            )
            self.assertNotIn("active_state", summary)
            self.assertEqual([item["path"] for item in summary["core_records"]], list(bootstrap.CORE_RECORDS))
            self.assertTrue(all(isinstance(item["bytes"], int) for item in summary["core_records"]))
            self.assertEqual(summary["current_references"], ["docs/memory/evidence/accepted.md"])
            self.assertEqual(summary["exceptional_handoff"], {"present": False})
            self.assertEqual(summary["git"], clean_git)
            self.assertEqual(summary["memory_health"], {"status": "pass", "returncode": 0})

            dirty_git = {
                "branch": "fixture",
                "head": "b" * 40,
                "worktree": {"available": True, "clean": False, "status_short": [" M docs/memory/CURRENT.md", "?? scratch.txt"]},
            }
            dirty = bootstrap.collect_summary(root, lambda _: (0, "MEMORY HEALTH: WARN"), lambda _: dirty_git)
            self.assertEqual(dirty["git"]["worktree"], dirty_git["worktree"])
            self.assertEqual(dirty["memory_health"]["status"], "warn")
            unavailable = {"branch": None, "head": None, "worktree": {"available": False, "clean": None, "status_short": None}}
            self.assertEqual(bootstrap.worktree_label(unavailable["worktree"]), "UNAVAILABLE")
            with mock.patch.object(
                bootstrap,
                "git_command",
                side_effect=[(True, "fixture\n"), (True, "c" * 40 + "\n"), (True, " M docs/memory/CURRENT.md\n?? scratch.txt\n")],
            ):
                observed = bootstrap.observe_git(root)
            self.assertEqual(observed["worktree"], {"available": True, "clean": False, "status_short": [" M docs/memory/CURRENT.md", "?? scratch.txt"]})
            with mock.patch.object(bootstrap, "git_command", return_value=(False, None)):
                unavailable_observed = bootstrap.observe_git(root)
            self.assertEqual(unavailable_observed["worktree"], unavailable["worktree"])

            handoff = root / "docs/memory/handoffs/CURRENT_HANDOFF.md"
            handoff.write_text(handoff.read_text(encoding="utf-8").replace("No exceptional transfer state is recorded.", "Exceptional transfer state."), encoding="utf-8")
            exceptional = bootstrap.collect_summary(root, lambda _: (1, "MEMORY HEALTH: FAIL"), lambda _: clean_git)
            self.assertEqual(exceptional["exceptional_handoff"], {"present": True})
            self.assertEqual(exceptional["memory_health"]["status"], "fail")


if __name__ == "__main__":
    unittest.main()
