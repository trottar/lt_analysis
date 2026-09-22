# Memory M1 schema-v3 compatibility

## Status

`SOURCE REVIEWED` — Memory M1 adds preparatory dual-schema memory validation.
It is local pure-Python infrastructure only, not farm/runtime validation.

## Starting identity and scope

- Starting branch: `test`.
- Starting HEAD: `50858f2b6a85f1c3355569acef555fadfae93e19`.
- Starting worktree: clean.
- Working interpreter:
  `C:\Users\trott\AppData\Local\Programs\Python\Python312\python.exe -B`.

M1 defines schema-3 semantics and validates synthetic schema-3 fixtures while
the live repository remains schema 2.  CURRENT, MEMORY, CURRENT_HANDOFF, and
roadmap/CURRENT were not migrated.  `tools/memory_bootstrap.py` was not
modified.

## Changed files

- `docs/memory/decisions/memory-schema-v3-contract.md`
- `docs/memory/phases/memory-m1-schema-v3-compatibility.md`
- `tools/check_memory_health.py`
- `tools/update_memory_manifest.py`
- `testing/test_memory_health.py`
- regenerated `docs/memory/manifest.json`

## Compatibility results

Schema 2 retains its required frontmatter, three-surface equality, SHA syntax
and resolver behavior, ordered CURRENT headings, links, runtime-evidence rule,
authority checks, size limits, role uniqueness, and manifest-freshness checks.
The live schema-2 memory-health and unchanged bootstrap checks pass.

Schema-3 fixtures validate minimal CURRENT frontmatter, ordered/non-empty
CURRENT sections, one NEXT inside Next Action, runtime evidence links,
exceptional-only handoff structure/authority/resume semantics, and the future
five-file full-read/selective-expansion contract.  Negative fixtures exercise
extra or unsupported frontmatter, malformed state, headings, NEXT placement,
runtime links, handoff role violations, and bootstrap-contract drift.

The transitional manifest parser accepts only minimal schema-3 frontmatter and
returns `{"memory_schema": 3}`.  Manifest output remains schema version 2 and
retains the current manifest fields; this is not M8 redesign.

## Local validation

All commands used the stated interpreter with `-B`:

- `-m unittest testing.test_memory_health -v` — return 0, PASS, 33 tests.
- `tools/check_memory_health.py --self-test` — return 0, PASS.
- `tools/update_memory_manifest.py --self-test` — return 0, PASS.
- `tools/check_memory_health.py --root .` — return 0, PASS.
- `tools/update_memory_manifest.py --root . --check` — return 0, PASS.
- `tools/memory_bootstrap.py --self-test` — return 0, PASS.
- `tools/memory_bootstrap.py --root . --json` — return 0, PASS.
- `git diff --check` — return 0, PASS.

## Preserved boundaries

No science source, validation profile, collector, active memory record, runtime
evidence record, procedure-PDF source, farm artifact, lifecycle hook, or
Windows packaging workflow changed.  No farm command, commit, push, ROOT/PyROOT
validation, full `main.py` validation, or farm/runtime status upgrade occurred.

## NEXT

NEXT — Memory M2: five-file bootstrap ownership and creation of USER.md,
TOOLS.md, and LEARNINGS.md.
