# Memory M6 roadmap/source/history normalization

## Status

`SOURCE REVIEWED` — M6 is repository-memory/navigation normalization, not
scientific, farm, ROOT/PyROOT, or runtime validation.

## Starting identity

- Starting branch: `test`.
- Starting HEAD: `799416599671ff0bdc3c7d4024fbd1f7f2f456cc`.
- Starting worktree: clean.

## Scope

M6 renames the live schema-3 roadmap to STATUS, de-authoritizes
roadmap/source/history navigation, adds an immutable-import wrapper, narrows
the schema-3 checker/fixture logic, and advances CURRENT to M7. It changes no
scientific/runtime status, evidence, analysis source, collector, profile, or
production behavior.

## Pre-edit ambiguity inventory

| Match class | Classification | M6 result |
| --- | --- | --- |
| Live roadmap pointers in MEMORY, README, and MAINTENANCE | Active/control ambiguity | Updated to `roadmap/STATUS.md` |
| SOURCE_INDEX live identity heading and `b02316f83bf8ef18641fa7217b90c04bcdca10e3` | Active/control ambiguity | Replaced with dynamic-identity rule and historical observation label |
| PROJECT_HISTORY current-checkpoint wording and `b02316f83bf8ef18641fa7217b90c04bcdca10e3` | Legitimate history needing relabel | Retained as a historical 2026-09-11 repository observation |
| CURRENT and schema-v3 decision live-HEAD prohibitions | Active authority/contract semantics | Retained; neither asserts a stored live identity |
| M0 evidence and M3/M4 phase references | Frozen historical text | Left unchanged |
| Immutable import package matches | Immutable imported provenance | Left byte-for-byte unchanged |
| Generated manifest entry | Generated integrity metadata | Regenerated only after authored changes |

## Roadmap result

- Old live schema-3 path: `docs/memory/roadmap/CURRENT.md`.
- New live schema-3 path: `docs/memory/roadmap/STATUS.md`.
- Phase-F dependency/status semantics and role-labeled identities are
  preserved; STATUS has no frontmatter, `## NEXT`, or `NEXT —`.
- Schema-2 fixture compatibility retains `roadmap/CURRENT.md` until M7.

## Source/history audit

| Surface | Previous ambiguity | Final role | Fixed-SHA treatment | Active-authority result |
| --- | --- | --- | --- | --- |
| SOURCE_INDEX | Stored live migration head | Provenance/navigation | Full SHAs labeled historical reviewed-source, import-preparation, or migration-era observations | No active status or stored live HEAD |
| ARTIFACT_INDEX | Snapshot limitations could read as current | Historical 2026-09-11 artifact snapshot | Historical evaluated identities remain table facts | CURRENT/evidence own current acceptance |
| chats/CHAT_INDEX | Chat summaries could appear authoritative | Distilled historical navigation | No full SHA added | CURRENT/source/evidence/phase/decision outrank it |
| history/CHAT_INDEX | Reconstructed evidence tags could appear active | Historical 2026 era mapping | No stored live HEAD | Tags are historical classification labels |
| history/PROJECT_HISTORY | Current-checkpoint/live wording | Chronology through import checkpoint | Full SHA explicitly historical repository observation; F.1.Fix.5 SHA remains historical analysis-source role | Dynamic Git/CURRENT/evidence own present state |
| import README/package | Imported inputs lacked wrapper boundary | Immutable historical provenance navigation | No imported SHA relabeled or changed | Package is never active authority |

## Imported package integrity

Every pre-existing file under
`docs/memory/import/chatgpt_project_history_2026-09-11/` remained
byte-for-byte unchanged. `git diff --name-only HEAD --
docs/memory/import/chatgpt_project_history_2026-09-11` produced no output.

## Validation

Using the discovered local Python 3.12 interpreter represented as `<PYTHON>`
with `-B`, all pre-change M5 baseline commands returned 0: memory health,
manifest check, schema-3 bootstrap, 36 unit tests, and all three embedded
self-tests.

Final M6 commands and results:

- `<PYTHON> -B -m unittest testing.test_memory_health -v` — 38 tests passed.
- `<PYTHON> -B tools/check_memory_health.py --self-test` — `SELF-TEST: PASS`.
- `<PYTHON> -B tools/update_memory_manifest.py --self-test` — `SELF-TEST: PASS`.
- `<PYTHON> -B tools/memory_bootstrap.py --self-test` — `SELF-TEST: PASS`.
- `<PYTHON> -B tools/update_memory_manifest.py --root . --write` then
  `--check` — manifest regenerated and `MANIFEST: PASS`.
- `<PYTHON> -B tools/check_memory_health.py --root .` — `MEMORY HEALTH: PASS`.
- `<PYTHON> -B tools/memory_bootstrap.py --root . --json` — schema-3 summary
  with passing health.
- `git diff --check` — passed.

## Preserved scientific frontier

F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`; E.8.1 remains
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; F.6.3 remains `BLOCKED`
pending E.8.1; and F.6.4 remains `BLOCKED` pending F.6.3 evidence. The
accepted yield remains the baseline, Method B remains diagnostic/cross-check
only, and Method A remains detached and non-production. No farm command ran.

## NEXT

NEXT — Memory M7: strict semantic health enforcement.
