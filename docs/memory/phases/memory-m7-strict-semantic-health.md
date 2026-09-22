# Memory M7 — strict semantic health enforcement

## Status

`SOURCE REVIEWED` after the recorded local deterministic acceptance checks.
This is memory-control/source work only; it is not farm, ROOT/PyROOT, or full
analysis validation.

## Starting identity and scope

M7 began on branch `test` at
`c6b7f758a13060a650f036d25758d1ee8049c66f` with a clean worktree
(`## test...origin/test`). It changes only the M7 memory-control allowlist and
does not alter scientific source, validation profiles/collectors, evidence, or
runtime artifacts.

M7 is the schema-2 compatibility cutoff. Live health and the active-state
parser now accept only minimal schema-3 frontmatter. The manifest intentionally
retains its transitional version-2 envelope and `{"memory_schema": 3}` active
state until M8.

## Pre-edit strictness inventory

Before M7, `check_memory_health.py` retained schema-2 `ACTIVE_SURFACES`,
mirrored-state/commit resolution, schema-2 fixtures, and a generic 16/32 KiB
policy exposed through global size-limit CLI flags. Its README navigation list
covered only nine control links, duplicate-role checking covered seven
top-level roles, local-link checking covered only README/CURRENT, and only
CURRENT runtime-closed claims were checked. It had no deterministic Markdown
representation gate, roadmap runtime-evidence-link rule, or roadmap status
format rule.

`update_memory_manifest.py` retained schema-2 active-state keys and a
successful schema-2 parse/self-test path. `memory_bootstrap.py` retained a
schema-2 self-test fixture and `active_status` assertion. The M6 roadmap had
bare runtime-closed/active/development/blocked labels and runtime-closed blocks
without direct Markdown evidence links. No non-CURRENT live/control/navigation
surface contained `NEXT —` in the pre-edit inventory.

## Enforced M7 invariants

- Schema 3 is the only accepted live active-state representation; CURRENT is
  the sole ordinary active-state authority and owns exactly one `NEXT —`.
- The exact ordered five-file startup core, selective-expansion constraints,
  exceptional-only handoff semantics, durable MEMORY role, sole STATUS roadmap,
  approved work-state vocabulary, and duplicate-role exclusions are checked.
- Required control/navigation files have UTF-8/NUL/final-newline/one-H1/fence
  balance/local-link deterministic representation checks. Runtime-closed
  CURRENT and STATUS blocks require direct canonical evidence links.
- Final byte thresholds are CURRENT 8/16 KiB, handoff 6/10 KiB, and MEMORY
  30/50 KiB for soft warning/hard failure respectively.

## Local validation

Using the established Python 3.12 interpreter with `-B`, the pushed-M6
baseline passed: health, manifest check, bootstrap JSON, and all three
embedded self-tests returned 0; `testing.test_memory_health -v` returned 0
with 38 tests.

After M7, `testing.test_memory_health -v` returned 0 with 32 strict schema-3
tests. `check_memory_health.py --self-test`,
`update_memory_manifest.py --self-test`, and
`memory_bootstrap.py --self-test` each returned 0. Manifest write/check,
live memory health, and bootstrap JSON each returned 0 after manifest
regeneration. The final diff checks and allowlist audit are recorded in this
same source-review task.

## Preserved scientific/runtime frontier

F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`; E.8.1 remains
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; F.6.3 remains `BLOCKED`
pending E.8.1; and F.6.4 remains `BLOCKED` pending F.6.3 evidence. The
accepted yield remains baseline, Method B remains diagnostic/cross-check only,
and Method A remains detached/non-production. No farm command ran.

## Next

NEXT — Memory M8: manifest and bootstrap identity semantics.
