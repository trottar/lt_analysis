# Memory M8 — manifest and bootstrap identity semantics

## Status

`SOURCE REVIEWED` after the recorded deterministic local acceptance checks.
This is memory tooling/source work only, not farm, ROOT/PyROOT, or full
analysis validation.

## Starting identity

M8 began on branch `test` at
`eb1086c38d93a1f284123213025da9a247eae659` with a clean worktree
(`## test...origin/test`).

## Scope and ownership inventory

M8 removes manifest semantic state, moves strict schema-3 CURRENT frontmatter
validation into health, redesigns bootstrap as a dynamic startup observation,
and advances CURRENT to M9. It does not change science or runtime behavior.

Before M8, `update_memory_manifest.parse_active_state()` owned CURRENT parsing
and was called by health and bootstrap. The version-2 manifest held
`active_state`, `files`, `generated_date_utc`, and `observed_git_head`, so date,
stored HEAD, and copied active state could make it stale. Bootstrap emitted
three `active_records`, parsed `active_state`, queried only branch/HEAD, treated
empty Git output as unavailable through a generic helper, and reported no
exceptional-handoff presence. The M7 tests encoded the transitional parser,
manifest envelope, and bootstrap shape.

## Before/after manifest contract

Before M8 the manifest used schema version 2 with `active_state`, `files`,
`generated_date_utc`, and `observed_git_head`. After M8 it uses schema version
3 and exactly `files` plus `schema_version`. It remains path-sorted file
integrity metadata, excludes itself, and has no active state, timestamp, or
repository identity field.

## Bootstrap contract

Bootstrap now reports dynamic Git branch/HEAD/worktree state, the ordered
five-core records and byte sizes, direct CURRENT references, memory-health
result, and whether the handoff transfer state is exceptional. It reports no
active state or scientific/current-state synthesis; humans still read the core.

## Validation

The pushed-M7 baseline health, manifest check, bootstrap JSON, 32-test suite,
and all three embedded self-tests each returned 0. After M8, the strict suite
returned 0 with 35 tests; health, manifest, and bootstrap self-tests returned
0. The completed commands were `<PYTHON> -B -m unittest
testing.test_memory_health -v`, `<PYTHON> -B tools/check_memory_health.py
--self-test`, `<PYTHON> -B tools/update_memory_manifest.py --self-test`,
`<PYTHON> -B tools/memory_bootstrap.py --self-test`, `<PYTHON> -B
tools/update_memory_manifest.py --root . --write`, `<PYTHON> -B
tools/update_memory_manifest.py --root . --check`, `<PYTHON> -B
tools/check_memory_health.py --root .`, and `<PYTHON> -B
tools/memory_bootstrap.py --root . --json`; each returned 0. Live bootstrap
reported dynamic branch/HEAD, a dirty M8 working tree with preserved
short-status lines, five core records, resolved CURRENT references, health
pass, and `exceptional_handoff.present = false`. The final JSON shape audit,
diff check, and exact allowlist audit are recorded after final regeneration.

## Preserved scientific frontier

F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`; E.8.1 remains
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; F.6.3 remains `BLOCKED`
pending E.8.1; and F.6.4 remains `BLOCKED` pending F.6.3 evidence. The
accepted yield remains baseline, Method B remains diagnostic/cross-check only,
and Method A remains detached/non-production. No farm command ran.

## Next

NEXT — Memory M9: fresh-session integration audit.
