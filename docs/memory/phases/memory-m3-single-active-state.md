# Memory M3 schema-3 activation and single active-state authority

## Status

`SOURCE REVIEWED` — M3 is a local representation and ownership migration, not
farm/runtime validation.

## Starting identity and purpose

- Starting branch: `test`.
- Starting HEAD: `92289ff1dadc294eb42ac52b5e716be67ccadf98`.
- Starting worktree: clean.

M3 activates the reviewed schema-3 representation, makes `CURRENT.md` the
sole ordinary active-state authority, and preserves the scientific/runtime
frontier without changing calculation, production behavior, evidence, or farm
status.

## Changed files

- `docs/memory/CURRENT.md`
- `docs/memory/handoffs/CURRENT_HANDOFF.md`
- `docs/memory/roadmap/CURRENT.md`
- `docs/memory/AGENTS.md`
- `docs/memory/README.md`
- `docs/memory/MAINTENANCE.md`
- `tools/check_memory_health.py`
- `testing/test_memory_health.py`
- `docs/memory/phases/memory-m3-single-active-state.md`
- regenerated `docs/memory/manifest.json`

## Ownership result

CURRENT now has minimal schema-3 frontmatter and is the sole ordinary
active-state representation. The handoff has no frontmatter and its normal
Transfer State is exactly `No exceptional transfer state is recorded.` The
roadmap retains Phase-F dependency/status structure but no active-state
frontmatter or exact active action.

## Enforcement and transitional tools

Memory health now checks that a schema-3 roadmap exists but has no frontmatter,
`## NEXT`, or active `NEXT —` statement. Deterministic unit coverage exercises
the valid target and all three prohibited roadmap forms. The manifest remains
schema version 2 with minimal schema-3 active state; the unchanged bootstrap
reports that transitional summary and does not replace reading CURRENT.

## Preserved frontier and validation

F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`; E.8.1 remains
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; F.6.3 and F.6.4 remain
blocked in their existing dependency order. The accepted yield remains the
baseline, Method B remains diagnostic only, and Method A remains detached and
non-production.

Using the discovered local Python 3.12 interpreter with `-B`, the M2 baseline
health, manifest, bootstrap, 33-test unit suite, and three embedded self-tests
passed before editing, all with return code 0:

- `<PYTHON> -B tools/check_memory_health.py --root .`
- `<PYTHON> -B tools/update_memory_manifest.py --root . --check`
- `<PYTHON> -B tools/memory_bootstrap.py --root . --json`
- `<PYTHON> -B -m unittest testing.test_memory_health -v`
- the three corresponding `--self-test` commands.

After the coherent migration, all final commands returned 0:

- `<PYTHON> -B -m unittest testing.test_memory_health -v` — 36 tests passed.
- Each of the three `--self-test` commands — `SELF-TEST: PASS`.
- `<PYTHON> -B tools/update_memory_manifest.py --root . --write` and the
  subsequent `--check` — manifest written and `MANIFEST: PASS`.
- `<PYTHON> -B tools/check_memory_health.py --root .` — `MEMORY HEALTH: PASS`.
- `<PYTHON> -B tools/memory_bootstrap.py --root . --json` — schema-3 summary
  and passing health.

No farm command, ROOT/PyROOT or full-main validation, commit, push,
lifecycle-hook work, science change, or runtime claim occurred.

## NEXT

NEXT — Memory M4: normalize MEMORY.md into durable cross-phase knowledge only,
moving detailed chronology/current-gate material to verified canonical records
without losing unique scientific provenance.
