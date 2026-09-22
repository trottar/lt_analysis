# Memory M5 Linux/JLab operational-role separation

## Status

`SOURCE REVIEWED` — M5 is operational-memory normalization only. It is not
scientific, farm, ROOT/PyROOT, or runtime validation.

## Starting identity

- Starting branch: `test`.
- Starting HEAD: `186998d2dd31c1269c4eff81baa0180f58795bcb`.
- Starting worktree: clean.

## Scope

M5 separates operational record ownership and advances CURRENT to M6. It
changes no scientific/runtime state, analysis behavior, evidence, farm
procedure, profile, collector, or source.

## Ownership map

| Subject | Canonical owner | Pointer/reference surfaces | Duplicate procedure removed/avoided | Result |
| --- | --- | --- | --- | --- |
| Startup order | AGENTS | README; MAINTENANCE | Kept only as intentional health-critical repetition | PASS |
| Git identity commands | TOOLS | AGENTS startup requirement | No command blocks elsewhere | PASS |
| Memory-control commands | TOOLS | MAINTENANCE interpreter/sequence pointer | No maintenance command catalog | PASS |
| JLab shell | TOOLS | COMMUNICATION requires valid `tcsh` | No competing environment definition | PASS |
| Farm paths | TOOLS | COMMUNICATION link | Removed duplicated literal path definitions | PASS |
| Detailed farm packaging | farm-validation-bundle decision | TOOLS; COMMUNICATION links | No phase command or detached-worktree sequence copied | PASS |
| Farm request/return communication | COMMUNICATION | TOOLS and decision links | Kept interaction pattern out of command catalog | PASS |
| Presentation-only versus bundle-only | COMMUNICATION | Detailed decision link | Kept request-level distinction without procedure copy | PASS |
| Source-changing Codex workflow | CODEX | AGENTS high-level authority; README link | No task-class procedure in AGENTS | PASS |
| Actor sequence | CODEX | AGENTS high-level authority | No competing workflow policy | PASS |
| Contract field list | CODEX | Frozen template instantiates policy | No full list in other operational records | PASS |
| Memory-maintenance sequence | MAINTENANCE | TOOLS command reference | No operational decision logic in TOOLS | PASS |
| User delivery preferences | USER | README navigation | No delivery-policy duplication | PASS |
| Lifecycle-hook status | M0-M9 plan / CURRENT | No operational dependency | Remains `BLOCKED / DEFERRED`; no hook work added | PASS |

## File-role results

- **AGENTS:** high-level evidence, startup, scientific boundaries, execution
  authority, and pointers to specialized operational owners.
- **USER:** stable Hall C/Linux-JLab context and collaboration/delivery
  preferences only.
- **TOOLS:** generic Linux/POSIX repository commands, memory commands, JLab
  shell/path facts, and provenance inspection references.
- **COMMUNICATION:** concise farm request/return interaction and requested
  operation distinctions.
- **CODEX:** source-changing task classes, contract policy, actor sequence, and
  local-versus-farm boundary.
- **MAINTENANCE:** memory roles, semantic/size triggers, startup health, and
  maintenance sequence.
- **README:** startup overview and navigation only.

## Linux/JLab result

No PowerShell, `.ffpkg`, Windows transport workflow, or durable local Windows
Python path was introduced. TOOLS owns the JLab `tcsh` and canonical path facts;
the detailed farm package procedure remains in its decision record.

## Validation

Using the discovered local Python 3.12 interpreter represented as `<PYTHON>`
with `-B`, the M4 baseline passed before editing: memory health, manifest
check, schema-3 bootstrap, 36 unit tests, and all three embedded self-tests.
An initial post-edit live health run exposed required README navigation and
concise authority-marker compatibility text; the affected owned records were
corrected without modifying a frozen tool or test.

Final M5 commands and results:

- `<PYTHON> -B -m unittest testing.test_memory_health -v` — 36 tests passed.
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

NEXT — Memory M6: roadmap/source/history normalization.
