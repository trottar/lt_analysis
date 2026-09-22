# Memory M4 durable MEMORY normalization

## Status

`SOURCE REVIEWED` — M4 is repository-memory normalization, not scientific or
farm/runtime validation.

## Starting identity and scope

- Starting branch: `test`.
- Starting HEAD: `79ab4c7b900c08a0a477e4019e0977f10873a2f9`.
- Starting worktree: clean.

M4 normalizes `MEMORY.md` into durable cross-phase knowledge and advances
CURRENT to M5. It changes no science, runtime evidence, production behavior,
analysis source, validation profile, or collector.

## Before and after metrics

- Starting MEMORY: 14,387 bytes and 235 physical lines.
- Normalized MEMORY: 5,041 bytes and 93 physical lines.

## Retained durable categories

- post-M3 authority and evidence precedence;
- provenance-role semantics and dynamic identity rules;
- scientific ownership and production ordering;
- Method-A/Method-B and detached-production boundaries;
- diagnostic provenance and farm/runtime validation rules; and
- canonical-record ownership pointers.

## Removed chronology categories and canonical ownership audit

| Removed MEMORY category | Canonical retained owner(s) | Evidence strength preserved | Key provenance preserved | Verification result |
| --- | --- | --- | --- | --- |
| Phase-C five-setting closure | `evidence/phase-c-five-setting-closure.md` | Direct recovered runtime closure | Gate scope, settings, Method-B non-promotion | PASS |
| E.3.Fix.2 Left-low closure | `evidence/e3-fix2-left-low-runtime.md` | Direct farm gate | Exact setting, presentation-only boundary | PASS |
| F.1.Fix.5 closure and identities | `evidence/f1-fix5-v4-bundle-inspection-2026-09-14.md`; `phases/phase-f1-method-a-acceptance-farm-gate.md` | Owner-accepted runtime closure | Source/runtime/collector roles and detached populations | PASS |
| F.2 representation and accepted basis | `evidence/f2-fix1-runtime-closure.md`; `phases/phase-f2-method-a-acceptance-representation.md` | Accepted runtime closure | `hgcer3`, support, diagnostic-only alternatives | PASS |
| F.3 response-map closure | `evidence/f3-runtime-closure.md`; `phases/phase-f3-method-a-support-aware-acceptance-map.md` | Accepted runtime closure | Parent support continuity and detached map boundary | PASS |
| F.4 parent-preserving correction | `evidence/f4-runtime-closure.md`; `phases/phase-f4-method-a-parent-preserving-correction.md` | Accepted runtime closure | Parent-only normalization and no production application | PASS |
| F.5 `(t,phi)` propagation | `evidence/f5-runtime-closure.md`; `phases/phase-f5-method-a-tphi-propagation.md` | Accepted runtime closure | Cell population, closure, detached aggregate boundary | PASS |
| F.5.Fix.1 compatibility repair | `phases/phase-f5-method-a-tphi-propagation.md` | Source-reviewed repair history | `zip(strict=True)` limitation and unchanged gates | PASS |
| F.5.2 presentation rerender | `evidence/f5-2-runtime-closure.md`; `phases/phase-f5-method-a-tphi-propagation.md` | Accepted runtime closure | Presentation-only status and frozen propagation payload | PASS |
| F.6.1 and Validation.1 | `evidence/f6-1-runtime-closure.md`; `phases/phase-f6-method-a-production-promotion.md` | Accepted runtime closure | Detached aggregate-only validation and population facts | PASS |
| F.6.2 scientific closure | `evidence/f6-2-scientific-runtime-closure.md`; `phases/phase-f6-method-a-production-promotion.md` | Accepted scientific runtime closure | Detached scientific scope and no promotion | PASS |
| F.6.2.Fix.5 presentation closure | `evidence/f6-2-fix5-presentation-runtime-closure.md`; `phases/phase-f6-method-a-production-promotion.md` | Accepted presentation runtime closure | Scientific versus renderer/profile provenance | PASS |
| E.8/E.8.1 presentation state | `CURRENT.md`; `roadmap/CURRENT.md` | Active state plus approved roadmap | Reviewed source/profile roles and presentation-only boundary | PASS |
| F.6.3/F.6.4 dependency ordering | `CURRENT.md`; `roadmap/CURRENT.md` | Active state plus approved roadmap | Blocker chain and non-automatic promotion | PASS |
| F.2.Fix.1 authority/provenance repair | `phases/phase-f2-method-a-acceptance-representation.md` | Source/contract record | Reconstruction, deterministic failure, and output-safety rules | PASS |
| Generic validation collector/profile description | `phases/phase-f2-method-a-acceptance-representation.md`; `evidence/f2-fix1-runtime-closure.md` | Source and accepted-bundle records | Profile-declared artifacts and no analyzer execution | PASS |

## Identifier and non-hash preservation audit

The starting committed MEMORY contained 22 distinct full 40-character tokens
and 4 distinct full 64-character tokens. Normalized MEMORY retains zero of
either form, so 26 identifiers were removed from it. A deterministic
starting-HEAD search excluded MEMORY and the manifest and found at least one
pre-existing canonical Markdown owner for every removed identifier: evidence
for accepted artifacts/runtime identities, phase records for implementation
history, and CURRENT or roadmap for active-frontier source roles. The new M4
record was not used as ownership proof.

The non-hash audit compared role, accepted/rejected status, detached versus
production boundary, setting/parent/cell population, accepted basis,
normalization/closure rule, diagnostic-only restriction, presentation-only
restriction, provenance role, and dependency ordering against the owners in
the table. No substantive unique fact was removed without equal or stronger
canonical retention.

## Validation

Using the discovered local Python 3.12 interpreter represented as `<PYTHON>`
with `-B`, all pre-change M3 baseline commands returned 0:

- `<PYTHON> -B tools/check_memory_health.py --root .` — `MEMORY HEALTH: PASS`.
- `<PYTHON> -B tools/update_memory_manifest.py --root . --check` — `MANIFEST: PASS`.
- `<PYTHON> -B tools/memory_bootstrap.py --root . --json` — schema-3 summary
  with passing health.
- `<PYTHON> -B -m unittest testing.test_memory_health -v` — 36 tests passed.
- Each corresponding `--self-test` command — `SELF-TEST: PASS`.

After normalization, all final commands returned 0:

- `<PYTHON> -B -m unittest testing.test_memory_health -v` — 36 tests passed.
- Each of the three `--self-test` commands — `SELF-TEST: PASS`.
- `<PYTHON> -B tools/update_memory_manifest.py --root . --write` and the
  subsequent `--check` — manifest written and `MANIFEST: PASS`.
- `<PYTHON> -B tools/check_memory_health.py --root .` — `MEMORY HEALTH: PASS`.
- `<PYTHON> -B tools/memory_bootstrap.py --root . --json` — schema-3 summary
  and passing health.
- `git diff --check` — passed.

## Preserved scientific frontier

F.6.2 and F.6.2.Fix.5 remain `CLOSED / RUNTIME VALIDATED`; E.8.1 remains
`DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; F.6.3 remains blocked pending
E.8.1; F.6.4 remains blocked pending F.6.3 evidence. The accepted yield
remains the baseline, Method B remains diagnostic-only, and Method A remains
detached and non-production.

## NEXT

NEXT — Memory M5: Linux/JLab operational-role separation.
