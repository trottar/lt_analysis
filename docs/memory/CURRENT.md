---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

Complete the staged KaonLT repository-memory M0-M9 refinement so a fresh
ChatGPT/Codex session can recover authoritative state and the exact task
without a giant continuation prompt, while preserving the accepted scientific
and runtime frontier.

## Current Work Item

Memory M4 — durable MEMORY normalization. M3 activated schema 3, made
`CURRENT.md` the sole ordinary active-state authority, and removed ordinary
active-state mirroring; see the [M3 phase record](phases/memory-m3-single-active-state.md).

## Verified State

- M0, M1, and M2 are `SOURCE REVIEWED` repository-memory work.
- M3 is `SOURCE REVIEWED`, contingent on the local acceptance criteria in its
  [phase record](phases/memory-m3-single-active-state.md) and the subsequent
  ChatGPT actual-diff audit/user push workflow; it is not runtime validation.
- F.6.2 is `CLOSED / RUNTIME VALIDATED`; see the direct
  [scientific closure evidence](evidence/f6-2-scientific-runtime-closure.md).
- F.6.2.Fix.5 is `CLOSED / RUNTIME VALIDATED`; see the direct
  [presentation closure evidence](evidence/f6-2-fix5-presentation-runtime-closure.md).
- E.8.1 is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. The accepted yield
  remains the baseline; Method B remains diagnostic/cross-check only; Method A
  remains detached and non-production.

## Source / Evidence Identity

- E.8.1 reviewed procedure-PDF source:
  `0bf1281ebf44868b68ff808090436653eea6ed60`.
- E.8.1 canonical-five-setting bundle/profile source:
  `c25f9d8248f2acfb8b0443a482fc76781eb405ca`.
- F.6.2 reviewed scientific source:
  `0b37af2a2927b08bdeaf897c545f290b55329cea`; its accepted closure owners are
  the [scientific evidence](evidence/f6-2-scientific-runtime-closure.md) and
  [Fix.5 presentation evidence](evidence/f6-2-fix5-presentation-runtime-closure.md).
- The [schema-v3 contract](decisions/memory-schema-v3-contract.md) defines
  representation ownership; the [M3 phase record](phases/memory-m3-single-active-state.md)
  records this migration. These are role-labeled historical identities, not
  assertions about the dynamically queried live HEAD.

## Blockers

No known memory-migration blocker prevents M4 after M3 source acceptance.
Lifecycle-hook dispatch remains `BLOCKED / DEFERRED` and is not a dependency.
Scientifically, E.8.1 remains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
F.6.3 remains `BLOCKED` pending E.8.1, and F.6.4 remains `BLOCKED` pending
F.6.3 evidence. The preserved scientific continuation is the five-setting
`Q4p4W2p74` E.8.1 farm gate, with `Left / highe` inspected first, but it is
not the active repository-development action during this migration.

## Next Action

NEXT — Memory M4: normalize MEMORY.md into durable cross-phase knowledge only,
moving detailed chronology/current-gate material to verified canonical
phase/evidence/decision records without losing unique scientific provenance.

## Success Criteria

M4 must preserve durable scientific ownership and evidence-precedence rules,
trace detailed chronology to canonical deep records before removing it from
MEMORY, preserve every unique provenance hash and identity, avoid changing
scientific status or production behavior, and leave `CURRENT.md` as the sole
ordinary active-state authority.

## Do Not Reopen Without New Evidence

Do not reopen accepted F.1 through F.6.2 science without concrete regression
evidence, alter the accepted baseline yield, promote Method A, or turn Method B
into production. E.8 remains presentation-only and consumes frozen accepted
F.6.2 science. Later memory phases must not rewrite canonical runtime evidence
merely for style.

## Relevant References

- [M0-M9 refinement plan](decisions/memory-system-refinement-m0-m9-plan.md)
- [Schema-v3 contract](decisions/memory-schema-v3-contract.md)
- [M3 single-active-state record](phases/memory-m3-single-active-state.md)
- [F.6.2 scientific runtime closure](evidence/f6-2-scientific-runtime-closure.md)
- [F.6.2.Fix.5 presentation runtime closure](evidence/f6-2-fix5-presentation-runtime-closure.md)
- [F.6 production-promotion phase record](phases/phase-f6-method-a-production-promotion.md)
