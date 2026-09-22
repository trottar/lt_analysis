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

Memory M6 — roadmap/source/history normalization. M5 established the
Linux/JLab operational ownership split among the repository control records.
See the [M5 phase record](phases/memory-m5-linux-jlab-operational-ownership.md).

## Verified State

- M0 through M4 are `SOURCE REVIEWED` repository-memory work.
- M5 is `SOURCE REVIEWED`, contingent on the local acceptance criteria in its
  [phase record](phases/memory-m5-linux-jlab-operational-ownership.md) and the
  subsequent ChatGPT actual-diff audit/user push workflow; it is not runtime
  validation.
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
  representation ownership; the [M5 phase record](phases/memory-m5-linux-jlab-operational-ownership.md)
  records the current operational-memory normalization. These are role-labeled
  historical identities, not assertions about the dynamically queried live HEAD.

## Blockers

No known memory-migration blocker prevents M6 after M5 source acceptance.
Lifecycle-hook dispatch remains `BLOCKED / DEFERRED` and is not a dependency.
Scientifically, E.8.1 remains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
F.6.3 remains `BLOCKED` pending E.8.1, and F.6.4 remains `BLOCKED` pending
F.6.3 evidence. The preserved scientific continuation is the five-setting
`Q4p4W2p74` E.8.1 farm gate, with `Left / highe` inspected first, but it is
not the active repository-development action during this migration.

## Next Action

NEXT — Memory M6: normalize roadmap/source/history ownership so roadmap
records dependencies/status only, source indexes contain stable provenance
rather than stored "live HEAD" claims, and chat/import/history surfaces remain
historical navigation rather than active authority.

## Success Criteria

M6 must retain roadmap dependency/status structure without a competing exact
active action; remove static "live HEAD" semantics; use explicit role labels
for fixed SHAs; keep chat/import/history as historical navigation; and change
no scientific or runtime status.

## Do Not Reopen Without New Evidence

Do not reopen accepted F.1 through F.6.2 science without concrete regression
evidence, alter the accepted baseline yield, promote Method A, or turn Method B
into production. E.8 remains presentation-only and consumes frozen accepted
F.6.2 science. Later memory phases must not rewrite canonical runtime evidence
merely for style.

## Relevant References

- [M0-M9 refinement plan](decisions/memory-system-refinement-m0-m9-plan.md)
- [Schema-v3 contract](decisions/memory-schema-v3-contract.md)
- [M5 operational-ownership record](phases/memory-m5-linux-jlab-operational-ownership.md)
- [F.6.2 scientific runtime closure](evidence/f6-2-scientific-runtime-closure.md)
- [F.6.2.Fix.5 presentation runtime closure](evidence/f6-2-fix5-presentation-runtime-closure.md)
- [F.6 production-promotion phase record](phases/phase-f6-method-a-production-promotion.md)
