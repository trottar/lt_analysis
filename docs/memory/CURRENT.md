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

Memory M9 — fresh-session integration audit. M8 finalized integrity-only
manifest and dynamic bootstrap semantics without changing scientific/runtime
state. See the [M8 phase record](phases/memory-m8-manifest-bootstrap-semantics.md).

## Verified State

- M0 through M7 are `SOURCE REVIEWED` repository-memory work.
- M8 is `SOURCE REVIEWED`, contingent on its local acceptance criteria and the
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
  representation ownership; the [M8 phase record](phases/memory-m8-manifest-bootstrap-semantics.md)
  records integrity-only manifest and dynamic bootstrap semantics. These are role-labeled
  historical identities, not assertions about the dynamically queried live HEAD.

## Blockers

No known memory blocker prevents M9 after M8 source acceptance.
Lifecycle-hook dispatch remains `BLOCKED` / `DEFERRED` and is not a dependency.
Scientifically, E.8.1 remains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`;
F.6.3 remains `BLOCKED` pending E.8.1, and F.6.4 remains `BLOCKED` pending
F.6.3 evidence. The preserved scientific continuation is the five-setting
`Q4p4W2p74` E.8.1 farm gate, with `Left / highe` inspected first, but it is
not the active repository-development action during this migration.

## Next Action

NEXT — Memory M9: perform a genuine fresh-session integration audit using
repository state only, proving recovery of repository identity requirements,
active objective/NEXT, scientific frontier/evidence boundaries, Method-A/
Method-B ownership, actor sequence, Linux/JLab context, and handoff semantics
without a giant chat-continuation prompt.

## Success Criteria

M9 must use a genuinely fresh ChatGPT/Codex-style session with no giant
continuation prompt; recover dynamic branch/HEAD/worktree requirements, the
active objective and exact NEXT from CURRENT, accepted scientific
frontier/evidence, Method B diagnostic-only and Method A detached boundaries,
production/diagnostic/presentation separation, source-review versus
farm-runtime distinction, Codex/ChatGPT/user/farm actor sequence, Linux/JLab
environment, and exceptional-handoff semantics. It must use task-directed
expansion rather than whole-tree loading and state NO CHAT CONTINUATION REQUIRED
only if that audit actually passes.

## Do Not Reopen Without New Evidence

Do not reopen accepted F.1 through F.6.2 science without concrete regression
evidence, alter the accepted baseline yield, promote Method A, or turn Method B
into production. E.8 remains presentation-only and consumes frozen accepted
F.6.2 science. Later memory phases must not rewrite canonical runtime evidence
merely for style.

## Relevant References

- [M0-M9 refinement plan](decisions/memory-system-refinement-m0-m9-plan.md)
- [Schema-v3 contract](decisions/memory-schema-v3-contract.md)
- [M8 manifest/bootstrap semantics record](phases/memory-m8-manifest-bootstrap-semantics.md)
- [F.6.2 scientific runtime closure](evidence/f6-2-scientific-runtime-closure.md)
- [F.6.2.Fix.5 presentation runtime closure](evidence/f6-2-fix5-presentation-runtime-closure.md)
- [F.6 production-promotion phase record](phases/phase-f6-method-a-production-promotion.md)
