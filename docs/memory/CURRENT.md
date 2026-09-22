---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

Complete the E.8.1 presentation-only farm-validation gate for the frozen
F.6.2 science across the canonical five `Q4p4W2p74` settings, while preserving
the accepted baseline yield and detached Method-A/Method-B boundaries.

## Current Work Item

E.8.1.Debug.1 is `ACTIVE`: a narrow local `-d` launcher implementation is the
immediate operational prerequisite before the existing E.8.1 farm gate. It
preserves paired canonical preflight, performs only the full `Left / lowe`
debug analysis, and stops before full high epsilon. It awaits independent
ChatGPT source review and the user-controlled Git handoff; it is not
`SOURCE REVIEWED` or farm/runtime validated. The scientific E.8.1 gate remains
unchanged behind this prerequisite; inspect `Q4p4W2p74 / Left / highe` first
before broader visual/runtime review. See the
[Debug.1 local phase record](phases/e8-1-debug-1-left-low-debug-mode.md).

## Verified State

- M0 through M9 are `SOURCE REVIEWED` repository-memory work. M9 completed the
  repository-only fresh-session integration audit; it is not farm/runtime
  validation. See the [M9 audit evidence](evidence/memory-m9-fresh-session-audit.md)
  and [M9 phase record](phases/memory-m9-fresh-session-integration.md).
- F.6.2 is `CLOSED / RUNTIME VALIDATED`; see the direct
  [scientific closure evidence](evidence/f6-2-scientific-runtime-closure.md).
- F.6.2.Fix.5 is `CLOSED / RUNTIME VALIDATED`; see the direct
  [presentation closure evidence](evidence/f6-2-fix5-presentation-runtime-closure.md).
- E.8.1 is `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`. The accepted yield
  remains the baseline; Method B remains diagnostic/cross-check only; Method A
  remains detached and non-production.
- E.8.1.Debug.1 is `ACTIVE` local-only launcher orchestration work. Its
  `-d` path has deterministic source-contract coverage only; no ROOT/PyROOT,
  full-analysis, farm, or runtime evidence exists.
- The parameterized farm validation-bundle wrapper is `SOURCE REVIEWED` at
  `31abff9b8ab02f715cf7d6285a8dd2f53855a5e4`. This is source review only, not
  farm, ROOT/PyROOT, or E.8.1 runtime validation; it changes neither E.8.1's
  source/profile identity nor any scientific/runtime status. See the
  [post-push state reconciliation](phases/farm-validation-bundle-wrapper-fix1-state-reconciliation.md).

## Source / Evidence Identity

- E.8.1 reviewed procedure-PDF source:
  `0bf1281ebf44868b68ff808090436653eea6ed60`.
- E.8.1 canonical-five-setting bundle/profile source:
  `c25f9d8248f2acfb8b0443a482fc76781eb405ca`.
- Parameterized farm validation-bundle wrapper source (`SOURCE REVIEWED` only):
  `31abff9b8ab02f715cf7d6285a8dd2f53855a5e4`.
- F.6.2 reviewed scientific source:
  `0b37af2a2927b08bdeaf897c545f290b55329cea`; its accepted closure owners are
  the [scientific evidence](evidence/f6-2-scientific-runtime-closure.md) and
  [Fix.5 presentation evidence](evidence/f6-2-fix5-presentation-runtime-closure.md).
- The [schema-v3 contract](decisions/memory-schema-v3-contract.md) defines
  representation ownership; the [M8 phase record](phases/memory-m8-manifest-bootstrap-semantics.md)
  records integrity-only manifest and dynamic bootstrap semantics. These are role-labeled
  historical identities, not assertions about the dynamically queried live HEAD.

## Blockers

E.8.1 remains `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`; no E.8.1
farm-evaluated commit exists. F.6.3 remains `BLOCKED` pending E.8.1, and F.6.4
remains `BLOCKED` pending F.6.3 evidence. Lifecycle-hook dispatch remains
`BLOCKED` / `DEFERRED`, is non-required, and is not a dependency.

## Next Action

NEXT — E.8.1.Debug.1: obtain independent ChatGPT source review of the narrow
local `-d` implementation, then complete the user-controlled Git handoff.
After this operational prerequisite is source-reviewed, pushed, and any
user-requested debug run has fresh evidence, resume the unchanged E.8.1
canonical-five-setting `Q4p4W2p74` farm-validation gate using the reviewed
wrapper; package and inspect `Left / highe` first, and broaden only after that
detailed visual/runtime review passes.

## Success Criteria

E.8.1 requires direct fresh Jefferson Lab farm evidence for the canonical
five-setting package, with the reviewed source/profile provenance, frozen F.6.2
JSON identity, checker gates, requested structured artifacts, and readable
ordinary procedure-PDF/pages independently inspected. `Left / highe` is the
first detailed gate. A completed ZIP or local source review alone is not farm
acceptance; do not broaden after a failed or incomplete first gate.

## Do Not Reopen Without New Evidence

Do not reopen accepted F.1 through F.6.2 science without concrete regression
evidence, alter the accepted baseline yield, promote Method A, or turn Method B
into production. E.8 remains presentation-only and consumes frozen accepted
F.6.2 science. Later memory phases must not rewrite canonical runtime evidence
merely for style.

## Relevant References

- [M9 fresh-session audit evidence](evidence/memory-m9-fresh-session-audit.md)
- [M9 fresh-session integration record](phases/memory-m9-fresh-session-integration.md)
- [Approved dependency/status roadmap](roadmap/STATUS.md)
- [F.6.2 scientific runtime closure](evidence/f6-2-scientific-runtime-closure.md)
- [F.6.2.Fix.5 presentation runtime closure](evidence/f6-2-fix5-presentation-runtime-closure.md)
- [F.6 production-promotion phase record](phases/phase-f6-method-a-production-promotion.md)
- [Farm validation-bundle procedure](decisions/farm-validation-bundle-procedure.md)
- [Parameterized wrapper phase record](phases/farm-validation-bundle-wrapper.md)
- [Wrapper post-push state reconciliation](phases/farm-validation-bundle-wrapper-fix1-state-reconciliation.md)
- [E.8.1.Debug.1 local debug-mode record](phases/e8-1-debug-1-left-low-debug-mode.md)
