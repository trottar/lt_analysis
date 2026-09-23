---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

Complete the E.8.1 presentation-only farm-validation gate for the frozen
F.6.2 science across the canonical five `Q4p4W2p74` settings, while preserving
the accepted baseline yield and detached Method-A/Method-B boundaries.

## Current Work Item

E.8.1 is `ACTIVE`: the fresh targeted `Q4p4W2p74 / Left / lowe` debug bundle
reached ordinary procedure rendering but exposed a fail-closed E.8 reader/test
blocker. The frozen accepted F.6.2 JSON SHA-256 remained
`5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`, but the
reader rejected its persisted parent-setting schema with
`frozen_f6_2_authority_rejected:parent_setting_identity_invalid`; its
`complete=false` page manifest ended at `full_background.e8.unavailable`.
E.8.1.Fix.1 is `SOURCE REVIEWED` and user-pushed at
`8985d9a212799c021c4ad1a759a689ea3826e0ea`: ChatGPT inspected its actual
complete diff and found the narrow reader/test repair PASS. E.8.1.Fix.2 is
`SOURCE REVIEWED`: ChatGPT inspected its actual diff and found the
profile/test provenance repair PASS. Codex checks were `NOT RUN by ChatGPT`;
no ROOT/PyROOT, farm, or runtime validation is claimed for either repair. No
Fix.2 commit exists. This is neither F.6.2 science reopening nor E.8.1 runtime
acceptance. See the
[E.8.1.Fix.1 reader record](phases/e8-1-fix1-frozen-f6-2-setting-reader.md)
and [E.8.1.Fix.2 profile record](phases/e8-1-fix2-bundle-profile-repin.md).

## Verified State

- M0 through M9 are `SOURCE REVIEWED` repository-memory work. M9 completed the
  repository-only fresh-session integration audit; it is not farm/runtime
  validation. See the [M9 audit evidence](evidence/memory-m9-fresh-session-audit.md)
  and [M9 phase record](phases/memory-m9-fresh-session-integration.md).
- F.6.2 is `CLOSED / RUNTIME VALIDATED`; see the direct
  [scientific closure evidence](evidence/f6-2-scientific-runtime-closure.md).
- F.6.2.Fix.5 is `CLOSED / RUNTIME VALIDATED`; see the direct
  [presentation closure evidence](evidence/f6-2-fix5-presentation-runtime-closure.md).
- E.8.1 is `ACTIVE`: farm evidence exposed a concrete presentation-reader/test
  blocker. The accepted yield remains the baseline; Method B remains
  diagnostic/cross-check only; Method A remains detached and non-production.
- E.8.1.Fix.1 is `SOURCE REVIEWED`: ChatGPT inspected the actual complete diff
  and found the reader/test repair PASS; user pushed it at
  `8985d9a212799c021c4ad1a759a689ea3826e0ea`. Codex tests were `NOT RUN by
  ChatGPT`; no ROOT/PyROOT, farm, or runtime validation is claimed.
- E.8.1.Fix.2 is `SOURCE REVIEWED`: ChatGPT inspected its complete actual diff
  and found the profile/test provenance repair PASS. Codex checks were `NOT RUN
  by ChatGPT`; no ROOT/PyROOT, farm, or runtime validation is claimed; no Fix.2
  commit exists until user commit/push.
- E.8.1.Debug.1 at `4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0` was reviewed
  with a Diamond-cut blocker: its pre-Step-2 selector changed the cut owner to
  a Left-low fallback. It is superseded by Fix.1 and has no farm/runtime
  evidence.
- E.8.1.Debug.1.Fix.1 is `SOURCE REVIEWED`; its fresh Left/lowe output reached
  procedure rendering, but the incomplete bundle does not establish full
  launcher provenance or the intentional full-high-epsilon skip.
- The parameterized farm validation-bundle wrapper is `SOURCE REVIEWED` at
  `31abff9b8ab02f715cf7d6285a8dd2f53855a5e4`. This is source review only, not
  farm, ROOT/PyROOT, or E.8.1 runtime validation; it changes neither E.8.1's
  source/profile identity nor any scientific/runtime status. See the
  [post-push state reconciliation](phases/farm-validation-bundle-wrapper-fix1-state-reconciliation.md).

## Source / Evidence Identity

- E.8.1 reviewed procedure-PDF source before the reader repair:
  `0bf1281ebf44868b68ff808090436653eea6ed60`.
- E.8.1.Fix.1 user-pushed, `SOURCE REVIEWED` reader/test source:
  `8985d9a212799c021c4ad1a759a689ea3826e0ea`.
- E.8.1 canonical-five-setting bundle/profile source:
  `c25f9d8248f2acfb8b0443a482fc76781eb405ca`.
- E.8.1.Debug.1 source (reviewed-with-blocker original implementation,
  superseded by Fix.1):
  `4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0`.
- Parameterized farm validation-bundle wrapper source (`SOURCE REVIEWED` only):
  `31abff9b8ab02f715cf7d6285a8dd2f53855a5e4`.
- F.6.2 reviewed scientific source:
  `0b37af2a2927b08bdeaf897c545f290b55329cea`; its accepted closure owners are
  the [scientific evidence](evidence/f6-2-scientific-runtime-closure.md) and
  [Fix.5 presentation evidence](evidence/f6-2-fix5-presentation-runtime-closure.md).
- Fresh E.8.1.Debug.1.Fix.1 Left/lowe bundle:
  `KaonLT_E8_1_Debug1_Fix1_Q4p4W2p74_Left_lowe_20260923-012755.zip`, SHA-256
  `490036e8d6c53d6e025da99a5c6c4ddcd529b6573007534efb7aaceac1acd002`,
  `complete=false`; it preserved F.6.2 identity but exposed the reader/test
  blocker in the [E.8.1.Fix.1 reader record](phases/e8-1-fix1-frozen-f6-2-setting-reader.md).
- The [schema-v3 contract](decisions/memory-schema-v3-contract.md) defines
  representation ownership; the [M8 phase record](phases/memory-m8-manifest-bootstrap-semantics.md)
  records integrity-only manifest and dynamic bootstrap semantics. These are role-labeled
  historical identities, not assertions about the dynamically queried live HEAD.

## Blockers

E.8.1 is `ACTIVE`: the reader/test repair and profile re-pin are source-reviewed,
but Fix.2 remains unpushed and no E.8.1 farm-evaluated commit exists. F.6.3
remains `BLOCKED` pending E.8.1, and F.6.4 remains `BLOCKED` pending F.6.3
evidence. Lifecycle-hook dispatch remains `BLOCKED` / `DEFERRED`, is
non-required, and is not a dependency.

## Next Action

NEXT — E.8.1.Fix.2: user commit/push the `SOURCE REVIEWED` profile-only
provenance re-pin, then ChatGPT inspect the pushed commit. That resulting Fix.2
commit, not the Fix.1 analysis source commit, becomes the wrapper
`--bundle-commit` for a fresh `Q4p4W2p74 / Left / lowe` bundle-only farm gate.
Only after that narrow gate passes return to detailed `Left / highe` review and
then broaden.

## Success Criteria

E.8.1 requires direct fresh Jefferson Lab farm evidence for the canonical
five-setting package, with the reviewed source/profile provenance, frozen F.6.2
JSON identity, checker gates, requested structured artifacts, and readable
ordinary procedure-PDF/pages independently inspected. Following the required
profile re-pin, `Left / lowe` is the first rerun gate; `Left / highe` is the
next detailed gate. A completed ZIP or local source review alone is not farm
acceptance; do not broaden after a failed or incomplete gate.

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
- [E.8.1.Debug.1.Fix.1 local repair record](phases/e8-1-debug-1-fix1-preserve-diamond-cut.md)
- [E.8.1.Fix.1 frozen-F.6.2 reader record](phases/e8-1-fix1-frozen-f6-2-setting-reader.md)
- [E.8.1.Fix.2 bundle-profile re-pin record](phases/e8-1-fix2-bundle-profile-repin.md)
