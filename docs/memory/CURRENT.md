---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

Complete the E.8.1 presentation-only farm-validation gate for frozen F.6.2
science across the canonical five `Q4p4W2p74` settings without changing the
accepted baseline yield or detached Method-A/Method-B boundaries.

## Current Work Item

E.8.1 is `ACTIVE`. The fresh targeted `Q4p4W2p74 / Left / lowe` bundle at
`9f7094b19f3cdb0887c95f1f972d98471f2ecfdc` completed its provenance and
artifact gate: the frozen F.6.2 JSON remained
`5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`, the
ordinary procedure PDF had the expected available 47-page inventory, and the
farm module passed. It validates the narrow Fix.1 reader and Fix.2 profile
provenance defects, but independent PDF inspection found context/handoff text
clipping, context em-dash mojibake, and overlay header/grid clipping. Fix.3 is
`SOURCE REVIEWED` and user-pushed at
`350c34c55b2de33ad01011559dc6d8ed84d9c8a7`; no Fix.3 ROOT/PyROOT, farm, or
runtime result is claimed. Fix.4 is `SOURCE REVIEWED`: ChatGPT inspected its
complete actual profile/test diff and found it
PASS. Codex-reported checks were `NOT RUN by ChatGPT`; no Fix.4 ROOT/PyROOT,
farm, or runtime result is claimed, and no Fix.4 commit exists yet. See the
[Left/lowe blocker evidence](evidence/e8-1-fix2-left-lowe-layout-blocker.md).

## Verified State

- M0 through M9 are `SOURCE REVIEWED` repository-memory work; M9 is not
  farm/runtime validation. See [M9 audit evidence](evidence/memory-m9-fresh-session-audit.md).
- F.6.2 and F.6.2.Fix.5 are `CLOSED / RUNTIME VALIDATED`; see the
  [scientific](evidence/f6-2-scientific-runtime-closure.md) and
  [presentation](evidence/f6-2-fix5-presentation-runtime-closure.md) closure evidence.
- E.8.1 is `ACTIVE`, presentation-only, and still preserves the baseline yield;
  Method B is diagnostic/cross-check only and Method A is detached/non-production.
- E.8.1.Fix.1 is `CLOSED / RUNTIME VALIDATED` only for its frozen persisted
  parent-setting reader defect. Its source was independently `SOURCE REVIEWED`
  before the successful fresh Left/lowe artifact result; this does not accept E.8.1.
- E.8.1.Fix.2 is `CLOSED / RUNTIME VALIDATED` only for its profile provenance
  re-pin. Its source was independently `SOURCE REVIEWED` before the complete,
  fail-closed fresh bundle result; this does not accept E.8.1.
- E.8.1.Fix.3 is `SOURCE REVIEWED` and user-pushed at
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`; its pushed source review remains
  distinct from ROOT/PyROOT, farm, and runtime validation.
- E.8.1.Fix.4 is `SOURCE REVIEWED`: ChatGPT inspected its complete actual
  profile/test diff and found the provenance re-pin PASS. Codex-reported checks
  were `NOT RUN by ChatGPT`; no ROOT/PyROOT, farm, or runtime claim is made,
  and no Fix.4 commit exists yet.
- E.8.1.Debug.1 at `4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0` was reviewed
  with a Diamond-cut blocker and is superseded; E.8.1.Debug.1.Fix.1 remains
  `SOURCE REVIEWED`. Its earlier incomplete output does not prove launcher
  provenance or the intentional full-high-epsilon skip.
- The parameterized validation-bundle wrapper is `SOURCE REVIEWED` at
  `31abff9b8ab02f715cf7d6285a8dd2f53855a5e4`; that is source review only.

## Source / Evidence Identity

- E.8.1 procedure-PDF source before Fix.1:
  `0bf1281ebf44868b68ff808090436653eea6ed60`.
- E.8.1.Fix.1 reader/test source: `8985d9a212799c021c4ad1a759a689ea3826e0ea`.
- E.8.1.Fix.2 pushed profile-repin source:
  `9f7094b19f3cdb0887c95f1f972d98471f2ecfdc`; it requires Fix.1 source.
- E.8.1.Fix.3 pushed procedure-PDF source:
  `350c34c55b2de33ad01011559dc6d8ed84d9c8a7`.
- E.8.1.Fix.4 `SOURCE REVIEWED` profile re-pin requires that exact Fix.3 source;
  its future profile/bundle commit is not yet known and must remain distinct.
- E.8.1 canonical-five-setting bundle/profile source:
  `c25f9d8248f2acfb8b0443a482fc76781eb405ca`.
- F.6.2 reviewed scientific source: `0b37af2a2927b08bdeaf897c545f290b55329cea`.
- Fresh Fix.2 Left/lowe bundle:
  `KaonLT_E8_1_Fix2_Q4p4W2p74_Left_lowe_20260923-040427.zip`, SHA-256
  `52e0d283acdb57b56ed4c34c59f98dfe077f2a89e47fd0ced60fdcc91d5a5fd4`; its
  PDF SHA-256 is `0745172508e2131e9256425b05bc44b19ffa6bcaa649224328679a75b7e18007`
  and page-manifest SHA-256 is
  `d4bc89f49ce374f63055afd974d734bf782cea6e55cefcc26caf3c2f83028dc1`.
- The wrapper source is `31abff9b8ab02f715cf7d6285a8dd2f53855a5e4`.

## Blockers

E.8.1 remains `ACTIVE` because the reviewed Fix.4 profile must be
user-committed/pushed, independently reviewed in its pushed state, and followed
by a fresh farm PDF inspection. The current farm record is a presentation
blocker, not E.8.1 acceptance. F.6.3 remains `BLOCKED` pending E.8.1; F.6.4
remains `BLOCKED` pending F.6.3 evidence. Lifecycle-hook dispatch is
non-required and `BLOCKED` / `DEFERRED`.

## Next Action

NEXT — User-controlled commit/push of the `SOURCE REVIEWED` E.8.1.Fix.4 change
set; ChatGPT then inspects the actual pushed profile commit. That future Fix.4
commit, distinct from the profile's required Fix.3 source
`350c34c55b2de33ad01011559dc6d8ed84d9c8a7`, becomes the wrapper
`--bundle-commit`. Only after pushed-state review rerun `Q4p4W2p74 / Left /
lowe`, return its bundle/PDF for independent visual review, and inspect the
repaired clipping, mojibake, handoff, and overlay-header defects before
returning to detailed `Left / highe` or broadening.

## Success Criteria

E.8.1 requires direct fresh Jefferson Lab farm evidence for the canonical-five
package: reviewed source/profile provenance, frozen F.6.2 identity, checker
gates, required structured artifacts, and independently readable procedure-PDF
pages. Source review, a ZIP, and local tests alone do not provide farm acceptance.

## Do Not Reopen Without New Evidence

Do not reopen accepted F.1 through F.6.2 science, alter the accepted baseline
yield, promote Method A, or make Method B production. E.8 consumes frozen F.6.2
science and remains presentation-only; later memory must not rewrite runtime
evidence for style.

## Relevant References

- [M9 audit evidence](evidence/memory-m9-fresh-session-audit.md)
- [F.6.2 scientific closure](evidence/f6-2-scientific-runtime-closure.md)
- [F.6.2.Fix.5 presentation closure](evidence/f6-2-fix5-presentation-runtime-closure.md)
- [Farm validation-bundle procedure](decisions/farm-validation-bundle-procedure.md)
- [Wrapper reconciliation](phases/farm-validation-bundle-wrapper-fix1-state-reconciliation.md)
- [E.8.1.Debug.1 local debug record](phases/e8-1-debug-1-left-low-debug-mode.md)
- [E.8.1.Debug.1.Fix.1 local repair record](phases/e8-1-debug-1-fix1-preserve-diamond-cut.md)
- [E.8.1.Fix.1 reader record](phases/e8-1-fix1-frozen-f6-2-setting-reader.md)
- [E.8.1.Fix.2 profile record](phases/e8-1-fix2-bundle-profile-repin.md)
- [E.8.1.Fix.3 layout record](phases/e8-1-fix3-procedure-pdf-layout.md)
- [E.8.1.Fix.4 profile re-pin record](phases/e8-1-fix4-bundle-profile-repin.md)
- [E.8.1.Fix.2 Left/lowe blocker evidence](evidence/e8-1-fix2-left-lowe-layout-blocker.md)
