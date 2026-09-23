---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

Complete the E.8.1 presentation-only farm-validation gate for frozen F.6.2
science across the canonical five Q4p4W2p74 settings without changing the
accepted baseline yield or detached Method-A/Method-B boundaries.

## Current Work Item

E.8.1 is ACTIVE. The fresh targeted Q4p4W2p74 / Left / lowe Fix.4 bundle at
fdd368f2084c9b9508e7ce679b8f7391f5b556f1 passed provenance, frozen-artifact,
checker, page-inventory, context, handoff, and map-page gates. It uses the
reviewed Fix.3 procedure source 350c34c55b2de33ad01011559dc6d8ed84d9c8a7 and
the accepted frozen F.6.2 JSON
5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1.

Independent PDFium and Poppler inspection still blocks only persisted-overlay
pages 38/41/44: their header and L/B/A legend are absent and the first plot row
is clipped. E.8.1.Fix.5 is SOURCE REVIEWED and user-pushed at
53fd262b730af8f1254e411a38231aebeb6a1da3. ChatGPT independently confirmed
that exact pushed source retains the reviewed narrow renderer/test repair;
Codex-reported deterministic checks were NOT RUN by ChatGPT. No Fix.5
ROOT/PyROOT, farm, or runtime result is claimed. E.8.1.Fix.6 is SOURCE
REVIEWED: independent ChatGPT actual-diff review found the profile/test
provenance re-pin PASS. Codex checks were NOT RUN by ChatGPT. No Fix.6
ROOT/PyROOT/farm/runtime result is claimed; the profile requires the distinct Fix.5 source
53fd262b730af8f1254e411a38231aebeb6a1da3, and no Fix.6 commit exists yet.
See the [Fix.4 overlay blocker evidence](evidence/e8-1-fix4-left-lowe-overlay-blocker.md).

## Verified State

- M0 through M9 are SOURCE REVIEWED repository-memory work; M9 is not
  farm/runtime validation. See [M9 audit evidence](evidence/memory-m9-fresh-session-audit.md).
- F.6.2 and F.6.2.Fix.5 are CLOSED / RUNTIME VALIDATED; see the
  [scientific](evidence/f6-2-scientific-runtime-closure.md) and
  [presentation](evidence/f6-2-fix5-presentation-runtime-closure.md) closure evidence.
- E.8.1 remains ACTIVE, presentation-only, and preserves the baseline yield;
  Method B is diagnostic/cross-check only and Method A is detached/non-production.
- E.8.1.Fix.1 is CLOSED / RUNTIME VALIDATED only for its frozen persisted
  parent-setting reader defect. E.8.1.Fix.2 is CLOSED / RUNTIME VALIDATED only
  for its profile provenance re-pin. Neither closure accepts E.8.1.
- E.8.1.Fix.3 is SOURCE REVIEWED and user-pushed at
  350c34c55b2de33ad01011559dc6d8ed84d9c8a7. The fresh bundle confirms its
  context/handoff repair; its earlier persisted-overlay subrepair is superseded
  by Fix.5.
- E.8.1.Fix.4 is CLOSED / RUNTIME VALIDATED only for its narrow profile
  provenance re-pin at fdd368f2084c9b9508e7ce679b8f7391f5b556f1. This is not
  E.8.1 acceptance.
- E.8.1.Fix.5 is SOURCE REVIEWED and user-pushed at
  53fd262b730af8f1254e411a38231aebeb6a1da3. ChatGPT independently confirmed
  that exact pushed source retains the reviewed renderer/test repair;
  Codex-reported deterministic checks were NOT RUN by ChatGPT. No ROOT/PyROOT,
  farm, or runtime validation is claimed.
- E.8.1.Fix.6 is SOURCE REVIEWED: ChatGPT independently inspected its complete
  actual diff and found the profile/test provenance re-pin PASS. Codex-reported
  deterministic checks were NOT RUN by ChatGPT. No ROOT/PyROOT, farm, or
  runtime validation is claimed; the user has not committed/pushed Fix.6.
- E.8.1.Debug.1 at 4c4a18a098d184fe0b13dd1fc57010ac7f5b30c0 was reviewed
  with a Diamond-cut blocker and is superseded; E.8.1.Debug.1.Fix.1 remains
  SOURCE REVIEWED. Its earlier incomplete output does not prove launcher
  provenance or the intentional full-high-epsilon skip.
- The parameterized validation-bundle wrapper is SOURCE REVIEWED at
  31abff9b8ab02f715cf7d6285a8dd2f53855a5e4; that is source review only.
- F.6.3 remains BLOCKED pending E.8.1.

## Source / Evidence Identity

- E.8.1 procedure-PDF source before Fix.1:
  0bf1281ebf44868b68ff808090436653eea6ed60.
- E.8.1.Fix.1 reader/test source: 8985d9a212799c021c4ad1a759a689ea3826e0ea.
- E.8.1.Fix.2 pushed profile-repin source: 9f7094b19f3cdb0887c95f1f972d98471f2ecfdc.
- E.8.1.Fix.3 pushed procedure-PDF source:
  350c34c55b2de33ad01011559dc6d8ed84d9c8a7.
- E.8.1.Fix.4 pushed profile/bundle source:
  fdd368f2084c9b9508e7ce679b8f7391f5b556f1; it requires the distinct Fix.3
  analysis/procedure source above.
- E.8.1.Fix.5 pushed analysis/procedure source:
  53fd262b730af8f1254e411a38231aebeb6a1da3. It remains SOURCE REVIEWED only;
  the future Fix.6 profile/bundle commit is unknown and must remain distinct
  from this required analysis source.
- The SOURCE REVIEWED E.8.1.Fix.6 profile requires that distinct Fix.5 source
  53fd262b730af8f1254e411a38231aebeb6a1da3. No Fix.6 commit exists yet, and
  its future profile/bundle identity remains unknown.
- E.8.1 canonical-five-setting bundle/profile source:
  c25f9d8248f2acfb8b0443a482fc76781eb405ca.
- F.6.2 reviewed scientific source: 0b37af2a2927b08bdeaf897c545f290b55329cea.
- The fresh Fix.4 bundle PDF SHA-256 is
  9a00fc6179891e936523613293dc40ebe5790f7ac49acf37aa942aea0333b98c; its
  page-manifest SHA-256 is
  d4bc89f49ce374f63055afd974d734bf782cea6e55cefcc26caf3c2f83028dc1.

## Blockers

E.8.1 remains ACTIVE because fresh farm PDF evidence blocks the persisted-overlay
presentation on pages 38/41/44 only. The accepted F.6.2 artifact, reader,
profile provenance, context/handoff pages, and map pages are not blockers and
must not be redesigned. F.6.3 remains BLOCKED pending E.8.1; F.6.4 remains
BLOCKED pending F.6.3 evidence. Lifecycle-hook dispatch is non-required and
BLOCKED / DEFERRED.

## Next Action

NEXT —

1. User-controlled commit/push of the reviewed Fix.6 profile/test/memory set.
2. ChatGPT pushed-state review.
3. Use the resulting pushed Fix.6 commit as the wrapper `--bundle-commit`.
4. Keep the profile `required_analysis_commit` equal to the distinct Fix.5
   source 53fd262b730af8f1254e411a38231aebeb6a1da3.
5. Run fresh targeted Q4p4W2p74 / Left / lowe.
6. Collect a fresh validation bundle and procedure PDF.
7. Independently inspect pages 37--47.
8. Require pages 38/41/44 to show the E.8 header, L/B/A legend, first
   canonical child title, and unclipped first row.
9. Only after Left/lowe visual PASS return to detailed Left / highe.
10. Broaden only after that detailed gate passes.

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
- [E.8.1.Fix.5 overlay geometry record](phases/e8-1-fix5-overlay-pdf-geometry.md)
- [E.8.1.Fix.6 profile re-pin record](phases/e8-1-fix6-bundle-profile-repin.md)
- [E.8.1.Fix.4 Left/lowe overlay blocker](evidence/e8-1-fix4-left-lowe-overlay-blocker.md)
