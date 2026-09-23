---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

E.8 is ACTIVE: build a presentation-only, complete visual audit of the kaon
missing-mass and signal-region yield chain while preserving the accepted
baseline production branch, frozen F.6.2 science, and detached
Method-A/Method-B boundaries.

## Current Work Item

The fresh `Q4p4W2p74 / Left / lowe` Fix.6 package
`KaonLT_E8_1_Fix6_Q4p4W2p74_Left_lowe_20260923-090850.zip` is direct, reviewed
Jefferson Lab runtime evidence. Its bundle/profile commit is
0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b, its distinct required
analysis/procedure source is 53fd262b730af8f1254e411a38231aebeb6a1da3, its
captured worktree is clean, and its complete manifest has no errors. The
frozen F.6.2 JSON and both accepted fingerprints are unchanged.

This evidence closes E.8.1.Fix.5 (`CLOSED / RUNTIME VALIDATED`) for its
persisted-overlay geometry repair and E.8.1.Fix.6 (`CLOSED / RUNTIME VALIDATED`)
for its narrow profile/provenance gate. Real-PyROOT regression and pages 37--47
passed review, including the context, overlays 38/41/44, maps, and handoff.
It does not close E.8.1 as a canonical-five setting program: the remaining
canonical-five expansion is `DEFERRED` by user decision, not failed.
See the [Fix.6 Left/lowe runtime closure](evidence/e8-1-fix6-left-lowe-runtime-closure.md).

## Verified State

- F.6.2 and F.6.2.Fix.5 are `CLOSED / RUNTIME VALIDATED`; their accepted
  scientific JSON SHA-256, artifact fingerprint, and validation fingerprint
  remain frozen. See the [scientific closure](evidence/f6-2-scientific-runtime-closure.md)
  and [presentation closure](evidence/f6-2-fix5-presentation-runtime-closure.md).
- E.8 remains `ACTIVE` and presentation-only. The baseline production branch
  remains authoritative; Method A is detached/non-production and Method B is
  diagnostic/cross-check only.
- E.8.1.Fix.1 and E.8.1.Fix.2 retain their narrow `CLOSED / RUNTIME VALIDATED`
  closures. E.8.1.Fix.3 remains `SOURCE REVIEWED`; its context/handoff repair
  has farm evidence and its earlier overlay subrepair was superseded by Fix.5.
- E.8.1.Fix.4 remains `CLOSED / RUNTIME VALIDATED` only for its prior narrow
  profile-provenance gate. Fix.5 and Fix.6 have the distinct narrow runtime
  closures above; neither grants canonical-five E.8.1 closure.
- E.8.1.Debug.1 remains superseded after its Diamond-cut blocker;
  E.8.1.Debug.1.Fix.1 remains `SOURCE REVIEWED`. Its earlier incomplete output
  does not prove launcher provenance or the intentional full-high-epsilon skip.
- The parameterized validation-bundle wrapper remains `SOURCE REVIEWED` at
  31abff9b8ab02f715cf7d6285a8dd2f53855a5e4; that is source review only.
- F.6.3 is `BLOCKED` pending source-reviewed E.8.2 and E.8.3, not final E.8
  closure. E.8.4 is `BLOCKED` pending F.6.3; final E.8 is `BLOCKED` pending
  E.8.4 and its later runtime/visual gate; F.6.4 is `BLOCKED` pending that full
  production-impact evidence. Lifecycle-hook dispatch remains BLOCKED / DEFERRED.

## Source / Evidence Identity

- Fix.5 analysis/procedure source:
  53fd262b730af8f1254e411a38231aebeb6a1da3.
- Fix.6 pushed profile/bundle source:
  0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b. Its profile requires the distinct
  Fix.5 source above; source ancestry passed with no unexpected committed files.
- Fresh Fix.6 PDF SHA-256:
  d809a80fb1a46602dbb6cb930ed1d6c8932c2583ce8c5158ac82a78e6e6301e5.
  The 47-page manifest SHA-256 is
  d4bc89f49ce374f63055afd974d734bf782cea6e55cefcc26caf3c2f83028dc1 with
  `renderer_failures=[]`.
- Frozen F.6.2 JSON SHA-256:
  5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1;
  artifact fingerprint:
  ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0;
  validation fingerprint:
  7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b.

## Blockers

E.8.2 needs a true authoritative pre-proton snapshot or producer path; it must
not infer one from snapshots created after accepted proton cleaning. E.8.3 is
blocked on source-reviewed E.8.2; F.6.3 is blocked on source-reviewed E.8.2 and
E.8.3; E.8.4 and final E.8 remain downstream. The deferred E.8.1 settings are
not a blocker or a failure for this roadmap.

## Next Action

NEXT — E.8.2 baseline full-analysis stage audit: define the authoritative
before/component/after missing-mass and stage-yield presentation path for
prompt/random, dummy, slow proton, baseline pion, Fit 1, Fit 2, and final
canonical `(t,phi)` MM/yields. It must first obtain or consume a true
authoritative pre-proton snapshot; current yield-level filling already applies
accepted proton cleaning and must not be relabeled as pre-proton.

## Success Criteria

E.8.2 succeeds only by making the existing authoritative baseline stage chain
readable without recomputation. Subsequent E.8.3, F.6.3, and E.8.4 work must
retain source review, direct runtime evidence, and production promotion as
separate milestones.

E.8 presentation consumes authoritative upstream runtime objects, persisted
snapshots, accepted detached Method-A artifacts, or later F.6.3 branch outputs.
It must never recompute a fit, factor, correction, normalization, or yield for
plotting. F.6.3 alone constructs the parallel full-analysis Method-A branch by
changing only `w0_j` to `w0_j * C_j` in the pion-subtraction template; all
other production behavior stays frozen. The Method-A minus baseline shift is a
correction effect, not automatically an uncertainty or production promotion.

## Do Not Reopen Without New Evidence

Do not reopen accepted F.1 through F.6.2 science, alter the frozen F.6.2 JSON
or fingerprints, change the baseline production branch, promote Method A, or
make Method B numerical. The Fix.5/Fix.6 Left/lowe closure must not be expanded
into canonical-five acceptance without direct new runtime evidence.

## Relevant References

- [Fix.6 Left/lowe runtime closure](evidence/e8-1-fix6-left-lowe-runtime-closure.md)
- [F.6.2 scientific closure](evidence/f6-2-scientific-runtime-closure.md)
- [F.6.2.Fix.5 presentation closure](evidence/f6-2-fix5-presentation-runtime-closure.md)
- [E.8 full-analysis procedure roadmap](decisions/e8-full-analysis-procedure-roadmap.md)
- [Phase-F Method-A roadmap](phases/phase-f6-method-a-production-promotion.md)
- [Roadmap dependency status](roadmap/STATUS.md)
