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

F.4.Refresh.2.Validation.1 generic farm materialization bundle/profile is
`SOURCE REVIEWED`: independent ChatGPT actual-diff/source-provenance review of
`kaonlt_review(20260930-152818).diff` passed. Its five required global JSONs
pin pushed materializer source `141a3d04f9e5d07be21dba14e0e63212c3990bf1`.
Codex-reported local profile/collector/materializer tests passed (4/28/13,
0 skips); ChatGPT did not run those suites. F.4.Refresh.2/Fix.1 remain
`SOURCE REVIEWED`. No farm materialization, authority acceptance, or Method-A
promotion follows. See the [profile phase record](phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md).

F.4.Refresh.1 current-baseline Method-A authority comparison and its Fix.1
nested-delta repair are `CLOSED / RUNTIME VALIDATED` only for the detached
comparison gate. The supplied farm comparator
`KaonLT_F4_Refresh1_authority_comparison_Q4p4W2p74_20260930-134654.json`
(SHA-256 `c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5`)
reports F.2/F.3 scientific equality, F.4 scientific difference, and first
changed stage `F4`. See [comparator evidence](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md).

Earlier independent ChatGPT actual-diff/source-runtime-path review of
`kaonlt_review(20260930-090154).diff` PASSED the comparator and Fix.1 source.
The detached comparator requires canonical-five current F.1 artifacts, records
raw-byte SHA-256 provenance, rebuilds candidates through the existing public
F.2/F.3/F.4 builders, hashes serialized candidate F.2/F.3, and uses an
in-memory candidate-F.3 authority override only for diagnostic F.4 rebuilding.
It compares scientific payloads with declared provenance exclusions and reports
recursive F.4 parent/source/phi deltas. Accepted authorities, scientific and
production source, and Method-A promotion are unchanged. Codex-reported local
checks passed (Fix.1 comparator 8, F.4 8, F.6.3 21, memory 35 tests, with no
required-suite skips); they were `NOT RUN by ChatGPT`. The original unchanged
F.2 and F.3 suites had reported 14 and 10 passing tests. The later farm
comparison closes only the detached diagnostic question; it does not validate
a refreshed F.4 authority, F.6.3/E.8.4 runtime, or production. See the
[phase record](phases/f4-refresh1-current-baseline-authority-comparison.md)
and [direct Fix.4 evidence](evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md).

E.8.4.Fix.4 validation-bundle profile re-pin is `SOURCE REVIEWED`: independent
ChatGPT actual-diff/source-provenance review of
`kaonlt_review(20260929-222602).diff` passed. The reviewed candidate changes only the
generic E.8.1 bundle profile and focused test's required analysis source from
Fix.3 `29d7b7f9635db899939efeb3508e941e994e8928` to pushed Fix.4
`6e2adf7a37ac9e79cad99242686804cf51701644`. Canonical-five settings,
artifact inventory, allowlists, collector, wrapper, and scientific/runtime
source remain unchanged. Codex-reported local checks passed (5 profile, 28
collector, and 19 alignment tests with 11 PyROOT-dependent skips); they were
`NOT RUN by ChatGPT`. No ROOT/PyROOT, full `main.py`, procedure-PDF, farm,
runtime, production, or Method-A promotion acceptance follows. See the
[profile re-pin phase record](phases/e8-4-fix4-bundle-profile-repin.md).

E.8.4.Fix.4 persisted-alignment semantic-version repair is `CLOSED / RUNTIME VALIDATED`
for the narrow stale-cache rejection/recomputation gate. Fresh farm evidence
shows `rejected_stale_then_created` under current semantics with an
`alignment_semantics_version mismatch` stale-record reason. The prior
source review records an independent ChatGPT actual-diff/source-runtime-path PASS of
`kaonlt_review(20260929-215203).diff`. The reviewed local candidate
starts from committed `test` HEAD `8d62dcfdf2fd8d08298c075940b4ab1b28d7f079`.
It rejects pre-Fix.3 schema-v2 alignment caches for one-time recomputation,
retains current-semantics reuse, and requires current semantics for direct
fine-bin parent validity. The two reviewed source/test diffs remain byte-identical
to that bundle. Codex-reported deterministic checks passed, including 19
focused alignment tests with 11 PyROOT-dependent skips; these checks were
`NOT RUN by ChatGPT`. This narrow Fix.4 farm result does not establish
full `main.py`, procedure-PDF, F.6.3/E.8.4 runtime, production, or Method-A
promotion acceptance. See the [Fix.4 phase
record](phases/e8-4-fix4-alignment-cache-semantics.md) and [fresh Left/lowe
blocker](evidence/e8-4-fix3-left-lowe-stale-alignment-cache-runtime-blocker.md).
Independent pushed-state review passed for user-controlled pushed Fix.4 source
`6e2adf7a37ac9e79cad99242686804cf51701644`, whose parent is the reviewed
starting HEAD above. That pushed-state review is source/provenance evidence;
the narrow runtime closure comes from the fresh alignment record.

E.8.4.Fix.3 pion-alignment determinism repair is `SOURCE REVIEWED`: independent
ChatGPT actual-diff/source review and pushed-state review passed for the
user-controlled pushed source `29d7b7f9635db899939efeb3508e941e994e8928`.
It makes persisted pion-control compatibility depend on checksum plus axis
rather than an ephemeral generated ROOT histogram name, and makes only the
renormalized template minimum-integral boundary machine-scale roundoff safe.
It preserves schema-v2 cache compatibility, all other provenance checks, the
configured threshold, and existing scan/production physics. This remains no
ROOT/PyROOT, full-analysis, procedure-PDF, farm, or runtime validation. See the
[Fix.3 determinism blocker evidence](evidence/e8-4-left-lowe-pion-alignment-determinism-blocker.md)
and [Fix.3 phase record](phases/e8-4-fix3-pion-alignment-determinism.md).

E.8.4.Fix.3 validation-bundle profile re-pin is `SOURCE REVIEWED`: final
independent ChatGPT actual-diff/source-provenance review of
`kaonlt_review_e8_4_fix3_profile_repin.diff` passed. It changes only the
generic E.8.1 bundle profile/test required analysis source to the pushed Fix.3
source; the collector, wrapper, and scientific/runtime source remain frozen.
Codex-reported deterministic checks were `NOT RUN by ChatGPT`. This is not
ROOT/PyROOT, procedure-PDF, farm, runtime, production, or Method-A-promotion
evidence. See the [Fix.3 bundle-profile re-pin](phases/e8-4-fix3-bundle-profile-repin.md).

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
- E.8 remains `ACTIVE` and presentation-only. E.8.2 is `SOURCE REVIEWED`:
  independent ChatGPT inspection of the complete cumulative E.8.2/Fix.1/Fix.2/
  Fix.3 actual diff passed, and the pushed set at
  `91bb7809d27d84d7709a6600dbb0dc9ab514a458` matches that reviewed candidate.
  The source-reviewed baseline audit captures the true same-traversal pre-proton
  diagnostic, the actual post-proton pre-/post-prune production snapshots, the
  exact baseline-pion before/template/after objects, and the existing `Y_0`.
  Production physics, `w0`, proton factors, pruning, random/dummy normalization,
  canonical binning, public yields/errors, Method A, Method B, and frozen F.6.2
  science remain unchanged. Codex-reported local checks were `NOT RUN by
  ChatGPT`; source review is not ROOT/PyROOT, full-analysis, farm, or runtime
  acceptance.
- The accepted active background profile is `no_empirical_residual`; it forces
  both empirical residual-background scales to zero. Legacy empirical residual
  Fit 1/Fit 2 are dormant historical source machinery and are outside the
  current E.8 scientific/presentation chain.
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
- E.8.3 is `SOURCE REVIEWED` from
  `kaonlt_review(20260924-064150).diff`. F.6.3 is `SOURCE REVIEWED` from
  independent ChatGPT actual-diff/source-runtime-path review of
  `kaonlt_review(20260928-161516).diff`. Codex-reported deterministic checks
  were `NOT RUN by ChatGPT`; these source reviews do not claim ROOT/PyROOT,
  full-analysis, farm, or runtime acceptance. E.8.4 is `SOURCE REVIEWED`:
  ChatGPT passed complete `kaonlt_review(20260928-233040).diff` actual-diff/
  source-runtime-path review. The consumer retains Fix.1's identity/geometry/
  common-input repairs and Fix.2's malformed-sidecar fail-closed boundary; it
  cannot construct either branch. Codex checks were `NOT RUN by ChatGPT`; no
  ROOT/PyROOT, full-analysis, farm, or runtime claim. Final E.8 is
  `BLOCKED` pending its later runtime/visual gate; F.6.4 is `BLOCKED` pending
  production-impact evidence. Lifecycle-hook dispatch remains BLOCKED /
  DEFERRED.
- E.8.4.Fix.3 is `SOURCE REVIEWED`: independent actual-diff/source and
  pushed-state reviews passed for user-controlled pushed commit
  `29d7b7f9635db899939efeb3508e941e994e8928`. It is a baseline
  pion-alignment determinism/provenance repair, not a reopened F.1--F.6.2
  scientific result or a Method-A promotion. Its direct farm evidence remains
  a blocker record, not runtime closure. Its narrow validation-bundle profile
  re-pin is `SOURCE REVIEWED` after final independent actual-diff/
  source-provenance review passed; that review is not runtime or production
  evidence.
- The E.8.4 validation-bundle profile provenance re-pin is `SOURCE REVIEWED`:
  independent ChatGPT reviewed `kaonlt_review(20260929-100154).diff` and
  returned a passing source/provenance result. It requires pushed source
  `1aa1fd4184a6f8b20043e00ebb1ed3e9505a4935`, with the exact canonical-five
  artifact inventory and profile/test-only range allowlist intact. The review
  is source/provenance only; Codex-reported checks were `NOT RUN by ChatGPT`,
  and it claims no ROOT/PyROOT, procedure-PDF, farm, or runtime acceptance.
  The non-fatal `check_memory_health.py` CURRENT soft-size warning remains
  distinct from a source or runtime failure.

## Source / Evidence Identity

- Fix.5 analysis/procedure source:
  53fd262b730af8f1254e411a38231aebeb6a1da3.
- Fix.6 pushed profile/bundle source:
  0ec29d4e1bb345eb37e8cca35b8b7e5cbe1b4d5b. Its profile requires the distinct
  Fix.5 source above; source ancestry passed with no unexpected committed files.
- E.8.2 pushed source/test/memory set:
  91bb7809d27d84d7709a6600dbb0dc9ab514a458. Pushed-state review passed; this
  is source/push provenance only, not farm or runtime evidence.
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

The fresh `Q4p4W2p74 / Left / lowe` F.6.3/E.8.4 runtime gate remains `BLOCKED` by
`f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`.
The fresh farm evidence confirms Fix.4 stale-cache rejection/recomputation.
Record-by-record current/frozen F.1 comparison found exact F.3 training
(52,397) and application (55,380) projections, but all 55,380 F.4/F.6.3
baseline projections differ in `analysis_MM`, `baseline_pion_weight_w0`, and
`signed_baseline_event_contribution`. F.4.Refresh.1 has established F.4 as
the first scientifically changed stage. F.4.Refresh.2 is source reviewed only;
farm materialization and authority review remain before any authority decision.

E.8.2, E.8.3, F.6.3, and E.8.4 are `SOURCE REVIEWED`, not runtime accepted.
Final E.8 remains downstream. The deferred E.8.1 settings are not a blocker or
failure.

## Next Action

NEXT — user-controlled commit/push of the independently reviewed F.4.Refresh.2.Validation.1 farm materialization bundle/profile candidate.

This memory reconciliation stops for independent final pre-push review. Source
review does not establish ROOT/PyROOT, procedure-PDF, farm, runtime, production,
or Method-A promotion.

F.6.3 alone owns the private parallel full-analysis branch; E.8.4 only consumes
its post-yield sidecar in the existing pair-safe procedure rerender. Method B
remains numerically absent and Fit 1/Fit 2 dormant. Local checks and source
review do not establish ROOT/PyROOT, full-analysis, farm, or runtime validation.

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
