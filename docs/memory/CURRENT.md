---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

E.8 is `ACTIVE`: present the complete kaon missing-mass and signal-region
yield chain using authoritative objects and snapshots. Preserve baseline
production, frozen F.6.2 science, and detached Method-A/Method-B boundaries.

## Current Work Item

F.4.Refresh.2.Validation.2/Fix.1 materialize -> verify -> package owner is
`SOURCE REVIEWED` for the candidate reviewed against
committed `test` HEAD `c406c138285727503b12115177dfb8bc7efcb7fe`.
Independent ChatGPT actual-diff/source-runtime-path review of
`kaonlt_review(20261001-092421).diff` passed. The reviewed materializer,
wrapper/collector, source pin, and exact hardening blobs remain unchanged;
see the [Validation.2 phase](phases/f4-refresh2-validation2-tracked-execution-owner.md)
and [Fix.1 phase](phases/f4-refresh2-validation2-fix1-post-hardening-source-allowlist.md).
Validation.1 remains `SOURCE REVIEWED`; farm execution remains `BLOCKED` until
user commit/push and pushed-state review. No farm output or authority update exists.
Memory-health hardening and push-stable continuity remain `SOURCE REVIEWED`;
see their [hardening](phases/memory-health-operational-completeness-hardening.md)
and [continuity](phases/post-push-current-continuity.md) records.

## Verified State

F.4.Refresh.1 and Fix.1 are `CLOSED / RUNTIME VALIDATED` only for the detached
comparison gate. Farm comparison found exact F.2/F.3 scientific equality and
F.4 as the first changed stage. See [comparator evidence](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md).

F.4.Refresh.2 and Fix.1 are `SOURCE REVIEWED` for detached candidate
materialization. Validation.1 is `SOURCE REVIEWED` for the profile. Accepted
F.2/F.3/F.4 authorities remain unchanged; see the [materializer phase](phases/f4-refresh2-current-baseline-candidate-materialization.md)
and [profile phase](phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md).

F.6.2 and Fix.5 remain `CLOSED / RUNTIME VALIDATED` for accepted science and
presentation respectively; see [scientific evidence](evidence/f6-2-scientific-runtime-closure.md)
and [presentation evidence](evidence/f6-2-fix5-presentation-runtime-closure.md).
Historical F.1-F.6.2 closures remain in the [roadmap](roadmap/STATUS.md).

E.8.4.Fix.4 is `CLOSED / RUNTIME VALIDATED` only for stale-alignment-cache
rejection and recomputation; it does not close F.6.3/E.8.4 runtime. See
[Fix.4 evidence](evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md).
E.8.1.Fix.5/Fix.6 retain narrow `CLOSED / RUNTIME VALIDATED` Left/lowe
overlay and profile closure; canonical-five expansion remains `DEFERRED`.
See [Left/lowe evidence](evidence/e8-1-fix6-left-lowe-runtime-closure.md).

E.8.2, E.8.3, F.6.3, and E.8.4 remain `SOURCE REVIEWED`, without full-analysis
runtime acceptance. Final E.8 and F.6.4 remain `BLOCKED`; Method A stays
detached and Method B diagnostic-only. See [roadmap status](roadmap/STATUS.md).

## Source / Evidence Identity

- Frozen F.6.2 JSON SHA-256:
  `5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1`;
  artifact fingerprint:
  `ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0`;
  validation fingerprint:
  `7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b`.
  See [F.6.2 evidence](evidence/f6-2-scientific-runtime-closure.md).
- F.4.Refresh.2 materializer reviewed/pushed source:
  `141a3d04f9e5d07be21dba14e0e63212c3990bf1`; Validation.1 profile
  source at observed base `86590fa655512926f2e4d0c50bf12b57d5198da5`.
  See [profile phase](phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md).
- F.4.Refresh.1 comparator artifact SHA-256:
  `c16287e9192288ed5f116b5eebf95f49261d7fd09755ee923070154e047cf7b5`.
  See [comparator evidence](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md).
- Active background profile is `no_empirical_residual`; legacy empirical
  residual Fit 1/Fit 2 are dormant, with both scales zero. See
  [durable knowledge](MEMORY.md).

## Blockers

F.4.Refresh.2 farm execution remains `BLOCKED` pending user-controlled
commit/push and independent pushed-state review. Source review did not itself
establish commit, push, farm execution, or runtime validation. See
[investigation](investigations/f4-refresh2-validation1-operational-readiness-failure.md).

The fresh `Q4p4W2p74 / Left / lowe` F.6.3/E.8.4 runtime gate remains
`BLOCKED` by `f6_3_f4_shared_reproduction_failed:f3_fingerprint_input_content_mismatch`.
F.3 training (52,397) and application (55,380) projections matched exactly;
all 55,380 F.4/F.6.3 baseline projections differed in `analysis_MM`,
`baseline_pion_weight_w0`, and `signed_baseline_event_contribution`.
See [direct divergence evidence](evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md)
and [Refresh.1 comparator evidence](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md).

## Next Action

NEXT — after user-controlled commit/push and pushed-state review of the SOURCE REVIEWED F.4.Refresh.2.Validation.2/Fix.1 execution owner, prepare the single narrow Q4p4W2p74 F.4.Refresh.2 materialize -> verify -> package farm gate; do not begin F.6.3/E.8.4 until the returned F.4.Refresh.2 evidence is reviewed.

No farm, accepted-authority update, production promotion, or full-analysis
runtime claim follows from this source review or memory/status reconciliation.

## Success Criteria

The reviewed hardening requires strict memory health with zero warnings, a concise
CURRENT, durable operational-readiness failure record, and a pre-farm audit of
every tracked executable owner. A multi-step operation without a reviewed,
pushed driver remains blocked. E.8 presentation consumes authoritative inputs
without recomputing science; F.6.3 alone owns the private parallel `w0 -> w0*C`
branch, and F.6.4 alone could decide production promotion after runtime review.

## Do Not Reopen Without New Evidence

Do not reopen F.1-F.6.2 accepted science, change frozen F.6.2 identities,
alter the baseline production branch, make Method B numerical, or promote
Method A. Narrow Fix.4 cache and E.8.1 Left/lowe closures do not imply
F.6.3/E.8.4 or canonical-five runtime closure.

## Relevant References

- [Hardening phase](phases/memory-health-operational-completeness-hardening.md)
- [Operational-readiness investigation](investigations/f4-refresh2-validation1-operational-readiness-failure.md)
- [F.4.Refresh.2 materializer phase](phases/f4-refresh2-current-baseline-candidate-materialization.md)
- [F.4.Refresh.2.Validation.1 profile phase](phases/f4-refresh2-validation1-farm-materialization-bundle-profile.md)
- [F.4.Refresh.1 comparator evidence](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md)
- [F.6.2 scientific closure](evidence/f6-2-scientific-runtime-closure.md)
- [E.8 procedure roadmap](decisions/e8-full-analysis-procedure-roadmap.md)
- [Roadmap status](roadmap/STATUS.md)
