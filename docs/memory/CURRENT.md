---
memory_schema: 3
---
# Current KaonLT development state

## Active Objective

E.8 is `ACTIVE`: present the complete kaon missing-mass and signal-region
yield chain using authoritative objects and snapshots. Preserve baseline
production, frozen F.6.2 science, and detached Method-A/Method-B boundaries.

## Current Work Item

F.6.3 [current-baseline candidate-lineage adoption](phases/f6-3-current-baseline-candidate-lineage-adoption.md)
is `ACTIVE`: local implementation and deterministic checks are complete at
starting `test` HEAD `b349967c0d4210a78b144ce6134d3c1f15970245`, awaiting
independent ChatGPT actual-diff review. It selects the validated candidate
F.3/F.4 files with F.6.3-only identity records, reproduces F.4 using its exact
candidate-construction sentinel, and preserves live-cache parity before the
private `w0 -> w0*C` branch. No production promotion or runtime closure follows.

F.4.Refresh.2.Validation.2/Fix.1 remains `SOURCE REVIEWED` at pushed source
`712ba32b062772d87fe44efa6865346e5c438827`; operational source is unchanged
through the observed starting HEAD. Validation.1 and memory-health hardening
remain `SOURCE REVIEWED`.

## Verified State

F.4.Refresh.1 and Fix.1 are `CLOSED / RUNTIME VALIDATED` only for the detached
comparison gate. Farm comparison found exact F.2/F.3 scientific equality and
F.4 as the first changed stage. See [comparator evidence](evidence/f4-refresh1-current-baseline-authority-comparator-runtime-closure.md).

This adoption uses the detached candidate inputs independently accepted from
`KaonLT_F4_Refresh2_materialization_Q4p4W2p74_20261001-120157.zip`.
The [adoption phase](phases/f6-3-current-baseline-candidate-lineage-adoption.md)
records the supplied gate evidence and its narrow status: complete/error-free,
non-authoritative materialization, exact F.2/F.3 equality, F.4 first changed.
Historical accepted F.2/F.3/F.4 authorities and production remain unchanged.

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
- Accepted detached Refresh.2 ZIP SHA-256:
  `4cbdbdd2c0403b961614a8ccbb91104f13fef59194cd48e5e3570aded0954898`;
  farm-evaluated HEAD `b349967c0d4210a78b144ce6134d3c1f15970245`.
  Exact candidate identities and reconstruction-only zero-head boundary are in
  the [adoption phase](phases/f6-3-current-baseline-candidate-lineage-adoption.md).
- Active background profile is `no_empirical_residual`; legacy empirical
  residual Fit 1/Fit 2 are dormant, with both scales zero. See
  [durable knowledge](MEMORY.md).

## Blockers

The historical baseline-divergence failure is documented in
[direct evidence](evidence/e8-4-fix4-left-lowe-runtime-closure-and-f4-baseline-divergence.md).
Accepted Refresh.2 evidence now provides the current candidate lineage; the
local F.6.3 adoption repairs its historical-path/default-authority source
blocker. Independent actual-diff review, user commit/push, and pushed-state
synchronization remain required before the narrow debug runtime gate.
Fresh F.6.3/E.8.4 evidence remains pending; local tests do not establish
ROOT/PyROOT integration, full-analysis acceptance, or production promotion.

## Next Action

After independent actual-diff review, user commit/push, and pushed-state synchronization:

NEXT — run Q4p4W2p74 / Left / lowe through the tracked -d debug full-analysis path and inspect baseline-versus-reweighted missing-mass spectra and per-(t,phi) yield changes.

No farm, accepted-authority update, production promotion, or full-analysis
runtime claim follows from this source review or memory/status reconciliation.

## Success Criteria

Memory checkpoints bracket substantive work; batch nonblocking drift at
milestones. Only material active-state/provenance ambiguity or hard integrity
failure blocks progress. Required synchronization verifies the complete farm
execution chain without recursively creating reconciliation phases.

Left/lowe evidence must show baseline versus reweighted missing-mass spectra
and per-t comparisons, per-(t,phi) baseline and reweighted yields, absolute and
fractional yield changes, parent-t preservation, and the effect in the
procedure PDF. No canonical-five expansion or unrelated hardening/presentation
cleanup precedes that evidence unless a concrete blocker requires it.
E.8 consumes authoritative outputs; F.6.3 owns the private `w0 -> w0*C` branch,
and only F.6.4 could decide production promotion after runtime review.

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
- [F.6.3 candidate-lineage adoption](phases/f6-3-current-baseline-candidate-lineage-adoption.md)
